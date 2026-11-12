#!/usr/bin/env python3
"""Step 10: does the sense-strand deficit survive a composition-preserving shuffle?

Step 09b found fewer potential start sites on the sense strand of protein-coding
introns (11 / 20 / 32% for weak / medium / strong sites, chr8+chr9). Two kinds of
cause are possible:
  - composition: the sense strand of introns is T- and G-rich (a known effect of
    transcription on mutation), and the model may simply predict fewer starts on
    such sequence;
  - arrangement: specific start-like sequences (motifs) are missing on the sense
    strand, beyond what composition implies.
This step separates the two. Inside every eligible intronic segment (same
definitions as steps 09/09b) the sequence is shuffled so that the counts of
single bases AND of neighbouring base pairs (dinucleotides) are kept exactly,
while the order is randomised (Altschul-Erickson shuffle, by a random Eulerian
path). Keeping dinucleotides of the plus strand also keeps them on the minus
strand, so the composition of the sense strand of every gene is unchanged.
Everything else (exons with their margins, promoters, other genes, gaps) keeps
its real sequence, so the context of each segment is unchanged too.

Shuffling is done in chunks of at most CHUNK nt within runs of the same gene and
the same repeat status, so regional composition (e.g. GC level along a gene) is
also preserved. Intergenic eligible segments are shuffled the same way and remain
the orientation-free control.

The shuffled chromosomes are scored with Puffin-D exactly as in step 08 (100 kb
windows every 50 kb, central 50 kb kept, forward and reverse-complement averaged),
and sites are counted exactly as in step 09b (above-threshold positions merged
within MERGE nt), with the thresholds of the REAL predictions, so a site means
the same strength before and after shuffling.

Reading the result, for protein-coding introns:
  - deficit gone after shuffling -> it sits in the arrangement of the sequence
    (specific start-like motifs), not in composition;
  - deficit unchanged -> composition alone carries it;
  - in between -> both contribute; the remaining fraction is reported.

Checks printed first: shuffling kept base and dinucleotide counts exactly (per
chromosome); the real-track counts reproduce step 09b.

Usage:
  python3 10_shuffle_test.py --tracks /scratch/.../08_tracks --chroms chr8,chr9 \
      --genome /path/hg38.fa --gtf ../02_istart_baseline/ref/gencode.v47.annotation.gtf.gz \
      --puffin-dir /path/tools/puffin --seeds 1,2 --outdir /scratch/.../10_shuffle \
      --out-prefix results/shuffle_chr8_chr9 [--thresholds a,b,c] [--no-cuda]
"""
import argparse
import gzip
import hashlib
import os
import re
import sys
import time

import numpy as np
import pandas as pd

SPLICE_PAD = 100
TSS_PAD = 1000
INTERGENIC_PAD = 10000
BLOCK = 1_000_000
TAILS = (1e-3, 3e-4, 1e-4)
MERGE = 50
CHUNK = 500
CATEGORIES = ("protein_coding", "lncRNA")
WINDOW = 100000
STRIDE = 50000
FLANK = (WINDOW - STRIDE) // 2
TRACK_PRO_CAP_PLUS = 4
TRACK_PRO_CAP_MINUS = 5
COMP = str.maketrans("ACGTN", "TGCAN")

RE_ID = re.compile(r'gene_id "([^"]+)"')
RE_TYPE = re.compile(r'gene_type "([^"]+)"')
RE_NAME = re.compile(r'gene_name "([^"]+)"')


# ---------- annotation and masks (as in steps 09 / 09b) ----------
def parse_gtf(path, chroms):
    genes, exons, tss = [], {c: [] for c in chroms}, {c: [] for c in chroms}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in exons:
                continue
            feat = f[2]
            s, e, st = int(f[3]) - 1, int(f[4]), f[6]
            if feat == "exon":
                exons[f[0]].append((s, e))
            elif feat == "transcript":
                tss[f[0]].append(s if st == "+" else e - 1)
            elif feat == "gene":
                a = f[8]
                nm = RE_NAME.search(a)
                genes.append({"chrom": f[0], "start": s, "end": e, "strand": st,
                              "gene_id": RE_ID.search(a).group(1), "gene_type": RE_TYPE.search(a).group(1),
                              "gene_name": nm.group(1) if nm else ""})
    return pd.DataFrame(genes), exons, tss


def coverage(intervals, L, pad=0):
    d = np.zeros(L + 1, dtype=np.int32)
    for s, e in intervals:
        d[max(0, s - pad)] += 1
        d[min(L, e + pad)] -= 1
    return np.cumsum(d[:-1])


def build_masks(c, L, ok, genes, exons, tss):
    """Eligible positions per gene (list of (gene row, start, end, sel)) and intergenic mask."""
    excl = coverage(exons[c], L, pad=SPLICE_PAD) > 0
    excl |= coverage([(t, t + 1) for t in tss[c]], L, pad=TSS_PAD) > 0
    g = genes[genes["chrom"] == c]
    body = {st: coverage(list(zip(g.loc[g.strand == st, "start"], g.loc[g.strand == st, "end"])), L)
            for st in ("+", "-")}
    gc_ = g[g["gene_type"].isin(CATEGORIES)]
    cat_cov = {st: coverage(list(zip(gc_.loc[gc_.strand == st, "start"], gc_.loc[gc_.strand == st, "end"])), L)
               for st in ("+", "-")}
    any_pad = coverage(list(zip(g["start"], g["end"])), L, pad=INTERGENIC_PAD) > 0
    glist = []
    for r in gc_.itertuples(index=False):
        s, e = r.start, r.end
        opp = "-" if r.strand == "+" else "+"
        sel = ok[s:e] & ~excl[s:e] & (body[opp][s:e] == 0) & (cat_cov[r.strand][s:e] == 1)
        if sel.any():
            glist.append((r, s, e, sel))
    ig = ok & ~excl & ~any_pad
    return glist, ig


# ---------- dinucleotide-preserving shuffle ----------
def dinuc_shuffle(s, rng):
    """Altschul-Erickson shuffle of a uint8 sequence: keeps every dinucleotide count,
    the first and the last base (random Eulerian path with a random last-edge tree)."""
    n = len(s)
    if n < 3:
        return s.copy()
    a, b = s[:-1], s[1:]
    t = int(s[-1])
    succ = {int(v): b[a == v] for v in np.unique(a)}
    while True:
        last = {v: int(rng.integers(len(arr))) for v, arr in succ.items() if v != t}
        good = True
        for v in last:
            seen, u = set(), v
            while u != t:
                if u in seen or u not in last:
                    good = False
                    break
                seen.add(u)
                u = int(succ[u][last[u]])
            if not good:
                break
        if good:
            break
    lists = {}
    for v, arr in succ.items():
        if v in last:
            li = last[v]
            lists[v] = rng.permutation(np.delete(arr, li)).tolist() + [int(arr[li])]
        else:
            lists[v] = rng.permutation(arr).tolist()
    ptr = dict.fromkeys(lists, 0)
    cur = int(s[0])
    out = [cur]
    for _ in range(n - 1):
        i = ptr[cur]
        nxt = lists[cur][i]
        ptr[cur] = i + 1
        out.append(nxt)
        cur = nxt
    return np.asarray(out, dtype=np.uint8)


def shuffle_chrom(up, is_rep, glist, ig, rng):
    """Shuffle eligible positions in runs of constant (gene, repeat status), in chunks."""
    L = len(up)
    label = np.full(L, -2, dtype=np.int64)
    label[ig] = -1
    for gi, (_, s, e, sel) in enumerate(glist):
        seg = label[s:e]
        seg[sel] = gi
    key = np.where(label >= -1, label * 2 + is_rep.astype(np.int64), -10)
    out = up.copy()
    change = np.flatnonzero(np.diff(key)) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, L]
    for s0, e0 in zip(starts, ends):
        if key[s0] == -10:
            continue
        for c0 in range(s0, e0, CHUNK):
            c1 = min(e0, c0 + CHUNK)
            if c1 - c0 >= 3:
                out[c0:c1] = dinuc_shuffle(up[c0:c1], rng)
    return out, key != -10


def composition(x):
    mono = np.bincount(x, minlength=256)
    di = np.bincount(x[:-1].astype(np.int64) * 256 + x[1:], minlength=65536)
    return mono, di


# ---------- prediction (as in step 08) ----------
def window_plan(clen):
    if clen < WINDOW:
        return []
    starts = list(range(0, clen - WINDOW + 1, STRIDE))
    if starts[-1] != clen - WINDOW:
        starts.append(clen - WINDOW)
    plan, prev_hi = [], 0
    for k, s in enumerate(starts):
        hi = clen if k == len(starts) - 1 else s + FLANK + STRIDE
        plan.append((s, prev_hi, hi))
        prev_hi = hi
    return plan


def tracks_of(arr):
    arr = np.asarray(arr, dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[0]
    return arr[TRACK_PRO_CAP_PLUS], arr[TRACK_PRO_CAP_MINUS]


def predict_chrom(model, seqstr, max_n=0.5):
    L = len(seqstr)
    P = np.full(L, np.nan, dtype=np.float32)
    M = np.full(L, np.nan, dtype=np.float32)
    for s, lo, hi in window_plan(L):
        seq = seqstr[s:s + WINDOW]
        if seq.count("N") / WINDOW > max_n:
            continue
        pf, mf = tracks_of(model.predict(seq))
        rp, rm = tracks_of(model.predict(seq.translate(COMP)[::-1]))
        P[lo:hi] = ((pf + rm[::-1]) / 2)[lo - s:hi - s]
        M[lo:hi] = ((mf + rp[::-1]) / 2)[lo - s:hi - s]
    return P, M


# ---------- counting (as in step 09b) ----------
def n_sites(hit):
    idx = np.flatnonzero(hit)
    return 0 if len(idx) == 0 else 1 + int((np.diff(idx) > MERGE).sum())


def count_into(row, tag, ts, ta, sel, rep, thr, tkeys):
    for cls, cm in (("unique", sel & ~rep), ("all", sel)):
        for t, k in zip(TAILS, tkeys):
            row[f"{tag}_{cls}_{k}_sense"] = n_sites(cm & (ts >= thr[t]))
            row[f"{tag}_{cls}_{k}_anti"] = n_sites(cm & (ta >= thr[t]))
        row[f"{tag}_{cls}_sum_sense"] = float(ts[cm].sum())
        row[f"{tag}_{cls}_sum_anti"] = float(ta[cm].sum())


def paired_boot(So, Ao, Ss, As, n_boot, rng):
    """Deficits (%) before and after shuffling, and their difference, bootstrapped over units."""
    So, Ao, Ss, As = (np.asarray(x, float) for x in (So, Ao, Ss, As))
    n = len(So)
    idx = rng.integers(0, n, size=(n_boot, n))
    ao, as_ = Ao[idx].sum(1), As[idx].sum(1)
    good = (ao > 0) & (as_ > 0)
    do = 100 * (1 - So[idx].sum(1)[good] / ao[good])
    ds = 100 * (1 - Ss[idx].sum(1)[good] / as_[good])
    pct = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    return pct(do), pct(ds), pct(do - ds)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks", required=True, help="real tracks, output directory of step 08")
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--puffin-dir", required=True)
    ap.add_argument("--seeds", default="1,2")
    ap.add_argument("--outdir", required=True, help="where shuffled tracks are written (scratch)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--thresholds", default="", help="override: comma-separated values for top 1e-3,3e-4,1e-4")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--no-cuda", action="store_true")
    a = ap.parse_args()
    chroms = a.chroms.split(",")
    seeds = [int(x) for x in a.seeds.split(",")]
    tracks_dir, genome, gtf = (os.path.abspath(x) for x in (a.tracks, a.genome, a.gtf))
    outdir, prefix, puffin_dir = (os.path.abspath(x) for x in (a.outdir, a.out_prefix, a.puffin_dir))
    os.makedirs(outdir, exist_ok=True)

    from pyfaidx import Fasta
    fa = Fasta(genome, as_raw=True, sequence_always_upper=False)
    print(f"reading {gtf} ...", flush=True)
    genes, exons, tss = parse_gtf(gtf, set(chroms))

    real = {}
    for c in chroms:
        p = np.load(os.path.join(tracks_dir, f"{c}.procap_plus.f16.npy")).astype(np.float32)
        m = np.load(os.path.join(tracks_dir, f"{c}.procap_minus.f16.npy")).astype(np.float32)
        if len(p) != len(fa[c]):
            sys.exit(f"{c}: track length {len(p)} != genome length {len(fa[c])}")
        real[c] = (p, m)
    tkeys = [f"top{t:g}" for t in TAILS]
    if a.thresholds:
        thr = dict(zip(TAILS, (float(x) for x in a.thresholds.split(","))))
    else:
        allv = np.concatenate([np.concatenate([p[np.isfinite(p)], m[np.isfinite(m)]]) for p, m in real.values()])
        thr = {t: float(np.quantile(allv, 1 - t)) for t in TAILS}
        del allv
    print("thresholds (from the real predictions): " + ", ".join(f"{t:g} -> {v:.4f}" for t, v in thr.items()))

    cwd = os.getcwd()
    os.chdir(puffin_dir)
    sys.path.insert(0, puffin_dir)
    from puffin_D import Puffin_D
    model = Puffin_D(use_cuda=not a.no_cuda)
    os.chdir(cwd)

    gene_rows, block_rows = [], []
    for c in chroms:
        t0 = time.time()
        p, m = real[c]
        L = len(p)
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        is_rep = raw >= 97
        up = np.where(is_rep, raw - 32, raw).astype(np.uint8)
        del raw
        ok = np.isfinite(p) & np.isfinite(m) & (up != ord("N"))
        glist, ig = build_masks(c, L, ok, genes, exons, tss)
        grows = [{"chrom": c, "gene_id": r.gene_id, "gene_name": r.gene_name, "gene_type": r.gene_type,
                  "strand": r.strand, "start": s, "end": e, "eligible_nt": int(sel.sum())}
                 for r, s, e, sel in glist]
        bspans = [(b0, min(L, b0 + BLOCK)) for b0 in range(0, L, BLOCK) if ig[b0:min(L, b0 + BLOCK)].any()]
        brows = [{"chrom": c, "block_start": b0, "eligible_nt": int(ig[b0:b1].sum())} for b0, b1 in bspans]

        def count_all(tag, P, Mn):
            for row, (r, s, e, sel) in zip(grows, glist):
                ts, ta = (P[s:e], Mn[s:e]) if r.strand == "+" else (Mn[s:e], P[s:e])
                count_into(row, tag, ts, ta, sel, is_rep[s:e], thr, tkeys)
            for row, (b0, b1) in zip(brows, bspans):
                count_into(row, tag, P[b0:b1], Mn[b0:b1], ig[b0:b1], is_rep[b0:b1], thr, tkeys)

        count_all("real", p, m)
        print(f"[{c}] masks and real counts ready ({time.time() - t0:.0f} s); "
              f"{len(glist)} genes, intergenic {int(ig.sum()):,} nt", flush=True)

        for seed in seeds:
            t1 = time.time()
            rng = np.random.default_rng([seed, chroms.index(c)])
            shuf, shuffled_mask = shuffle_chrom(up, is_rep, glist, ig, rng)
            m0, d0 = composition(up)
            m1, d1 = composition(shuf)
            same = bool(np.array_equal(m0, m1) and np.array_equal(d0, d1))
            changed = float((shuf != up)[shuffled_mask].mean())
            md5 = hashlib.md5(shuf.tobytes()).hexdigest()
            print(f"  [{c} seed {seed}] shuffled {int(shuffled_mask.sum()):,} nt in {time.time() - t1:.0f} s; "
                  f"base and dinucleotide counts identical: {same}; positions changed: {100 * changed:.1f}%; "
                  f"md5 {md5}", flush=True)
            if not same:
                sys.exit("composition check failed")
            t2 = time.time()
            os.chdir(puffin_dir)
            P, Mn = predict_chrom(model, shuf.tobytes().decode("ascii"))
            os.chdir(cwd)
            sd = os.path.join(outdir, f"seed{seed}")
            os.makedirs(sd, exist_ok=True)
            np.save(os.path.join(sd, f"{c}.procap_plus.f16.npy"), P.astype(np.float16))
            np.save(os.path.join(sd, f"{c}.procap_minus.f16.npy"), Mn.astype(np.float16))
            print(f"  [{c} seed {seed}] predicted in {time.time() - t2:.0f} s", flush=True)
            count_all(f"shuf{seed}", P, Mn)
            del P, Mn, shuf
        gene_rows += grows
        block_rows += brows
        del up, is_rep, ok, glist, ig

    G = pd.DataFrame(gene_rows)
    B = pd.DataFrame(block_rows)
    G.to_csv(prefix + "_genes.tsv", sep="\t", index=False)
    B.to_csv(prefix + "_blocks.tsv", sep="\t", index=False)

    rng = np.random.default_rng(1)
    rows = []
    sets = [("intron_protein_coding", G[G["gene_type"] == "protein_coding"], "gene"),
            ("intron_lncRNA", G[G["gene_type"] == "lncRNA"], "gene"),
            ("intergenic_plus_vs_minus", B, "1Mb_block")]
    for label, df, unit in sets:
        for cls in ("unique", "all"):
            for k in tkeys + ["sum"]:
                So, Ao = df[f"real_{cls}_{k}_sense"], df[f"real_{cls}_{k}_anti"]
                for seed in seeds:
                    Ss, As = df[f"shuf{seed}_{cls}_{k}_sense"], df[f"shuf{seed}_{cls}_{k}_anti"]
                    if min(Ao.sum(), As.sum()) == 0:
                        continue
                    dr = 100 * (1 - So.sum() / Ao.sum())
                    dsh = 100 * (1 - Ss.sum() / As.sum())
                    (rlo, rhi), (slo, shi), (dlo, dhi) = paired_boot(So, Ao, Ss, As, a.n_boot, rng)
                    rows.append({"set": label, "unit": unit, "n_units": len(df), "seq": cls, "metric": k,
                                 "seed": seed, "real_sense": So.sum(), "real_anti": Ao.sum(),
                                 "shuf_sense": Ss.sum(), "shuf_anti": As.sum(),
                                 "real_deficit_pct": dr, "real_lo95": rlo, "real_hi95": rhi,
                                 "shuf_deficit_pct": dsh, "shuf_lo95": slo, "shuf_hi95": shi,
                                 "drop_pct_points": dr - dsh, "drop_lo95": dlo, "drop_hi95": dhi,
                                 "remaining_fraction": dsh / dr if dr != 0 else np.nan,
                                 "sites_shuf_vs_real": (Ss.sum() + As.sum()) / (So.sum() + Ao.sum())})
    R = pd.DataFrame(rows)
    R["test"] = ("deficit = 100*(1 - total sense / total antisense) over units; real and shuffled bootstrapped "
                 "together over units (paired); drop = real - shuffled")
    R.to_csv(prefix + "_summary.tsv", sep="\t", index=False)

    fmt = lambda x: f"{x:.3g}"
    cols = ["set", "seq", "metric", "seed", "real_sense", "real_anti", "real_deficit_pct", "shuf_sense",
            "shuf_anti", "shuf_deficit_pct", "shuf_lo95", "shuf_hi95", "drop_pct_points", "drop_lo95",
            "drop_hi95", "sites_shuf_vs_real"]
    print("\n== sites: deficit on the sense strand, real sequence vs composition-preserving shuffle ==")
    print("   (intergenic: deficit of plus vs minus; should be ~0 in both)")
    print(R[R["metric"] != "sum"][cols].to_string(index=False, float_format=fmt))
    print("\n== sum of predictions (low background; the 'more background on the sense strand' puzzle) ==")
    print(R[R["metric"] == "sum"][cols].to_string(index=False, float_format=fmt))
    print(f"\nCIs: percentile bootstrap over units, {a.n_boot} resamples; sites merged within {MERGE} nt; "
          f"shuffle chunks <= {CHUNK} nt.")
    print(f"written: {prefix}_genes.tsv, _blocks.tsv, _summary.tsv; shuffled tracks in {outdir}")


if __name__ == "__main__":
    main()
