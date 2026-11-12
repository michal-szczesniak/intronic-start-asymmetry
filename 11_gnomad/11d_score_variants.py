#!/usr/bin/env python3
"""Step 11d: score gnomAD SNVs for their effect on predicted initiation (both strands).

For every PASS SNV at an eligible position (definitions of steps 09/09b/10: protein-
coding and lncRNA introns, intergenic control), the predicted PRO-cap signal with
the reference and with the alternative allele is read within +/-HW nt of the
variant, on the plus and the minus strand (forward and reverse-complement averaged,
as in step 08). Two stages, because step 11c showed that the model has weak
long-range cross-talk:

  Stage 1, screen: many variants share one 100 kb window, spaced >= SPACING nt apart
    (bfloat16, fast). Error vs single-variant runs at 1 kb spacing was small
    (99th percentile ~2e-3), so this is used only to decide which variants need
    an exact value.
  Stage 2, exact: every variant whose screened change exceeds FLAG on either strand,
    or whose readout comes within NEAR of a site threshold, is scored alone
    (float32), with the reference computed in the same batch shape. A random
    AUDIT fraction of the other variants is also scored alone, to measure what the
    screen misses.

Windows are tiled as in step 08 (100 kb every 50 kb, central part kept), so each
variant is scored once, with >= 25 kb of context on both sides.

Variants scored: all protein-coding intronic SNVs, and random subsets of lncRNA
intronic (--lnc-frac) and intergenic (--ig-frac) SNVs. The reference allele is
checked against the genome.

Output (scratch), one gzipped table per chromosome, one row per variant:
  chrom pos ref alt AC AN set gene_id strand ref_plus ref_minus scr_alt_plus
  scr_alt_minus exact alt_plus alt_minus
  - strand: gene strand ('+' for intergenic, where plus is arbitrary);
  - exact: 0 = alt_* are screen values, 1 = flagged and scored alone,
           2 = audit sample scored alone;
  The job can be restarted: finished windows are listed in <out>.done and skipped.

Usage:
  python3 11d_score_variants.py --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --tracks /scratch/.../08_tracks --gnomad /scratch/.../gnomad_v4.1 --chrom chr8 \
      --puffin-dir /path/tools/puffin --out /scratch/.../11d_scores/chr8.scores.tsv.gz \
      [--limit-windows 20] [--batch 16] [--no-cuda]
"""
import argparse
import gzip
import os
import re
import sys
import time

import numpy as np
import pandas as pd

SPLICE_PAD = 100
TSS_PAD = 1000
INTERGENIC_PAD = 10000
CATEGORIES = ("protein_coding", "lncRNA")
WINDOW = 100000
STRIDE = 50000
FLANK = (WINDOW - STRIDE) // 2
HW = 50
TRACK_PRO_CAP_PLUS = 4
TRACK_PRO_CAP_MINUS = 5
SPACING = 1000
THRESHOLDS = (0.0439, 0.1432)
FLAG = 0.003
NEAR = 0.003
COMP = str.maketrans("ACGTN", "TGCAN")
LUT = np.full((256, 4), 0.25, dtype=np.float32)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = 0
    LUT[ord(_b), _i] = 1
    LUT[ord(_b.lower())] = LUT[ord(_b)]
OFFS = np.arange(-HW, HW + 1)
RE_ID = re.compile(r'gene_id "([^"]+)"')
RE_TYPE = re.compile(r'gene_type "([^"]+)"')


# ---------- annotation and eligible positions (as in steps 09b / 10) ----------
def parse_gtf(path, chrom):
    genes, exons, tss = [], [], []
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom:
                continue
            s, e, st = int(f[3]) - 1, int(f[4]), f[6]
            if f[2] == "exon":
                exons.append((s, e))
            elif f[2] == "transcript":
                tss.append(s if st == "+" else e - 1)
            elif f[2] == "gene":
                genes.append({"start": s, "end": e, "strand": st, "gene_id": RE_ID.search(f[8]).group(1),
                              "gene_type": RE_TYPE.search(f[8]).group(1)})
    return pd.DataFrame(genes), exons, tss


def coverage(intervals, L, pad=0):
    d = np.zeros(L + 1, dtype=np.int32)
    for s, e in intervals:
        d[max(0, s - pad)] += 1
        d[min(L, e + pad)] -= 1
    return np.cumsum(d[:-1])


def eligible_labels(L, ok, g, exons, tss):
    """label[i] = index into gl (gene list) for eligible intronic positions, -1 for
    eligible intergenic positions, -2 otherwise."""
    excl = coverage(exons, L, pad=SPLICE_PAD) > 0
    excl |= coverage([(t, t + 1) for t in tss], L, pad=TSS_PAD) > 0
    body = {st: coverage(list(zip(g.loc[g.strand == st, "start"], g.loc[g.strand == st, "end"])), L)
            for st in ("+", "-")}
    gc_ = g[g["gene_type"].isin(CATEGORIES)].reset_index(drop=True)
    cat_cov = {st: coverage(list(zip(gc_.loc[gc_.strand == st, "start"], gc_.loc[gc_.strand == st, "end"])), L)
               for st in ("+", "-")}
    any_pad = coverage(list(zip(g["start"], g["end"])), L, pad=INTERGENIC_PAD) > 0
    label = np.full(L, -2, dtype=np.int32)
    label[ok & ~excl & ~any_pad] = -1
    for gi, r in enumerate(gc_.itertuples(index=False)):
        s, e = r.start, r.end
        opp = "-" if r.strand == "+" else "+"
        sel = ok[s:e] & ~excl[s:e] & (body[opp][s:e] == 0) & (cat_cov[r.strand][s:e] == 1)
        seg = label[s:e]
        seg[sel] = gi
    return label, gc_


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


# ---------- model ----------
class Scorer:
    def __init__(self, model, use_cuda, batch):
        self.model, self.cuda, self.B = model, use_cuda, batch
        self.torch = None
        try:
            import torch
            if isinstance(model, torch.nn.Module):
                self.torch = torch
                torch.backends.cudnn.benchmark = False
                torch.backends.cudnn.deterministic = True
        except ImportError:
            pass

    def raw(self, seqs, bf16):
        """(n, 2, L) array of PRO-cap plus/minus tracks for each input sequence as given.
        Batches are always padded to the same size, so every pass uses the same kernels."""
        n, L = len(seqs), len(seqs[0])
        out = np.empty((n, 2, L), dtype=np.float32)
        if self.torch is None:
            for k, s in enumerate(seqs):
                a = np.asarray(self.model.predict(s), dtype=np.float32).reshape(-1, L)
                out[k] = a[[TRACK_PRO_CAP_PLUS, TRACK_PRO_CAP_MINUS]]
            return out
        torch = self.torch
        for i in range(0, n, self.B):
            chunk = seqs[i:i + self.B]
            m = len(chunk)
            chunk = chunk + [chunk[-1]] * (self.B - m)
            x = torch.from_numpy(np.stack([LUT[np.frombuffer(s.encode(), dtype=np.uint8)].T for s in chunk]))
            if self.cuda:
                x = x.cuda()
            with torch.no_grad():
                if bf16:
                    with torch.autocast("cuda", dtype=torch.bfloat16):
                        y = self.model(x)
                else:
                    y = self.model(x)
            out[i:i + m] = y[:m, [TRACK_PRO_CAP_PLUS, TRACK_PRO_CAP_MINUS]].float().cpu().numpy()
        return out

    def readouts(self, seqs, pos_lists, bf16):
        """For each sequence k and its positions pos_lists[k]: max within +/-HW of the
        forward+RC averaged plus and minus tracks. Returns list of (n_k, 2) arrays."""
        n = len(seqs)
        r = self.raw(list(seqs) + [s.translate(COMP)[::-1] for s in seqs], bf16)
        L = len(seqs[0])
        res = []
        for k in range(n):
            idx = np.asarray(pos_lists[k])[:, None] + OFFS[None, :]
            plus = (r[k, 0][idx] + r[n + k, 1][L - 1 - idx]) / 2
            minus = (r[k, 1][idx] + r[n + k, 0][L - 1 - idx]) / 2
            res.append(np.stack([plus.max(axis=1), minus.max(axis=1)], axis=1))
        return res


def chunked(S, seqs, poslists, bf16, chunk=32):
    """S.readouts over chunks of sequences, to bound memory."""
    out = []
    for i in range(0, len(seqs), chunk):
        out.extend(S.readouts(seqs[i:i + chunk], poslists[i:i + chunk], bf16))
    return out


def mutate(seq, offs, alts):
    b = bytearray(seq.encode())
    for o, a in zip(offs, alts):
        b[o] = ord(a)
    return b.decode()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--tracks", required=True, help="step 08 tracks (only used for the eligible-position mask)")
    ap.add_argument("--gnomad", required=True, help="directory with <chrom>.snv_pass.tsv.gz from step 11b")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--puffin-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lnc-frac", type=float, default=0.3)
    ap.add_argument("--ig-frac", type=float, default=0.2)
    ap.add_argument("--audit-frac", type=float, default=0.005)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--limit-windows", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-cuda", action="store_true")
    a = ap.parse_args()
    c = a.chrom
    genome, gtf, tracks, gdir, out, puffin_dir = (os.path.abspath(x) for x in
                                                  (a.genome, a.gtf, a.tracks, a.gnomad, a.out, a.puffin_dir))
    rng = np.random.default_rng([a.seed, int(re.sub(r"\D", "", c) or 0)])
    t0 = time.time()

    from pyfaidx import Fasta
    fa = Fasta(genome, as_raw=True, sequence_always_upper=False)
    L = len(fa[c])
    up = np.frombuffer(str(fa[c][:]).upper().encode(), dtype=np.uint8)
    p8 = np.load(os.path.join(tracks, f"{c}.procap_plus.f16.npy"), mmap_mode="r")
    m8 = np.load(os.path.join(tracks, f"{c}.procap_minus.f16.npy"), mmap_mode="r")
    ok = np.isfinite(p8) & np.isfinite(m8) & (up != ord("N"))
    g, exons, tss = parse_gtf(gtf, c)
    label, gl = eligible_labels(L, ok, g, exons, tss)
    del ok, p8, m8
    print(f"[{c}] eligible positions ready in {time.time() - t0:.0f} s", flush=True)

    v = pd.read_csv(os.path.join(gdir, f"{c}.snv_pass.tsv.gz"), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "alt", "AC", "AN"],
                    dtype={"chrom": "category", "pos": np.int64, "ref": "category", "alt": "category",
                           "AC": np.int64, "AN": np.int64})
    n0 = len(v)
    i0 = v["pos"].to_numpy() - 1
    refok = up[i0] == np.frombuffer("".join(v["ref"].astype(str)).encode(), dtype=np.uint8)
    lab = label[i0]
    keep_p = rng.random(len(v))
    is_pc = (lab >= 0) & (gl["gene_type"].to_numpy()[np.maximum(lab, 0)] == "protein_coding")
    is_lnc = (lab >= 0) & ~is_pc
    is_ig = lab == -1
    sel = refok & (is_pc | (is_lnc & (keep_p < a.lnc_frac)) | (is_ig & (keep_p < a.ig_frac)))
    v = v[sel].reset_index(drop=True)
    lab = lab[sel]
    v["set"] = np.where(lab == -1, "intergenic", np.where(is_pc[sel], "intron_protein_coding", "intron_lncRNA"))
    v["gene_id"] = np.where(lab >= 0, gl["gene_id"].to_numpy()[np.maximum(lab, 0)], "")
    v["strand"] = np.where(lab >= 0, gl["strand"].to_numpy()[np.maximum(lab, 0)], "+")
    v["audit"] = rng.random(len(v)) < a.audit_frac
    vpos = v["pos"].to_numpy() - 1
    print(f"[{c}] {n0:,} SNVs; reference allele matches the genome for {100 * refok.mean():.3f}%; "
          f"to score: {len(v):,} (" + ", ".join(f"{k} {n:,}" for k, n in v["set"].value_counts().items()) + ")",
          flush=True)

    os.chdir(puffin_dir)
    sys.path.insert(0, puffin_dir)
    from puffin_D import Puffin_D
    model = Puffin_D(use_cuda=not a.no_cuda)
    S = Scorer(model, not a.no_cuda, a.batch)
    print(f"[{c}] batched torch path: {S.torch is not None}, batch {a.batch}", flush=True)

    done_path = out + ".done"
    done = set()
    if os.path.exists(done_path):
        done = {int(x) for x in open(done_path).read().split()}
    plan = window_plan(L)
    if a.limit_windows:
        plan = plan[:a.limit_windows]
    cols = ["chrom", "pos", "ref", "alt", "AC", "AN", "set", "gene_id", "strand", "ref_plus", "ref_minus",
            "scr_alt_plus", "scr_alt_minus", "exact", "alt_plus", "alt_minus"]
    if not os.path.exists(out):
        with gzip.open(out, "wt") as fh:
            fh.write("\t".join(cols) + "\n")
    print(f"[{c}] {len(plan)} windows, {len(done)} already done", flush=True)

    t1 = time.time()
    n_scored = n_flag = n_aud = 0
    for wi, (s, lo, hi) in enumerate(plan):
        if s in done:
            continue
        a_, b_ = np.searchsorted(vpos, [lo, hi])
        if b_ == a_:
            with open(done_path, "a") as fh:
                fh.write(f"{s}\n")
            continue
        offs = vpos[a_:b_] - s
        inside = (offs >= HW) & (offs < WINDOW - HW)   # readout must fit in the window
        w = v.iloc[a_:b_][inside]
        offs = offs[inside]
        if len(offs) == 0:
            with open(done_path, "a") as fh:
                fh.write(f"{s}\n")
            continue
        alts = w["alt"].astype(str).to_numpy()
        seq = str(fa[c][s:s + WINDOW]).upper()
        n = len(offs)
        # stage 1: screen, variants spread over K sets so that members of a set are >= SPACING apart
        K = int((np.searchsorted(offs, offs + SPACING) - np.arange(n)).max())
        sets = [np.arange(j, n, K) for j in range(K)]
        ref_bf = S.readouts([seq], [offs], bf16=True)[0]
        scr = np.empty((n, 2), dtype=np.float32)
        seqs = [mutate(seq, offs[ix], alts[ix]) for ix in sets]
        for ix, r in zip(sets, chunked(S, seqs, [offs[ix] for ix in sets], bf16=True)):
            scr[ix] = r
        del seqs
        d = np.abs(scr - ref_bf)
        flag = (d > FLAG).any(axis=1)
        for T in THRESHOLDS:
            flag |= ((np.minimum(scr, ref_bf) < T + NEAR) & (np.maximum(scr, ref_bf) > T - NEAR)).any(axis=1)
        aud = ~flag & w["audit"].to_numpy()
        # stage 2: exact, one variant per sequence, float32, reference in the same batch shape
        ref32 = S.readouts([seq], [offs], bf16=False)[0]
        alt = scr.copy()
        exact = np.zeros(n, dtype=np.int8)
        todo = np.flatnonzero(flag | aud)
        if len(todo):
            ex = chunked(S, [mutate(seq, [offs[i]], [alts[i]]) for i in todo], [[offs[i]] for i in todo], bf16=False)
            for i, r in zip(todo, ex):
                alt[i] = r[0]
            exact[flag] = 1
            exact[aud] = 2
        res = pd.DataFrame({
            "chrom": c, "pos": w["pos"].to_numpy(), "ref": w["ref"].astype(str).to_numpy(), "alt": alts,
            "AC": w["AC"].to_numpy(), "AN": w["AN"].to_numpy(), "set": w["set"].to_numpy(),
            "gene_id": w["gene_id"].to_numpy(), "strand": w["strand"].to_numpy(),
            "ref_plus": ref32[:, 0], "ref_minus": ref32[:, 1],
            "scr_alt_plus": scr[:, 0] - ref_bf[:, 0] + ref32[:, 0],
            "scr_alt_minus": scr[:, 1] - ref_bf[:, 1] + ref32[:, 1],
            "exact": exact, "alt_plus": np.where(exact > 0, alt[:, 0], scr[:, 0] - ref_bf[:, 0] + ref32[:, 0]),
            "alt_minus": np.where(exact > 0, alt[:, 1], scr[:, 1] - ref_bf[:, 1] + ref32[:, 1])})
        with gzip.open(out, "at") as fh:
            res.to_csv(fh, sep="\t", header=False, index=False, float_format="%.6g")
        with open(done_path, "a") as fh:
            fh.write(f"{s}\n")
        n_scored += n
        n_flag += int(flag.sum())
        n_aud += int(aud.sum())
        if (wi + 1) % 20 == 0 or wi == len(plan) - 1:
            el = time.time() - t1
            print(f"  window {wi + 1}/{len(plan)}: {n_scored:,} variants, flagged {100 * n_flag / max(1, n_scored):.2f}%, "
                  f"audit {n_aud:,}; {el / 60:.1f} min, {1000 * el / max(1, n_scored):.2f} ms per variant", flush=True)
    print(f"[{c}] finished in {(time.time() - t0) / 60:.1f} min; written {out}")


if __name__ == "__main__":
    main()
