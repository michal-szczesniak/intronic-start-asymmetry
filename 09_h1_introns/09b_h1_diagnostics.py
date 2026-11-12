#!/usr/bin/env python3
"""Step 09b: is the intronic sense-strand deficit of step 09 robust?

Step 09 found ~18% fewer weak potential start sites (top 1e-3 of predictions) on
the sense strand of protein-coding introns, but the intergenic control, which has
no orientation, was itself off by -31% (plus strand richer). This step asks four
questions with the same position definitions as step 09 (code copied, not changed):

  A. Sites instead of positions. Step 09 counted nucleotides above a threshold,
     so one broad peak or a run inside a simple repeat counted many times.
     Here consecutive above-threshold positions closer than MERGE nt are merged
     into one site, and sites are counted.
  B. Concentration. Is a result carried by a handful of genes or 1 Mb blocks?
     Reported: share of all sites in the top 1% of units, and the deficit after
     dropping the DROP_TOP units with the largest |sense - antisense|.
  C. Gene strand. Protein-coding introns split into plus- and minus-strand genes.
     A chromosome-scale excess on one strand (as the intergenic control hints)
     would push the two halves in opposite directions; a real sense effect pushes
     both the same way.
  D. Composition. Genes split into tertiles of sense-strand AT skew and GC skew.
     If the deficit grows with the skew, it travels with composition, which is
     what transcription-associated mutation produces (the rival explanation).

Thresholds: top 1e-3, 3e-4 and 1e-4 of all predictions on the chosen chromosomes
(both strands pooled), computed exactly as in step 09.

Outputs: per-gene table, per-block table (intergenic, 1 Mb) and a summary table.

Usage:
  python3 09b_h1_diagnostics.py --tracks /scratch/.../08_tracks --chroms chr8,chr9 \
      --genome /path/hg38.fa --gtf ../02_istart_baseline/ref/gencode.v47.annotation.gtf.gz \
      --out-prefix results/h1diag_chr8_chr9
"""
import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

SPLICE_PAD = 100
TSS_PAD = 1000
INTERGENIC_PAD = 10000
BLOCK = 1_000_000
TAILS = (1e-3, 3e-4, 1e-4)
MERGE = 50
DROP_TOP = 10
CATEGORIES = ("protein_coding", "lncRNA")

RE_ID = re.compile(r'gene_id "([^"]+)"')
RE_TYPE = re.compile(r'gene_type "([^"]+)"')
RE_NAME = re.compile(r'gene_name "([^"]+)"')


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


def n_sites(hit):
    """Number of sites: above-threshold positions merged when closer than MERGE nt."""
    idx = np.flatnonzero(hit)
    if len(idx) == 0:
        return 0
    return 1 + int((np.diff(idx) > MERGE).sum())


def ratio_ci(S, A, n_boot, rng):
    S, A = np.asarray(S, float), np.asarray(A, float)
    if len(S) < 2 or A.sum() == 0:
        return np.nan, np.nan, np.nan
    R = S.sum() / A.sum()
    idx = rng.integers(0, len(S), size=(n_boot, len(S)))
    a = A[idx].sum(axis=1)
    r = S[idx].sum(axis=1)[a > 0] / a[a > 0]
    return R, float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    chroms = a.chroms.split(",")

    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)
    print(f"reading {a.gtf} ...", flush=True)
    genes, exons, tss = parse_gtf(a.gtf, set(chroms))

    tracks = {}
    for c in chroms:
        p = np.load(os.path.join(a.tracks, f"{c}.procap_plus.f16.npy")).astype(np.float32)
        m = np.load(os.path.join(a.tracks, f"{c}.procap_minus.f16.npy")).astype(np.float32)
        if len(p) != len(fa[c]):
            sys.exit(f"{c}: track length {len(p)} != genome length {len(fa[c])}")
        tracks[c] = (p, m)
    allv = np.concatenate([np.concatenate([p[np.isfinite(p)], m[np.isfinite(m)]]) for p, m in tracks.values()])
    thr = {t: float(np.quantile(allv, 1 - t)) for t in TAILS}
    del allv
    tkeys = [f"top{t:g}" for t in TAILS]
    print("thresholds: " + ", ".join(f"{t:g} -> {v:.4f}" for t, v in thr.items()), flush=True)

    def count(row, cls, ts, ta, cm):
        row[f"{cls}_nt"] = int(cm.sum())
        for t, k in zip(TAILS, tkeys):
            hs, ha = cm & (ts >= thr[t]), cm & (ta >= thr[t])
            row[f"{cls}_{k}_pos_sense"], row[f"{cls}_{k}_pos_anti"] = int(hs.sum()), int(ha.sum())
            row[f"{cls}_{k}_site_sense"], row[f"{cls}_{k}_site_anti"] = n_sites(hs), n_sites(ha)

    gene_rows, block_rows = [], []
    for c in chroms:
        p, m = tracks[c]
        L = len(p)
        seq = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        is_rep = seq >= 97
        up = np.where(is_rep, seq - 32, seq)
        ok = np.isfinite(p) & np.isfinite(m) & (up != ord("N"))
        excl = coverage(exons[c], L, pad=SPLICE_PAD) > 0
        excl |= coverage([(t, t + 1) for t in tss[c]], L, pad=TSS_PAD) > 0
        g = genes[genes["chrom"] == c]
        body = {st: coverage(list(zip(g.loc[g.strand == st, "start"], g.loc[g.strand == st, "end"])), L)
                for st in ("+", "-")}
        gc_ = g[g["gene_type"].isin(CATEGORIES)]
        cat_cov = {st: coverage(list(zip(gc_.loc[gc_.strand == st, "start"], gc_.loc[gc_.strand == st, "end"])), L)
                   for st in ("+", "-")}
        any_pad = coverage(list(zip(g["start"], g["end"])), L, pad=INTERGENIC_PAD) > 0

        for r in gc_.itertuples(index=False):
            s, e = r.start, r.end
            opp = "-" if r.strand == "+" else "+"
            sel = ok[s:e] & ~excl[s:e] & (body[opp][s:e] == 0) & (cat_cov[r.strand][s:e] == 1)
            if not sel.any():
                continue
            ts, ta = (p[s:e], m[s:e]) if r.strand == "+" else (m[s:e], p[s:e])
            rep = is_rep[s:e]
            row = {"chrom": c, "gene_id": r.gene_id, "gene_name": r.gene_name, "gene_type": r.gene_type,
                   "strand": r.strand, "start": s, "end": e}
            for cls, cm in (("unique", sel & ~rep), ("all", sel)):
                count(row, cls, ts, ta, cm)
            b = up[s:e][sel & ~rep]
            cnt = {x: int((b == ord(x)).sum()) for x in "ACGT"}
            if r.strand == "-":
                cnt = {"A": cnt["T"], "C": cnt["G"], "G": cnt["C"], "T": cnt["A"]}
            row.update({f"sense_{x}": v for x, v in cnt.items()})
            gene_rows.append(row)

        ig = ok & ~excl & ~any_pad
        for b0 in range(0, L, BLOCK):
            b1 = min(L, b0 + BLOCK)
            sel = ig[b0:b1]
            if not sel.any():
                continue
            rep = is_rep[b0:b1]
            row = {"chrom": c, "block_start": b0}
            for cls, cm in (("unique", sel & ~rep), ("all", sel)):
                count(row, cls, p[b0:b1], m[b0:b1], cm)
            bb = up[b0:b1][sel & ~rep]
            row.update({f"sense_{x}": int((bb == ord(x)).sum()) for x in "ACGT"})
            block_rows.append(row)
        del seq, is_rep, up, ok, excl, body, cat_cov, any_pad

    G = pd.DataFrame(gene_rows)
    B = pd.DataFrame(block_rows)
    for df in (G, B):
        df["at_skew"] = (df["sense_T"] - df["sense_A"]) / (df["sense_T"] + df["sense_A"]).replace(0, np.nan)
        df["gc_skew"] = (df["sense_G"] - df["sense_C"]) / (df["sense_G"] + df["sense_C"]).replace(0, np.nan)
    G.to_csv(a.out_prefix + "_genes.tsv", sep="\t", index=False)
    B.to_csv(a.out_prefix + "_blocks.tsv", sep="\t", index=False)

    rows = []

    def add(section, label, df, cls, k, unit):
        S = df[f"{cls}_{k}_sense"].to_numpy(float)
        A = df[f"{cls}_{k}_anti"].to_numpy(float)
        R, lo, hi = ratio_ci(S, A, a.n_boot, rng)
        tot = S + A
        top1 = np.sort(tot)[::-1][:max(1, int(np.ceil(0.01 * len(tot))))].sum() / tot.sum() if tot.sum() else np.nan
        keep = np.argsort(-np.abs(S - A))[DROP_TOP:]
        Rd = S[keep].sum() / A[keep].sum() if A[keep].sum() else np.nan
        rows.append({"section": section, "set": label, "unit": unit, "n_units": len(df), "seq": cls, "count": k,
                     "sense": S.sum(), "antisense": A.sum(),
                     "deficit_pct": 100 * (1 - R), "deficit_lo95": 100 * (1 - hi), "deficit_hi95": 100 * (1 - lo),
                     "share_in_top1pct_units": top1, f"deficit_pct_without_top{DROP_TOP}": 100 * (1 - Rd)})

    pc = G[G["gene_type"] == "protein_coding"]
    lnc = G[G["gene_type"] == "lncRNA"]
    for cls in ("unique", "all"):
        for k in tkeys:
            for what in ("pos", "site"):
                kk = f"{k}_{what}"
                add("A", "intron_protein_coding", pc, cls, kk, "gene")
                add("A", "intron_lncRNA", lnc, cls, kk, "gene")
                add("A", "intergenic_plus_vs_minus", B, cls, kk, "1Mb_block")
    for st in ("+", "-"):
        for k in tkeys:
            add("C", f"intron_protein_coding_gene{st}", pc[pc["strand"] == st], "unique", f"{k}_site", "gene")
    for skew in ("at_skew", "gc_skew"):
        q = pc[skew].quantile([1 / 3, 2 / 3]).to_numpy()
        for i, (lo_, hi_) in enumerate([(-np.inf, q[0]), (q[0], q[1]), (q[1], np.inf)]):
            sub = pc[(pc[skew] > lo_) & (pc[skew] <= hi_)]
            add("D", f"intron_protein_coding_{skew}_tertile{i + 1}(median {sub[skew].median():+.3f})",
                sub, "unique", f"{tkeys[0]}_site", "gene")
    R = pd.DataFrame(rows)
    R.to_csv(a.out_prefix + "_summary.tsv", sep="\t", index=False)

    fmt = lambda x: f"{x:.3g}"
    cols = ["set", "n_units", "seq", "count", "sense", "antisense", "deficit_pct", "deficit_lo95",
            "deficit_hi95", "share_in_top1pct_units", f"deficit_pct_without_top{DROP_TOP}"]
    print("\n== A+B: positions vs merged sites; concentration (deficit > 0 = fewer on sense / on minus for intergenic) ==")
    sA = R[R["section"] == "A"]
    print(sA[cols].to_string(index=False, float_format=fmt))
    print("\n== C: protein-coding introns by gene strand (unique sequence, sites) ==")
    print(R[R["section"] == "C"][cols].to_string(index=False, float_format=fmt))
    print(f"\n== D: protein-coding introns by tertile of sense-strand composition skew (unique, sites, {tkeys[0]}) ==")
    print(R[R["section"] == "D"][cols].to_string(index=False, float_format=fmt))
    print(f"\nCIs: percentile bootstrap over units, {a.n_boot} resamples, seed {a.seed}; sites = above-threshold "
          f"positions merged within {MERGE} nt.")
    print(f"written: {a.out_prefix}_genes.tsv, _blocks.tsv, _summary.tsv")


if __name__ == "__main__":
    main()
