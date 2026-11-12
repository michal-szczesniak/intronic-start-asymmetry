#!/usr/bin/env python3
"""Step 09: H1 on whole introns - is predicted initiation weaker on the sense strand?

Question. Along the introns of a gene, does the DNA sequence (as read by Puffin-D,
strand-symmetric tracks from step 08) carry fewer potential start sites on the
gene's own strand (sense) than on the opposite strand (antisense)?

Which positions count as intronic for gene g (strand s):
  - inside the gene body of g, a protein-coding or lncRNA gene;
  - not in any annotated exon (any gene, either strand), nor within SPLICE_PAD nt
    of one (keeps splice sites, branch points and polypyrimidine tracts out);
  - not within TSS_PAD nt of any annotated transcript start (known promoters out);
  - not inside any gene body on the opposite strand (sense must be unambiguous);
  - not inside another protein-coding or lncRNA gene on the same strand (no
    position counted twice);
  - not N, and with a prediction on both strands.
Positions are split into unique sequence and repeats (soft-masked in hg38).

Intergenic positions (the control, which has no orientation): further than
INTERGENIC_PAD from any gene body, same exon/TSS/N exclusions. With strand-symmetric
predictions and no gene to orient them, plus and minus should balance; any
imbalance there is the residual bias of the whole procedure.

Readouts, per gene, sense vs antisense:
  - number of positions whose prediction exceeds a threshold; thresholds are the
    top 1e-3, 1e-4 and 1e-5 of all predictions on the chosen chromosomes (both
    strands pooled), i.e. increasingly strong potential start sites;
  - sum of predictions (dominated by low background; reported, not the main one).
Main statistic: R = sum over genes of sense / sum over genes of antisense, and
deficit = 100 * (1 - R). CI: percentile bootstrap over genes (genes, not
positions, are the independent units). Also a Wilcoxon signed-rank test on the
per-gene differences (sense - antisense). For the intergenic control the same R
is plus / minus, bootstrapped over 1 Mb blocks.

Sanity check before results: at annotated protein-coding gene starts, predicted
initiation should be stronger on the gene's strand than on the opposite strand
(mRNA promoters are sense-biased). A strand mix-up anywhere in the chain would
show up here as a fraction near or below 50%.

Also written: sense-strand base composition of intronic unique sequence
(AT skew = (T-A)/(T+A), GC skew = (G-C)/(G+C)); transcription is known to skew
it, which is the mutational rival explanation tested in a later step.

Usage:
  python3 09_h1_intron_strand.py --tracks /scratch/.../08_tracks --chroms chr8,chr9 \
      --genome /path/hg38.fa --gtf ../02_istart_baseline/ref/gencode.v47.annotation.gtf.gz \
      --out-genes results/h1_genes_chr8_chr9.tsv --out-summary results/h1_summary_chr8_chr9.tsv
"""
import argparse
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

SPLICE_PAD = 100
TSS_PAD = 1000
INTERGENIC_PAD = 10000
BLOCK = 1_000_000
TAILS = (1e-3, 1e-4, 1e-5)
TSS_HW = 100
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
            s, e, st = int(f[3]) - 1, int(f[4]), f[6]   # GTF 1-based inclusive -> 0-based half-open
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


def boot_ratio(S, A, n_boot, rng):
    S, A = np.asarray(S, float), np.asarray(A, float)
    if len(S) < 2 or A.sum() == 0:
        return np.nan, np.nan
    idx = rng.integers(0, len(S), size=(n_boot, len(S)))
    a = A[idx].sum(axis=1)
    r = np.where(a > 0, S[idx].sum(axis=1) / np.where(a > 0, a, 1), np.nan)
    r = r[np.isfinite(r)]
    return float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks", required=True, help="output directory of step 08")
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-genes", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    chroms = a.chroms.split(",")

    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)

    print(f"reading {a.gtf} ...", flush=True)
    genes, exons, tss = parse_gtf(a.gtf, set(chroms))
    print(f"  {len(genes)} genes on {chroms}; "
          + ", ".join(f"{c}: {len(exons[c])} exon records, {len(tss[c])} transcripts" for c in chroms))

    tracks = {}
    for c in chroms:
        p = np.load(os.path.join(a.tracks, f"{c}.procap_plus.f16.npy")).astype(np.float32)
        m = np.load(os.path.join(a.tracks, f"{c}.procap_minus.f16.npy")).astype(np.float32)
        if len(p) != len(fa[c]):
            sys.exit(f"{c}: track length {len(p)} != genome length {len(fa[c])}")
        tracks[c] = (p, m)

    # thresholds from all predictions on these chromosomes, both strands pooled
    allv = np.concatenate([np.concatenate([p[np.isfinite(p)], m[np.isfinite(m)]]) for p, m in tracks.values()])
    thr = {t: float(np.quantile(allv, 1 - t)) for t in TAILS}
    del allv
    print("thresholds (top fraction of all predictions): "
          + ", ".join(f"{t:g} -> {v:.4f}" for t, v in thr.items()))
    tkeys = [f"top{t:g}" for t in TAILS]

    gene_rows, block_rows, tss_rows = [], [], []
    for c in chroms:
        p, m = tracks[c]
        L = len(p)
        seq = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        is_rep = seq >= 97                       # lowercase = soft-masked repeat
        up = np.where(is_rep, seq - 32, seq)     # upper-case codes
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

        # sanity check at protein-coding gene starts
        for r in g[g["gene_type"] == "protein_coding"].itertuples(index=False):
            t = r.start if r.strand == "+" else r.end - 1
            lo, hi = t - TSS_HW, t + TSS_HW + 1
            if lo < 0 or hi > L or not ok[lo:hi].all():
                continue
            ts, ta = (p, m) if r.strand == "+" else (m, p)
            tss_rows.append((float(ts[lo:hi].max()), float(ta[lo:hi].max())))

        # genes
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
            for cls, cm in (("unique", sel & ~rep), ("repeat", sel & rep)):
                vs, va = ts[cm], ta[cm]
                row[f"{cls}_nt"] = int(cm.sum())
                row[f"{cls}_sum_sense"], row[f"{cls}_sum_anti"] = float(vs.sum()), float(va.sum())
                for t, k in zip(TAILS, tkeys):
                    row[f"{cls}_{k}_sense"] = int((vs >= thr[t]).sum())
                    row[f"{cls}_{k}_anti"] = int((va >= thr[t]).sum())
            b = up[s:e][sel & ~rep]
            cnt = {x: int((b == ord(x)).sum()) for x in "ACGT"}
            if r.strand == "-":                  # sense strand of a minus-strand gene = complement
                cnt = {"A": cnt["T"], "C": cnt["G"], "G": cnt["C"], "T": cnt["A"]}
            row.update({f"sense_{x}": v for x, v in cnt.items()})
            gene_rows.append(row)

        # intergenic control, in 1 Mb blocks
        ig = ok & ~excl & ~any_pad
        for b0 in range(0, L, BLOCK):
            b1 = min(L, b0 + BLOCK)
            sel = ig[b0:b1]
            if not sel.any():
                continue
            rep = is_rep[b0:b1]
            row = {"chrom": c, "block_start": b0}
            for cls, cm in (("unique", sel & ~rep), ("repeat", sel & rep)):
                vp, vm = p[b0:b1][cm], m[b0:b1][cm]
                row[f"{cls}_nt"] = int(cm.sum())
                row[f"{cls}_sum_sense"], row[f"{cls}_sum_anti"] = float(vp.sum()), float(vm.sum())
                for t, k in zip(TAILS, tkeys):
                    row[f"{cls}_{k}_sense"] = int((vp >= thr[t]).sum())
                    row[f"{cls}_{k}_anti"] = int((vm >= thr[t]).sum())
            bb = up[b0:b1][sel & ~rep]
            row.update({f"sense_{x}": int((bb == ord(x)).sum()) for x in "ACGT"})
            block_rows.append(row)
        del seq, is_rep, up, ok, excl, body, cat_cov, any_pad

    G = pd.DataFrame(gene_rows)
    B = pd.DataFrame(block_rows)
    for df in (G, B):
        for k in tkeys + ["sum"]:
            for side in ("sense", "anti"):
                df[f"all_{k}_{side}"] = df[f"unique_{k}_{side}"] + df[f"repeat_{k}_{side}"]
    G.to_csv(a.out_genes, sep="\t", index=False)

    tr = np.array(tss_rows)
    print(f"\n== sanity check: protein-coding gene starts (max within +/-{TSS_HW} nt) ==")
    if len(tr):
        print(f"  n = {len(tr)}; gene strand stronger than opposite strand at {100 * (tr[:, 0] > tr[:, 1]).mean():.1f}% "
              f"(ties {100 * (tr[:, 0] == tr[:, 1]).mean():.1f}%); median sense {np.median(tr[:, 0]):.3f}, "
              f"antisense {np.median(tr[:, 1]):.3f}")

    srows = []

    def summarise(df, unit, label, chrom_label):
        for cls in ("unique", "repeat", "all"):
            for k in tkeys + ["sum"]:
                S, A = df[f"{cls}_{k}_sense"].to_numpy(float), df[f"{cls}_{k}_anti"].to_numpy(float)
                if A.sum() == 0:
                    continue
                R = S.sum() / A.sum()
                lo, hi = boot_ratio(S, A, a.n_boot, rng)
                diff = S - A
                nz = diff[diff != 0]
                p = float(wilcoxon(nz).pvalue) if unit == "gene" and len(nz) >= 10 else np.nan
                srows.append({"chroms": chrom_label, "set": label, "unit": unit, "n_units": len(df),
                              "seq": cls, "metric": k, "nt": int(df[f"{cls}_nt"].sum()) if cls != "all"
                              else int(df["unique_nt"].sum() + df["repeat_nt"].sum()),
                              "sense": S.sum(), "antisense": A.sum(), "ratio": R, "ratio_lo95": lo,
                              "ratio_hi95": hi, "deficit_pct": 100 * (1 - R),
                              "deficit_pct_lo95": 100 * (1 - hi), "deficit_pct_hi95": 100 * (1 - lo),
                              "pct_units_sense_lower": 100 * float((nz < 0).mean()) if len(nz) else np.nan,
                              "wilcoxon_p": p})

    for chrom_label, cs in [(",".join(chroms), chroms)] + ([(c, [c]) for c in chroms] if len(chroms) > 1 else []):
        for cat in CATEGORIES:
            sub = G[(G["gene_type"] == cat) & G["chrom"].isin(cs)]
            if len(sub):
                summarise(sub, "gene", f"intron_{cat}", chrom_label)
        subb = B[B["chrom"].isin(cs)]
        if len(subb):
            summarise(subb, "1Mb_block", "intergenic_plus_vs_minus", chrom_label)
    S = pd.DataFrame(srows)
    S["test"] = ("ratio = total sense / total antisense; CI percentile bootstrap over units "
                 "(genes, or 1 Mb blocks for intergenic); Wilcoxon signed-rank on per-gene differences")
    S.to_csv(a.out_summary, sep="\t", index=False)

    print("\n== sense-strand composition of unique intronic sequence (intergenic: plus strand) ==")
    for label, df in [(f"intron_{cat}", G[G["gene_type"] == cat]) for cat in CATEGORIES] + [("intergenic", B)]:
        n = {x: df[f"sense_{x}"].sum() for x in "ACGT"}
        if sum(n.values()):
            print(f"  {label:22s} AT skew {(n['T'] - n['A']) / (n['T'] + n['A']):+.4f}, "
                  f"GC skew {(n['G'] - n['C']) / (n['G'] + n['C']):+.4f}, "
                  f"GC {(n['G'] + n['C']) / sum(n.values()):.3f}")

    show = ["chroms", "set", "n_units", "seq", "metric", "nt", "sense", "antisense",
            "deficit_pct", "deficit_pct_lo95", "deficit_pct_hi95", "pct_units_sense_lower", "wilcoxon_p"]
    fmt = lambda x: f"{x:.3g}"
    main_mask = (S["chroms"] == ",".join(chroms)) & S["metric"].isin(tkeys)
    print("\n== main result: potential start sites, sense vs antisense (deficit > 0 = fewer on sense) ==")
    print(S[main_mask][show].to_string(index=False, float_format=fmt))
    if len(chroms) > 1:
        print("\n== per chromosome, unique sequence, protein-coding introns and intergenic ==")
        pm = (S["chroms"] != ",".join(chroms)) & (S["seq"] == "unique") & S["metric"].isin(tkeys) \
            & S["set"].isin(["intron_protein_coding", "intergenic_plus_vs_minus"])
        print(S[pm][show].to_string(index=False, float_format=fmt))
    print("\n== sum of predictions (background-dominated; secondary) ==")
    print(S[(S["chroms"] == ",".join(chroms)) & (S["metric"] == "sum")][show].to_string(index=False, float_format=fmt))
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out_genes}; {a.out_summary}")


if __name__ == "__main__":
    main()
