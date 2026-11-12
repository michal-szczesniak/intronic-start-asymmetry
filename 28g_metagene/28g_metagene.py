#!/usr/bin/env python3
"""Step 28g (12 Nov 2026; descriptive, for Figure 3F): strong measured start sites on the two strands of active
protein-coding genes in 1-kb bins of distance from the annotated 5' end of the gene.

No new site calling: the input is the site table that step 28b wrote (<prefix>_sites.tsv.gz: every strong site of
a protein-coding gene with its strand relative to the gene and its distance from the gene's 5' end) and the gene
table of step 28 (activity and span of every gene). Run for K562 PRO-cap and HCT116 PRO-cap (steps 28 and 28b)
and K562 GRO-cap (step 29b, the step-28 and step-28b scripts on GRO-cap).

For each 1-kb bin (1-2, 2-3, ..., 29-30 kb; eligible positions lie >= 1 kb from annotated transcript starts, so the
few sites below 1 kb are counted in the first bin, as in step 28b):
  sense, anti        strong sites on the gene strand and on the opposite strand
  genes_reaching     active genes whose annotated span is longer than the lower edge of the bin
  sense_per1000, anti_per1000   sites per 1,000 genes reaching the bin (the denominator of the density)
  ratio, lo, hi      sense / anti with a 95% CI from 2,000 resamples of genes within activity thirds (as step 28b)
  ratio3, lo3, hi3   the same in a 3-kb window centred on the bin (bins k-1, k, k+1; the first and last bins
                     use the two available bins)
A last row (bin_lo_kb = 30, bin_hi_kb empty) gives the sites at >= 30 kb: counts, genes longer than 30 kb, and the
ratio with its CI (no density, no window).
Checks (STOP otherwise): active genes as given (--expect-genes); all strong sites of active genes summed over
every distance equal step 28 (--expect-total S,A); sites below 5 kb equal the "1-5" stratum of step 28b
(--expect-first S,A).

Usage: python3 28g_metagene.py --label "K562 PRO-cap" --sites 28b_K562_sites.tsv.gz --genes 28_K562_genes.tsv \
          --expect-genes 10087 --expect-total 1224,1483 --expect-first S,A --out 28g_K562_bins.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd

ACT = ("A1", "A2", "A3")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--metric", default="strong")
    ap.add_argument("--max-kb", type=int, default=30)
    ap.add_argument("--expect-genes", type=int, required=True)
    ap.add_argument("--expect-total", required=True, help="sense,anti strong sites of active genes (step 28)")
    ap.add_argument("--expect-first", required=True, help="sense,anti below 5 kb (step 28b, stratum 1-5)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    G = pd.read_csv(a.genes, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str})
    G = G[(G["gene_type"] == "protein_coding") & G["activity"].isin(ACT)].reset_index(drop=True)
    print(f"== 28g {a.label}: {len(G):,} active protein-coding genes (expected {a.expect_genes:,})")
    if len(G) != a.expect_genes:
        sys.exit("STOP: number of active genes differs")
    if G.duplicated(["chrom", "gene_id"]).any():
        sys.exit("STOP: duplicated genes")
    key = {k: i for i, k in enumerate(zip(G["chrom"], G["gene_id"]))}
    glen_kb = (G["end"] - G["start"]).to_numpy(float) / 1e3

    S = pd.read_csv(a.sites, sep="\t", dtype={"gene_id": str})
    S = S[(S["metric"] == a.metric) & S["activity"].isin(ACT)]
    gi = np.array([key.get(k, -1) for k in zip(S["chrom"], S["gene_id"])])
    if (gi < 0).any():
        sys.exit(f"STOP: {int((gi < 0).sum())} sites in genes not in the gene table")
    sense = (S["strand"] == "sense").to_numpy()
    d = S["dist5_kb"].to_numpy(float)
    tot = (int(sense.sum()), int((~sense).sum()))
    first = (int((sense & (d < 5)).sum()), int((~sense & (d < 5)).sum()))
    exp_t = tuple(int(x) for x in a.expect_total.split(","))
    exp_f = tuple(int(x) for x in a.expect_first.split(","))
    print(f"  strong sites, sense / antisense: all {tot[0]:,} / {tot[1]:,} (step 28: {exp_t[0]:,} / {exp_t[1]:,}); "
          f"below 5 kb {first[0]:,} / {first[1]:,} (step 28b: {exp_f[0]:,} / {exp_f[1]:,})")
    if tot != exp_t or first != exp_f:
        sys.exit("STOP: steps 28 and 28b are not reproduced")
    below1 = int((d < 1).sum())
    beyond = int((d >= a.max_kb).sum())
    print(f"  sites below 1 kb (counted in the first bin): {below1}; sites at >= {a.max_kb} kb (not binned): {beyond:,}")

    nb = a.max_kb - 1
    b = np.floor(d).astype(int) - 1
    b[d < 1] = 0
    keep = d < a.max_kb
    CS = np.zeros((len(G), nb))
    CA = np.zeros((len(G), nb))
    np.add.at(CS, (gi[keep & sense], b[keep & sense]), 1)
    np.add.at(CA, (gi[keep & ~sense], b[keep & ~sense]), 1)

    act = G["activity"].to_numpy()
    parts = [np.flatnonzero(act == t) for t in ACT]
    Ws = []
    for ix in parts:
        idx = rng.integers(0, len(ix), size=(a.n_boot, len(ix)))
        W = np.zeros((a.n_boot, len(ix)))
        for k in range(a.n_boot):
            W[k] = np.bincount(idx[k], minlength=len(ix))
        Ws.append(W)
    bS = sum(W @ CS[ix] for W, ix in zip(Ws, parts))
    bA = sum(W @ CA[ix] for W, ix in zip(Ws, parts))

    def win(M, k):
        lo_, hi_ = max(k - 1, 0), min(k + 1, nb - 1)
        return M[..., lo_:hi_ + 1].sum(axis=-1)

    rows = []
    for k in range(nb):
        s, an = CS[:, k].sum(), CA[:, k].sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            r_b = bS[:, k] / bA[:, k]
            s3, a3 = win(CS.sum(axis=0), k), win(CA.sum(axis=0), k)
            r3_b = win(bS, k) / win(bA, k)
        r_b, r3_b = r_b[np.isfinite(r_b)], r3_b[np.isfinite(r3_b)]
        reach = int((glen_kb > k + 1).sum())
        rows.append({"data_set": a.label, "metric": a.metric, "bin_lo_kb": k + 1, "bin_hi_kb": k + 2,
                     "sense": int(s), "anti": int(an), "genes_reaching": reach,
                     "sense_per1000": 1000 * s / reach, "anti_per1000": 1000 * an / reach,
                     "ratio": s / an if an else np.nan,
                     "lo": np.percentile(r_b, 2.5), "hi": np.percentile(r_b, 97.5),
                     "sense3": int(s3), "anti3": int(a3), "ratio3": s3 / a3 if a3 else np.nan,
                     "lo3": np.percentile(r3_b, 2.5), "hi3": np.percentile(r3_b, 97.5)})
    far = d >= a.max_kb
    fs = np.zeros(len(G)); fa_ = np.zeros(len(G))
    np.add.at(fs, gi[far & sense], 1)
    np.add.at(fa_, gi[far & ~sense], 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rf = sum(W @ fs[ix] for W, ix in zip(Ws, parts)) / sum(W @ fa_[ix] for W, ix in zip(Ws, parts))
    rf = rf[np.isfinite(rf)]
    rows.append({"data_set": a.label, "metric": a.metric, "bin_lo_kb": a.max_kb, "bin_hi_kb": np.nan,
                 "sense": int(fs.sum()), "anti": int(fa_.sum()), "genes_reaching": int((glen_kb > a.max_kb).sum()),
                 "sense_per1000": np.nan, "anti_per1000": np.nan, "ratio": fs.sum() / fa_.sum(),
                 "lo": np.percentile(rf, 2.5), "hi": np.percentile(rf, 97.5), "sense3": np.nan, "anti3": np.nan,
                 "ratio3": np.nan, "lo3": np.nan, "hi3": np.nan})
    R = pd.DataFrame(rows)
    R.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(R[["bin_lo_kb", "sense", "anti", "genes_reaching", "ratio", "lo", "hi", "ratio3", "lo3", "hi3"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"written {a.out}")


if __name__ == "__main__":
    main()
