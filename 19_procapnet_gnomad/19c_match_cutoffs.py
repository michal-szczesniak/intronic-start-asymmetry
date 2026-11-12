#!/usr/bin/env python3
"""Step 19c: change cut-offs for the ProCapNet variant classes, matched to Puffin-D's class sizes.

Step 12c sorts variants by the change they cause in the predicted signal, with Puffin-D cut-offs:
no effect = change below 0.001 on both strands; moderate gain or loss >= 0.003; strong >= 0.01;
"the other strand changes little" = below 0.003. ProCapNet predicts counts, on another scale, so
its cut-offs are set where the same share of changes passes them: for each Puffin-D cut-off x,
s(x) = share of |change| values (plus and minus strand of every scored variant on --chroms) that
are >= x; the ProCapNet cut-off is the value that the same share s(x) of ProCapNet's |change|
values for the same variants reaches. The site thresholds were matched the same way in step 16d
(quantiles of the predictions on chr8 + chr9).

Printed, for intronic variants of protein-coding genes: the size of each class under both
models with these cut-offs, and how many variants the two models place in the same class; for
comparison, the cut-offs scaled by the ratio of the weak-site thresholds.

Output <out>: quantity, puffin_cutoff, share_at_or_above, procapnet_cutoff, procapnet_share,
scaled_cutoff (read by 19b --cutoffs-file).

Usage:
  python3 19c_match_cutoffs.py --puffin-dir $SCRATCH/11d_scores --procapnet-dir $SCRATCH/19a_scores \\
      --chroms chr8,chr9 --pc-thresholds 1.131,2.658 --out results/19c_cutoffs.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PUFFIN = {"neutral": 0.001, "mid": 0.003, "other": 0.003, "strong": 0.01}
PUFFIN_THR = (0.0439, 0.1432)
USE = ["pos", "alt", "set", "strand", "ref_plus", "ref_minus", "alt_plus", "alt_minus"]


def load(path):
    return pd.read_csv(path, sep="\t", usecols=USE, keep_default_na=False, na_values=[],
                       dtype={"pos": np.int64, "alt": str, "set": str, "strand": str})


def classes(d, thr, cut):
    """The main step-12c classes, sense/antisense relative to the gene strand (plus for intergenic)."""
    plus_gene = (d["strand"] == "+").to_numpy() | (d["set"] == "intergenic").to_numpy()
    rp, rm, ap_, am = (d[k].to_numpy(float) for k in ("ref_plus", "ref_minus", "alt_plus", "alt_minus"))
    rs, ra = np.where(plus_gene, rp, rm), np.where(plus_gene, rm, rp)
    as_, aa = np.where(plus_gene, ap_, am), np.where(plus_gene, am, ap_)
    ds, da = as_ - rs, aa - ra
    T = thr[0]
    return {"no effect": (np.abs(ds) < cut["neutral"]) & (np.abs(da) < cut["neutral"]),
            "creates sense start (weak)": (rs < T) & (as_ >= T), "creates antisense start (weak)": (ra < T) & (aa >= T),
            "creates sense start (medium)": (rs < thr[1]) & (as_ >= thr[1]),
            "removes sense start (weak)": (rs >= T) & (as_ < T),
            "strong gain, sense": (ds >= cut["strong"]) & (da < cut["other"]),
            "strong gain, antisense": (da >= cut["strong"]) & (ds < cut["other"]),
            "moderate gain, sense": (ds >= cut["mid"]) & (ds < cut["strong"]) & (np.abs(da) < cut["other"]),
            "strong loss, sense": (ds <= -cut["strong"]) & (np.abs(da) < cut["other"])}, ds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--puffin-dir", required=True, help="step-11d score files <chrom>.scores.tsv.gz")
    ap.add_argument("--procapnet-dir", required=True, help="step-19a score files <chrom>.scores.tsv.gz")
    ap.add_argument("--chroms", default="chr8,chr9")
    ap.add_argument("--pc-thresholds", required=True, help="ProCapNet weak and medium site thresholds (step 16d)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    pc_thr = tuple(float(x) for x in a.pc_thresholds.split(","))
    P, C = [], []
    for c in [x for x in a.chroms.split(",") if x]:
        fp, fc = (os.path.join(d, f"{c}.scores.tsv.gz") for d in (a.puffin_dir, a.procapnet_dir))
        for f in (fp, fc):
            if not os.path.isfile(f):
                sys.exit(f"STOP: missing {f}")
        p, q = load(fp), load(fc)
        if len(p) != len(q) or not (np.array_equal(p["pos"].to_numpy(), q["pos"].to_numpy())
                                     and (p["alt"].to_numpy() == q["alt"].to_numpy()).all()):
            sys.exit(f"STOP: {c}: the two score files do not list the same variants in the same order")
        P.append(p)
        C.append(q)
        print(f"  {c}: {len(p):,} variants in both files", flush=True)
    P, C = pd.concat(P, ignore_index=True), pd.concat(C, ignore_index=True)
    absP = np.abs(np.r_[P["alt_plus"] - P["ref_plus"], P["alt_minus"] - P["ref_minus"]])
    absC = np.abs(np.r_[C["alt_plus"] - C["ref_plus"], C["alt_minus"] - C["ref_minus"]])
    M = len(absC)
    srt = np.sort(absC)
    scale = pc_thr[0] / PUFFIN_THR[0]
    rows = []
    for q, x in PUFFIN.items():
        share = float((absP >= x).mean())
        k = int(round(share * M))
        y = float(srt[M - k]) if k > 0 else float(srt[-1]) * 1.000001
        rows.append({"quantity": q, "puffin_cutoff": x, "share_at_or_above": share, "procapnet_cutoff": y,
                     "procapnet_share": float((absC >= y).mean()), "scaled_cutoff": x * scale})
    R = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    R.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(f"\n== cut-offs matched on {a.chroms} ({M:,} strand changes per model) ==")
    print(R.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    cutP = PUFFIN
    cutC = dict(zip(R["quantity"], R["procapnet_cutoff"]))
    pc = (P["set"] == "intron_protein_coding").to_numpy()
    kP, dsP = classes(P[pc], PUFFIN_THR, cutP)
    kC, dsC = classes(C[pc], pc_thr, cutC)
    print(f"\n== classes of {int(pc.sum()):,} intronic variants of protein-coding genes ==")
    print(f"  {'class':32s} {'Puffin-D':>10s} {'ProCapNet':>10s} {'both':>10s}")
    for k in kP:
        print(f"  {k:32s} {int(kP[k].sum()):10,d} {int(kC[k].sum()):10,d} {int((kP[k] & kC[k]).sum()):10,d}")
    rng = np.random.default_rng(1)
    s = rng.choice(len(dsP), size=min(len(dsP), 1_000_000), replace=False)
    print(f"\n  Spearman of the sense-strand change, Puffin-D vs ProCapNet (sample of {len(s):,}): "
          f"{spearmanr(dsP[s], dsC[s])[0]:+.3f}")
    print(f"written {a.out}")


if __name__ == "__main__":
    main()
