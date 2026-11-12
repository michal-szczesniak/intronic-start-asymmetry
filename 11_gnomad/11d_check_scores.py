#!/usr/bin/env python3
"""Step 11d check: what a scores table from 11d_score_variants.py contains.
Prints counts, the size of predicted changes, how well the screen matched the exact
values (flagged and audit variants), how many audit variants the screen missed, and
how many variants create a site on the sense vs the antisense strand (a first look,
not the test). Usage: python3 11d_check_scores.py <scores.tsv.gz>"""
import sys

import numpy as np
import pandas as pd

FLAG = 0.003
THRESHOLDS = (0.0439, 0.1432)

d = pd.read_csv(sys.argv[1], sep="\t")
n0 = len(d)
d = d.drop_duplicates(["chrom", "pos", "alt"])
print(f"{n0:,} rows, {n0 - len(d):,} duplicates removed")
print("variants per set:", d["set"].value_counts().to_dict())
print("exact: 0 = screen only, 1 = flagged, 2 = audit:", d["exact"].value_counts().sort_index().to_dict())
dp, dm = d["alt_plus"] - d["ref_plus"], d["alt_minus"] - d["ref_minus"]
big = np.maximum(dp.abs(), dm.abs())
print("\n|change| (larger of the two strands), quantiles 50/90/99/99.9%:")
for s, g in big.groupby(d["set"]):
    print(f"  {s:22s} " + ", ".join(f"{np.quantile(g, q):.4f}" for q in (.5, .9, .99, .999)))
for k, lab in ((1, "flagged"), (2, "audit")):
    x = d[d["exact"] == k]
    if len(x):
        err = np.abs(np.r_[x["scr_alt_plus"] - x["alt_plus"], x["scr_alt_minus"] - x["alt_minus"]])
        print(f"\nscreen vs exact, {lab} variants (n = {len(x):,}): |error| 50/99% {np.median(err):.2e} / "
              f"{np.quantile(err, .99):.2e}, max {err.max():.2e}")
au = d[d["exact"] == 2]
if len(au):
    miss = (np.maximum((au["alt_plus"] - au["ref_plus"]).abs(), (au["alt_minus"] - au["ref_minus"]).abs()) > FLAG)
    print(f"audit: variants the screen did not flag but whose exact |change| > {FLAG}: {int(miss.sum())} of {len(au):,} "
          f"({100 * miss.mean():.2f}%)")
print("\nsite creation (reference below threshold, alternative at or above), first look:")
for T in THRESHOLDS:
    for s, g in d.groupby("set"):
        cp = (g["ref_plus"] < T) & (g["alt_plus"] >= T)
        cm = (g["ref_minus"] < T) & (g["alt_minus"] >= T)
        if s == "intergenic":
            print(f"  T={T} {s:22s} plus {int(cp.sum())}, minus {int(cm.sum())} (of {len(g):,})")
        else:
            plus_gene = g["strand"] == "+"
            sense = (cp & plus_gene) | (cm & ~plus_gene)
            anti = (cm & plus_gene) | (cp & ~plus_gene)
            print(f"  T={T} {s:22s} sense {int(sense.sum())}, antisense {int(anti.sum())} (of {len(g):,})")
