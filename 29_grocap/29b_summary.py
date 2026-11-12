#!/usr/bin/env python3
"""Step 29b summary: the three pre-specified questions for K562 GRO-cap, next to the same estimates from K562
PRO-cap (steps 28 and 28b), and the reading rule of the step-29 pre-specification.

  G1  deficit of strong sites in active protein-coding genes, all positions (step-28 table: subset active, metric
      strong, stratum all); prediction > 0
  G2  deficit in the 3' fifth minus the 5' fifth (step-28 table: stratum P5-P1, with two-sided bootstrap P);
      prediction < 0
  G3  deficit within 5 kb of the 5' end minus beyond (step-28b table: analysis sites_strong_bin5, stratum 1-5-rest,
      with two-sided bootstrap P); prediction > 0
Secondary rows: activity thirds, fifths, distance bins, inside minus outside H3K36me3 peaks, intergenic blocks.

Usage: python3 29b_summary.py <GRO-cap 28 strata> <GRO-cap 28b strata> <PRO-cap 28 strata> <PRO-cap 28b strata> <out prefix>
Writes <prefix>.tsv and prints the table and the reading.
"""
import os
import sys

import numpy as np
import pandas as pd

g28, g28b, p28, p28b, out = sys.argv[1:6]
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)


def pick28(df, subset, stratum, metric="strong"):
    m = df[(df["subset"] == subset) & (df["metric"] == metric) & (df["stratum"] == stratum)]
    if len(m) != 1:
        sys.exit(f"STOP: {subset}/{stratum}/{metric}: {len(m)} rows in a step-28 table")
    r = m.iloc[0]
    return float(r["deficit"]), float(r["lo"]), float(r["hi"]), float(r["p_boot"]) if pd.notna(r["p_boot"]) else np.nan, \
        (int(r["sense"]) if pd.notna(r["sense"]) else -1), (int(r["anti"]) if pd.notna(r["anti"]) else -1)


def pick28b(df, analysis, stratum):
    m = df[(df["analysis"] == analysis) & (df["stratum"] == stratum)]
    if len(m) != 1:
        sys.exit(f"STOP: {analysis}/{stratum}: {len(m)} rows in a step-28b table")
    r = m.iloc[0]
    return float(r["estimate"]), float(r["lo"]), float(r["hi"]), float(r["p_boot"]) if pd.notna(r["p_boot"]) else np.nan, \
        (int(r["sense"]) if pd.notna(r["sense"]) else -1), (int(r["anti"]) if pd.notna(r["anti"]) else -1)


T = {"GRO-cap": (pd.read_csv(g28, sep="\t"), pd.read_csv(g28b, sep="\t")),
     "PRO-cap": (pd.read_csv(p28, sep="\t"), pd.read_csv(p28b, sep="\t"))}
ITEMS = [("G1", "deficit, active genes, all positions", "28", ("active", "all")),
         ("G2", "3' fifth minus 5' fifth", "28", ("active", "P5-P1")),
         ("G3", "within 5 kb minus beyond", "28b", ("sites_strong_bin5", "1-5-rest"))]
SEC = [(f"third {t}", "28", (t, "all")) for t in ("A1", "A2", "A3")] + \
      [(f"fifth {p}", "28", ("active", p)) for p in ("P1", "P2", "P3", "P4", "P5")] + \
      [(f"bin {b} kb", "28b", ("sites_strong_bin5", b)) for b in ("1-5", "5-10", "10-20", "20-40", "40-80", ">=80")] + \
      [("inside minus outside H3K36me3 peaks", "28", ("active", "in-out")),
       ("silent genes", "28", ("silent", "all")), ("intergenic 1-Mb blocks", "28", ("intergenic", "all"))]
rows = []
for kind, items in (("primary", [(i[0] + ": " + i[1], i[2], i[3]) for i in ITEMS]), ("secondary", SEC)):
    for label, src, key in items:
        for data, (t28, t28b) in T.items():
            v = pick28(t28, *key) if src == "28" else pick28b(t28b, *key)
            rows.append({"kind": kind, "item": label, "data": data, "estimate": v[0], "lo95": v[1], "hi95": v[2],
                         "p_boot": v[3], "sense_sites": v[4], "antisense_sites": v[5]})
R = pd.DataFrame(rows)
R.to_csv(out + ".tsv", sep="\t", index=False, float_format="%.6g")
pd.set_option("display.width", 200)
print(R.round(2).to_string(index=False))


def g(k):
    r = R[(R["data"] == "GRO-cap") & R["item"].str.startswith(k + ":")].iloc[0]
    return r["estimate"], r["lo95"], r["hi95"]


g1, g2, g3 = g("G1"), g("G2"), g("G3")
ok1 = g1[1] > 0
ok2 = g2[2] < 0
ok3 = g3[1] > 0
print(f"\nG1 {'CI above 0' if ok1 else 'CI includes 0 or below'}; G2 {'CI below 0' if ok2 else 'CI includes 0 or above'}; "
      f"G3 {'CI above 0' if ok3 else 'CI includes 0 or below'}")
if ok1 and ok2 and ok3:
    print("READING (pre-specified rule): replicated by a second method (deficit and its position near the start)")
elif ok1:
    print("READING (pre-specified rule): deficit replicated; position not resolved")
else:
    print("READING (pre-specified rule): deficit not replicated")
print(f"written: {out}.tsv")
