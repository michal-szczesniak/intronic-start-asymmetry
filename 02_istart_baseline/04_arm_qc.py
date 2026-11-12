#!/usr/bin/env python3
"""Step 04: quality control of the recovered arm coordinates, before any modelling.

Two things must be known before the pairs are handed to a sequence model of
initiation:

  1. How far apart are the two arms? The minimum observed distance is 0 nt, and
     arms that nearly coincide would be scored in essentially the same sequence
     window, so they cannot contribute an asymmetry.
  2. How many pairs sit at an annotated transcription start site? Step 26 kept
     intronic pairs at least 1,000 nt from host exons, but an annotated
     alternative promoter of another gene can still fall inside a host intron.
     Such pairs are real promoters rather than accidental starts and will need to
     be analysed separately.

Distances are reported per arm (sense and antisense) to the nearest annotated TSS
of any transcript in the GTF, and separately to the nearest protein-coding one.
Arm positions are 0-based (BED, from PINTS); GTF coordinates are 1-based and are
converted here.

Usage:
  python3 04_arm_qc.py --arms results/arms_K562.tsv --gtf ref/gencode.v47.annotation.gtf.gz \
      --out results/arms_K562.qc.tsv
"""
import argparse
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

INTRON_CLASS = "intron genu"


def open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_tss(path):
    """chrom -> sorted array of 0-based transcript TSS positions; all and protein-coding."""
    allt, pc = defaultdict(list), defaultdict(list)
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1
            allt[f[0]].append(tss)
            if 'gene_type "protein_coding"' in f[8]:
                pc[f[0]].append(tss)
    return ({c: np.array(sorted(v)) for c, v in allt.items()},
            {c: np.array(sorted(v)) for c, v in pc.items()})


def nearest(positions, chrom, pos):
    arr = positions.get(chrom)
    if arr is None or not len(arr):
        return np.nan
    i = np.searchsorted(arr, pos)
    best = np.inf
    if i < len(arr):
        best = min(best, abs(arr[i] - pos))
    if i > 0:
        best = min(best, abs(pos - arr[i - 1]))
    return float(best)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    d = pd.read_csv(a.arms, sep="\t")
    if "pos_sense" not in d.columns or "status" not in d.columns:
        sys.exit(f"{a.arms}: not an output of 03_arm_coordinates.py")
    d = d[d["status"] == "ok"].copy()
    line = d["line"].iloc[0]
    allt, pc = read_tss(a.gtf)
    print(f"[{line}] {sum(len(v) for v in allt.values())} transcript TSS, "
          f"{sum(len(v) for v in pc.values())} protein-coding, from {a.gtf}")

    for arm in ("sense", "antisense"):
        d[f"dist_tss_any_{arm}"] = [nearest(allt, c, p) for c, p in zip(d["chrom"], d[f"pos_{arm}"])]
        d[f"dist_tss_pc_{arm}"] = [nearest(pc, c, p) for c, p in zip(d["chrom"], d[f"pos_{arm}"])]
    d["dist_tss_any_min"] = d[["dist_tss_any_sense", "dist_tss_any_antisense"]].min(axis=1)
    d.to_csv(a.out, sep="\t", index=False)

    print(f"\n[{line}] arm distance (nt), percentiles by class:")
    pct = [0, 1, 5, 25, 50, 75, 95, 100]
    for k, g in d.groupby("klasa"):
        vals = np.percentile(g["arm_distance"], pct)
        print(f"  {k:<14} n={len(g):>5}  " + "  ".join(f"p{q}={v:.0f}" for q, v in zip(pct, vals)))
    print(f"\n[{line}] pairs with arm distance below a threshold:")
    for t in (1, 10, 25, 50, 100):
        sub = d[d["arm_distance"] < t]
        by = sub["klasa"].value_counts().to_dict()
        print(f"  < {t:>3} nt: {len(sub):>5} total  " + "  ".join(f"{k}={v}" for k, v in by.items()))
    print(f"\n[{line}] distance from either arm to the nearest annotated TSS (any transcript):")
    for t in (100, 500, 1000, 5000):
        sub = d[d["dist_tss_any_min"] < t]
        by = sub["klasa"].value_counts().to_dict()
        print(f"  < {t:>5} nt: {len(sub):>5} total  " + "  ".join(f"{k}={v}" for k, v in by.items()))
    print(f"\n[{line}] median distance to the nearest annotated TSS, by class and arm:")
    for k, g in d.groupby("klasa"):
        print(f"  {k:<14} sense: any={g['dist_tss_any_sense'].median():.0f} "
              f"pc={g['dist_tss_pc_sense'].median():.0f}   "
              f"antisense: any={g['dist_tss_any_antisense'].median():.0f} "
              f"pc={g['dist_tss_pc_antisense'].median():.0f}")
    print(f"\nwritten: {a.out}")


if __name__ == "__main__":
    main()
