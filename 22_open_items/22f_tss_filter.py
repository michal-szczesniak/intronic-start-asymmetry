#!/usr/bin/env python3
"""Step 22f: write iSTART pair tables without pairs that sit at an annotated transcription start.

Robustness test for the depth-matched comparison of step 02c. Intergenic pairs lie near annotated
starts (mostly of non-coding transcripts) much more often than intronic pairs, so the reference group
could be enriched for promoters. This script removes, from both groups, every pair with either start
closer than a threshold to the start of any annotated GENCODE v47 transcript (dist_tss_any_min from
step 04), and writes tables in exactly the input format of 02c_depth_matched.py (same columns, same
row order), so that 02c can be run on them unchanged.

The step-04 table (arms_<line>.qc.tsv) and the 02c input (istart_directionality_<line>.tsv) hold the
same pairs in the same order; this is checked on chrom, mid, klasa, sens and anty before any row is used.

Usage:
  python3 22f_tss_filter.py --table K562=input/istart_directionality_K562.tsv --arms K562=results/arms_K562.qc.tsv \
      [--table HCT116=... --arms HCT116=...] --min-dist 100,1000 --outdir results/22f_tables
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]
KEYS = ["chrom", "mid", "klasa", "sens", "anty"]


def pairs(specs):
    out = {}
    for s in specs:
        line, path = s.split("=", 1)
        out[line] = path
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True)
    ap.add_argument("--arms", action="append", required=True)
    ap.add_argument("--min-dist", default="100,1000", help="comma-separated thresholds in nt")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    tables, arms = pairs(a.table), pairs(a.arms)
    if set(tables) != set(arms):
        sys.exit("--table and --arms must name the same lines")
    thresholds = [int(x) for x in a.min_dist.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for line in tables:
        d = pd.read_csv(tables[line], sep="\t")
        if list(d.columns) != EXPECTED_COLUMNS:
            sys.exit(f"{tables[line]}: unexpected columns {list(d.columns)}")
        q = pd.read_csv(arms[line], sep="\t")
        if len(q) != len(d):
            sys.exit(f"{line}: {len(d)} rows in the 02c input but {len(q)} in the step-04 table")
        for k in KEYS:
            x, y = d[k].to_numpy(), q[k].to_numpy()
            same = (x == y) | (pd.isna(x) & pd.isna(y)) if x.dtype.kind == "f" else (x == y)
            if not np.all(same):
                sys.exit(f"{line}: column {k} differs between the two tables at {int((~same).sum())} rows")
        dist = q["dist_tss_any_min"].to_numpy()
        if np.isnan(dist).any():
            sys.exit(f"{line}: missing distances")
        intr = (d["klasa"] == INTRON_CLASS).to_numpy()
        act = intr & (d["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL).to_numpy()
        print(f"[{line}] {len(d):,} pairs checked row by row against {arms[line]}")
        for t in thresholds:
            keep = dist >= t
            out = os.path.join(a.outdir, f"istart_directionality_{line}.min{t}.tsv")
            d.loc[keep].to_csv(out, sep="\t", index=False)
            for g, m in (("intron_active", act), ("intergenic", ~intr)):
                rows.append({"line": line, "min_dist_nt": t, "group": g, "n_all": int(m.sum()),
                             "n_kept": int((m & keep).sum()), "pct_removed": 100 * (1 - (m & keep).sum() / m.sum())})
            print(f"   min {t:>5} nt: kept {int(keep.sum()):,} of {len(d):,}; written {out}")
    s = pd.DataFrame(rows)
    s.to_csv(os.path.join(a.outdir, "22f_filter_counts.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\n" + s.to_string(index=False, float_format=lambda v: f"{v:.1f}"))


if __name__ == "__main__":
    main()
