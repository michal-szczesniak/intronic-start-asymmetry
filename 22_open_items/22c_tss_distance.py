#!/usr/bin/env python3
"""Step 22c: distance from the iSTART pairs to the nearest annotated transcription start.

Reads the step-04 tables (arms_<line>.qc.tsv), whose distance columns were computed by
04_arm_qc.py: for each arm, the distance (nt) to the nearest annotated TSS of any GENCODE v47
transcript (dist_tss_any_<arm>) and of a protein-coding transcript (dist_tss_pc_<arm>);
dist_tss_any_min is the smaller of the two arms. Here the protein-coding minimum over the two
arms is added, and every measure is summarised per cell line and group.

Groups: intron_all (klasa 'intron genu'), intron_active (and host_signal >= 5), intron_silent,
intergenic (all other pairs). Only pairs with resolved arm coordinates are in the step-04 tables.

Usage:
  python3 22c_tss_distance.py --arms K562=arms_K562.qc.tsv --arms HCT116=arms_HCT116.qc.tsv \
      --out results/22c_tss_distance.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
NEED = ["klasa", "host_signal", "status", "dist_tss_any_sense", "dist_tss_any_antisense",
        "dist_tss_pc_sense", "dist_tss_pc_antisense", "dist_tss_any_min"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", action="append", required=True, help="LINE=path/to/arms_LINE.qc.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for spec in a.arms:
        line, path = spec.split("=", 1)
        d = pd.read_csv(path, sep="\t")
        miss = [c for c in NEED if c not in d.columns]
        if miss:
            sys.exit(f"{path}: missing columns {miss}")
        if (d["status"] != "ok").any():
            sys.exit(f"{path}: rows with status other than 'ok' (expected only resolved pairs)")
        chk = np.nanmax(np.abs(d[["dist_tss_any_sense", "dist_tss_any_antisense"]].min(axis=1) - d["dist_tss_any_min"]))
        if chk > 0:
            sys.exit(f"{path}: dist_tss_any_min is not the minimum of the two arms (max difference {chk})")
        d["dist_tss_pc_min"] = d[["dist_tss_pc_sense", "dist_tss_pc_antisense"]].min(axis=1)
        intr = d["klasa"] == INTRON_CLASS
        act = intr & (d["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL)
        groups = {"intron_all": intr, "intron_active": act, "intron_silent": intr & ~act, "intergenic": ~intr}
        print(f"[{line}] {path}: {len(d):,} pairs; " + ", ".join(f"{g} {int(m.sum()):,}" for g, m in groups.items()))
        for g, m in groups.items():
            for meas in ["dist_tss_any_min", "dist_tss_pc_min", "dist_tss_any_sense", "dist_tss_any_antisense"]:
                x = d.loc[m, meas].dropna().to_numpy()
                rows.append({"line": line, "group": g, "measure": meas, "n": len(x),
                             "n_missing": int(m.sum()) - len(x),
                             "median_nt": float(np.median(x)) if len(x) else np.nan,
                             "q25_nt": float(np.percentile(x, 25)) if len(x) else np.nan,
                             "q75_nt": float(np.percentile(x, 75)) if len(x) else np.nan,
                             "n_within_100": int((x < 100).sum()), "pct_within_100": 100 * (x < 100).mean() if len(x) else np.nan,
                             "n_within_1000": int((x < 1000).sum()), "pct_within_1000": 100 * (x < 1000).mean() if len(x) else np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    show = out[out["measure"].isin(["dist_tss_any_min", "dist_tss_pc_min"])]
    print("\n" + show.to_string(index=False, float_format=lambda v: f"{v:.1f}"))
    print("\nwithin_100 / within_1000 = strictly less than 100 / 1,000 nt, as in 04_arm_qc.py")
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
