#!/usr/bin/env python3
"""Step 02: reproduce the iSTART directionality numbers before any new analysis.

Re-computes, from the frozen iSTART tables (directionality.tsv, made by
iSTART/scripts/26_directionality.py), the same summaries that script prints:
for each cell line and group, n, median DI, % of pairs with DI < 0,
sense/antisense ratio (2^median DI) and a two-sided Wilcoxon signed-rank
test of DI against 0.

Groups (definitions copied from 26_directionality.py):
  intergenic      klasa != 'intron genu' (no protein-coding gene within 5 kb;
                  arm orientation arbitrary)
  intron_silent   klasa == 'intron genu' and host_signal < 5
  intron_active   klasa == 'intron genu' and host_signal >= 5
  intron_active split into tertiles of host_tx (rank-based qcut, as in step 26)
  intron_active split by k36 (H3K36me3 peak within 1 kb)

DI = log2((sense + 1) / (antisense + 1)); sense = arm on the host strand.

Usage:
  python3 02_istart_baseline.py --table K562=input/istart_directionality_K562.tsv \
      --table HCT116=input/istart_directionality_HCT116.tsv --out results/istart_baseline_summary.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

ACTIVE_MIN_HOST_SIGNAL = 5
INTRON_CLASS = "intron genu"
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]


def to_bool(series):
    """k36 is written as True/False text; NaN stays NaN."""
    return series.map({True: True, False: False, "True": True, "False": False})


def summarise(line, group, df):
    n = len(df)
    row = {"line": line, "group": group, "n": n, "median_DI": np.nan, "pct_DI_negative": np.nan,
           "sense_to_antisense_ratio": np.nan, "sense_weaker_pct": np.nan, "wilcoxon_p": np.nan,
           "test": "two-sided Wilcoxon signed-rank, H0: median DI = 0 (zeros dropped)"}
    if n < 5:
        return row
    med = df["DI"].median()
    row["median_DI"] = med
    row["pct_DI_negative"] = 100 * (df["DI"] < 0).mean()
    row["sense_to_antisense_ratio"] = 2 ** med
    row["sense_weaker_pct"] = 100 * (1 - 2 ** med)
    if (df["DI"] != 0).any():
        row["wilcoxon_p"] = wilcoxon(df["DI"]).pvalue
    return row


def analyse(line, path):
    d = pd.read_csv(path, sep="\t")
    if list(d.columns) != EXPECTED_COLUMNS:
        sys.exit(f"{path}: unexpected columns {list(d.columns)}")
    counts = d["klasa"].value_counts()
    print(f"\n[{line}] {path}: {len(d)} rows; classes: " +
          ", ".join(f"{k}={v}" for k, v in counts.items()))
    if len(counts) != 2 or INTRON_CLASS not in counts:
        sys.exit(f"{path}: expected exactly two classes including '{INTRON_CLASS}'")
    d["k36"] = to_bool(d["k36"])
    gi = d[d["klasa"] == INTRON_CLASS]
    act = gi[gi["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL]
    sil = gi[gi["host_signal"] < ACTIVE_MIN_HOST_SIGNAL]
    rows = [summarise(line, "intergenic", d[d["klasa"] != INTRON_CLASS]),
            summarise(line, "intron_silent", sil),
            summarise(line, "intron_active", act)]
    extra = {}
    if len(act) >= 15:
        q = pd.qcut(act["host_tx"].rank(method="first"), 3, labels=["low", "mid", "high"])
        for lev in ["low", "mid", "high"]:
            rows.append(summarise(line, f"intron_active_host_tx_{lev}", act[q == lev]))
        r, p = spearmanr(act["host_tx"], act["DI"])
        extra = {"spearman_r_hosttx_DI": r, "spearman_p": p, "n_active": len(act)}
    rows.append(summarise(line, "intron_active_k36_yes", act[act["k36"] == True]))  # noqa: E712
    rows.append(summarise(line, "intron_active_k36_no", act[act["k36"] == False]))  # noqa: E712
    return rows, extra


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True, help="LINE=path/to/directionality.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    all_rows = []
    for spec in a.table:
        line, path = spec.split("=", 1)
        rows, extra = analyse(line, path)
        all_rows += rows
        if extra:
            print(f"[{line}] Spearman(host_tx, DI) in active hosts: r = {extra['spearman_r_hosttx_DI']:.2f}, "
                  f"p = {extra['spearman_p']:.1e}, n = {extra['n_active']}")
    out = pd.DataFrame(all_rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    show = out.drop(columns=["test"]).copy()
    print("\n" + show.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nTest for every row: {out['test'].iloc[0]}")
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
