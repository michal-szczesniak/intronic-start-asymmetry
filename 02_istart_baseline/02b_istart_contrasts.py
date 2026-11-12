#!/usr/bin/env python3
"""Step 02b: does the intronic sense/antisense asymmetry depend on host activity?

The project premise is that the asymmetry is present whether or not the host gene
is transcribed, which would argue for a sequence (not chromatin) explanation.
Step 02 reproduced the per-group numbers but did not test group differences.
Here we add, per cell line:

  * bootstrap 95% CI for the median DI and for "sense weaker by X%" in each group;
  * two-sided Mann-Whitney U tests between groups (active vs silent intronic,
    each intronic group vs intergenic, and within active hosts k36 yes vs no);
  * signal depth per group (median of sens + anty), because low counts inflate |DI|
    and reduce power, which matters for the small silent-host group.

Groups follow iSTART/scripts/26_directionality.py: intronic pairs are
klasa == 'intron genu'; a host is active when host_signal >= 5.

Usage:
  python3 02b_istart_contrasts.py --table K562=input/istart_directionality_K562.tsv \
      --table HCT116=input/istart_directionality_HCT116.tsv \
      --out-groups results/istart_contrasts_groups.tsv \
      --out-tests results/istart_contrasts_tests.tsv [--n-boot 10000] [--seed 1]
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ACTIVE_MIN_HOST_SIGNAL = 5
INTRON_CLASS = "intron genu"
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]


def boot_median_ci(x, n_boot, rng):
    """Percentile bootstrap CI for the median."""
    x = np.asarray(x, dtype=float)
    if len(x) < 5:
        return np.nan, np.nan
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    meds = np.median(x[idx], axis=1)
    return np.percentile(meds, 2.5), np.percentile(meds, 97.5)


def groups_of(d):
    gi = d[d["klasa"] == INTRON_CLASS]
    act = gi[gi["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL]
    return {
        "intergenic": d[d["klasa"] != INTRON_CLASS],
        "intron_silent": gi[gi["host_signal"] < ACTIVE_MIN_HOST_SIGNAL],
        "intron_active": act,
        "intron_active_k36_yes": act[act["k36"] == True],   # noqa: E712
        "intron_active_k36_no": act[act["k36"] == False],   # noqa: E712
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True, help="LINE=path/to/directionality.tsv")
    ap.add_argument("--out-groups", required=True)
    ap.add_argument("--out-tests", required=True)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    CONTRASTS = [("intron_active", "intron_silent"),
                 ("intron_active", "intergenic"),
                 ("intron_silent", "intergenic"),
                 ("intron_active_k36_yes", "intron_active_k36_no")]
    grows, trows = [], []
    for spec in a.table:
        line, path = spec.split("=", 1)
        d = pd.read_csv(path, sep="\t")
        if list(d.columns) != EXPECTED_COLUMNS:
            sys.exit(f"{path}: unexpected columns {list(d.columns)}")
        d["k36"] = d["k36"].map({True: True, False: False, "True": True, "False": False})
        g = groups_of(d)
        for name, df in g.items():
            lo, hi = boot_median_ci(df["DI"], a.n_boot, rng)
            grows.append({
                "line": line, "group": name, "n": len(df),
                "median_DI": df["DI"].median() if len(df) else np.nan,
                "median_DI_lo95": lo, "median_DI_hi95": hi,
                "sense_weaker_pct": 100 * (1 - 2 ** df["DI"].median()) if len(df) else np.nan,
                "sense_weaker_pct_lo95": 100 * (1 - 2 ** hi), "sense_weaker_pct_hi95": 100 * (1 - 2 ** lo),
                "pct_DI_zero": 100 * (df["DI"] == 0).mean() if len(df) else np.nan,
                "median_total_signal": (df["sens"] + df["anty"]).median() if len(df) else np.nan,
            })
        for x, y in CONTRASTS:
            dx, dy = g[x]["DI"], g[y]["DI"]
            if len(dx) < 5 or len(dy) < 5:
                continue
            u, p = mannwhitneyu(dx, dy, alternative="two-sided")
            trows.append({"line": line, "group_a": x, "group_b": y, "n_a": len(dx), "n_b": len(dy),
                          "median_DI_a": dx.median(), "median_DI_b": dy.median(),
                          "median_diff_a_minus_b": dx.median() - dy.median(),
                          "mannwhitney_p": p,
                          "test": "two-sided Mann-Whitney U on DI"})
    gout, tout = pd.DataFrame(grows), pd.DataFrame(trows)
    gout.to_csv(a.out_groups, sep="\t", index=False, float_format="%.4g")
    tout.to_csv(a.out_tests, sep="\t", index=False, float_format="%.4g")
    print(gout.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print()
    print(tout.drop(columns=["test"]).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples, seed {a.seed}. "
          f"Contrast test: {tout['test'].iloc[0]}.")
    print(f"written: {a.out_groups}; {a.out_tests}")


if __name__ == "__main__":
    main()
