#!/usr/bin/env python3
"""Step 02c: is the intronic DI deficit an artefact of lower signal depth?

Intronic pairs carry less capped signal than intergenic ones (medians 58 vs 87 in
K562, 96 vs 124 in HCT116), and depth differences between compared groups are a
standard confounder. This script asks whether the deficit survives when the two
groups are made comparable in depth. Depth = sens + anty (capped reads summed in
+/-50 nt around both arms of a pair).

Three checks per cell line:
  1. Does DI depend on depth at all? Spearman of depth against DI and against
     |DI| within the intergenic group, where no asymmetry is expected; plus
     median DI per depth quintile in each group.
  2. Stratified comparison: intron_active vs intergenic within depth quintiles
     (cut points from the pooled distribution of both groups).
  3. Matched comparison: 1:1 nearest-neighbour matching of each intron_active
     pair to an unused intergenic pair on depth, with a relative caliper;
     Mann-Whitney U plus a bootstrap CI for the median DI difference.

Groups follow iSTART/scripts/26_directionality.py (intronic = klasa 'intron genu';
active host = host_signal >= 5).

Usage:
  python3 02c_depth_matched.py --table K562=input/istart_directionality_K562.tsv \
      --table HCT116=input/istart_directionality_HCT116.tsv \
      --out-depth results/istart_depth_profile.tsv \
      --out-strata results/istart_depth_strata.tsv \
      --out-matched results/istart_depth_matched.tsv [--caliper 0.1] [--n-boot 10000] [--seed 1]
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ACTIVE_MIN_HOST_SIGNAL = 5
INTRON_CLASS = "intron genu"
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]
N_STRATA = 5


def boot_diff_ci(x, y, n_boot, rng):
    """Percentile bootstrap CI for median(x) - median(y), resampling each group."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 5 or len(y) < 5:
        return np.nan, np.nan
    dx = np.median(x[rng.integers(0, len(x), size=(n_boot, len(x)))], axis=1)
    dy = np.median(y[rng.integers(0, len(y), size=(n_boot, len(y)))], axis=1)
    d = dx - dy
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


def match_on_depth(case_depth, ctrl_depth, caliper, rng):
    """Greedy 1:1 nearest-neighbour matching without replacement.

    Returns index arrays into case_depth and ctrl_depth. A match is kept only if
    |case - ctrl| <= caliper * case depth. Cases are processed in random order so
    that the result does not depend on genomic order.
    """
    order = np.argsort(ctrl_depth)
    sorted_depth = ctrl_depth[order]
    available = np.ones(len(ctrl_depth), dtype=bool)
    case_idx, ctrl_idx = [], []
    for i in rng.permutation(len(case_depth)):
        d = case_depth[i]
        pos = np.searchsorted(sorted_depth, d)
        best, best_gap = -1, np.inf
        lo, hi = pos - 1, pos
        while lo >= 0 or hi < len(sorted_depth):
            gap = np.inf
            if hi < len(sorted_depth):
                gap = min(gap, abs(sorted_depth[hi] - d))
            if lo >= 0:
                gap = min(gap, abs(d - sorted_depth[lo]))
            if gap > best_gap or gap > caliper * max(d, 1.0):
                break
            if hi < len(sorted_depth) and abs(sorted_depth[hi] - d) <= gap and available[order[hi]]:
                best, best_gap = order[hi], abs(sorted_depth[hi] - d)
                break
            if lo >= 0 and abs(d - sorted_depth[lo]) <= gap and available[order[lo]]:
                best, best_gap = order[lo], abs(d - sorted_depth[lo])
                break
            if hi < len(sorted_depth) and abs(sorted_depth[hi] - d) <= abs(d - sorted_depth[lo] if lo >= 0 else np.inf):
                hi += 1
            else:
                lo -= 1
        if best >= 0:
            available[best] = False
            case_idx.append(i)
            ctrl_idx.append(best)
    return np.array(case_idx, int), np.array(ctrl_idx, int)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True, help="LINE=path/to/directionality.tsv")
    ap.add_argument("--out-depth", required=True)
    ap.add_argument("--out-strata", required=True)
    ap.add_argument("--out-matched", required=True)
    ap.add_argument("--caliper", type=float, default=0.1, help="max relative depth difference (default 0.1)")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    depth_rows, strata_rows, matched_rows = [], [], []

    for spec in a.table:
        line, path = spec.split("=", 1)
        d = pd.read_csv(path, sep="\t")
        if list(d.columns) != EXPECTED_COLUMNS:
            sys.exit(f"{path}: unexpected columns {list(d.columns)}")
        d["depth"] = d["sens"] + d["anty"]
        gi = d[d["klasa"] == INTRON_CLASS]
        act = gi[gi["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL].copy()
        inter = d[d["klasa"] != INTRON_CLASS].copy()

        # 1. does DI depend on depth, in the group with no expected asymmetry?
        r_di, p_di = spearmanr(inter["depth"], inter["DI"])
        r_abs, p_abs = spearmanr(inter["depth"], inter["DI"].abs())
        print(f"\n[{line}] intergenic pairs (n = {len(inter)}): Spearman(depth, DI) r = {r_di:.3f}, p = {p_di:.1e}; "
              f"Spearman(depth, |DI|) r = {r_abs:.3f}, p = {p_abs:.1e}")

        pooled = pd.concat([act.assign(group="intron_active"), inter.assign(group="intergenic")])
        edges = np.unique(np.quantile(pooled["depth"], np.linspace(0, 1, N_STRATA + 1)))
        pooled["stratum"] = pd.cut(pooled["depth"], edges, include_lowest=True, labels=False)
        for g, df in pooled.groupby("group"):
            for s, sub in df.groupby("stratum"):
                depth_rows.append({"line": line, "group": g, "stratum": int(s), "n": len(sub),
                                   "depth_lo": edges[int(s)], "depth_hi": edges[int(s) + 1],
                                   "median_depth": sub["depth"].median(), "median_DI": sub["DI"].median(),
                                   "median_abs_DI": sub["DI"].abs().median()})
        for s, sub in pooled.groupby("stratum"):
            x = sub.loc[sub["group"] == "intron_active", "DI"]
            y = sub.loc[sub["group"] == "intergenic", "DI"]
            if len(x) < 5 or len(y) < 5:
                continue
            lo, hi = boot_diff_ci(x, y, a.n_boot, rng)
            strata_rows.append({"line": line, "stratum": int(s), "depth_lo": edges[int(s)],
                                "depth_hi": edges[int(s) + 1], "n_intron_active": len(x), "n_intergenic": len(y),
                                "median_DI_intron": x.median(), "median_DI_intergenic": y.median(),
                                "median_diff": x.median() - y.median(), "diff_lo95": lo, "diff_hi95": hi,
                                "sense_weaker_pct": 100 * (1 - 2 ** (x.median() - y.median())),
                                "mannwhitney_p": mannwhitneyu(x, y, alternative="two-sided").pvalue})

        # 3. matched comparison
        ci, ti = match_on_depth(act["depth"].to_numpy(float), inter["depth"].to_numpy(float), a.caliper, rng)
        x = act["DI"].to_numpy()[ci]
        y = inter["DI"].to_numpy()[ti]
        lo, hi = boot_diff_ci(x, y, a.n_boot, rng)
        diff = np.median(x) - np.median(y)
        matched_rows.append({
            "line": line, "n_intron_active": len(act), "n_matched": len(ci),
            "pct_matched": 100 * len(ci) / len(act), "caliper_rel": a.caliper,
            "median_depth_intron": float(np.median(act["depth"].to_numpy()[ci])),
            "median_depth_intergenic": float(np.median(inter["depth"].to_numpy()[ti])),
            "median_DI_intron": float(np.median(x)), "median_DI_intergenic": float(np.median(y)),
            "median_diff": diff, "diff_lo95": lo, "diff_hi95": hi,
            "sense_weaker_pct": 100 * (1 - 2 ** diff),
            "sense_weaker_pct_lo95": 100 * (1 - 2 ** hi), "sense_weaker_pct_hi95": 100 * (1 - 2 ** lo),
            "mannwhitney_p": mannwhitneyu(x, y, alternative="two-sided").pvalue})

    dd, ds, dm = pd.DataFrame(depth_rows), pd.DataFrame(strata_rows), pd.DataFrame(matched_rows)
    dd.to_csv(a.out_depth, sep="\t", index=False, float_format="%.4g")
    ds.to_csv(a.out_strata, sep="\t", index=False, float_format="%.4g")
    dm.to_csv(a.out_matched, sep="\t", index=False, float_format="%.4g")
    print("\n== median DI by depth quintile ==")
    print(dd.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\n== intron_active vs intergenic within depth quintiles ==")
    print(ds.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\n== depth-matched 1:1 comparison ==")
    print(dm.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples, seed {a.seed}; tests two-sided Mann-Whitney U on DI.")
    print(f"written: {a.out_depth}; {a.out_strata}; {a.out_matched}")


if __name__ == "__main__":
    main()
