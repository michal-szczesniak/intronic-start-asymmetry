#!/usr/bin/env python3
"""Step 07b: how much of the observed strand asymmetry does Puffin-D reproduce?

Step 07 compared medians of predicted and observed DI directly. That comparison
is unfair to the model: Puffin was trained on log-transformed, averaged profiles
with a shape-based loss, so its outputs are compressed and their absolute scale
is not calibrated. A tiny predicted median says nothing until it is set against
the spread of the predictions. This step therefore asks the same question in
two scale-free ways, for observed DI and for each predicted DI metric in
parallel:

  1. Rank effect. AUC = P(case DI < control DI) + 0.5 P(tie). 0.5 = no shift;
     above 0.5 = the sense arm is relatively weaker in the case group.
     frac_rank = (AUC_pred - 0.5) / (AUC_obs - 0.5): how strongly the case
     group stands out in the predictions, relative to how strongly it stands out
     in the measurements. It is a ratio of standardised effects, not a share:
     measurements carry extra noise, so it runs higher than the share below.
  2. Calibration in observed units. Within control pairs, observed DI is
     regressed on predicted DI (linear, and by 20 quantile bins as a check).
     The fitted relation turns each case pair's predicted DI into an expected
     observed DI. frac_linear / frac_binned = expected mean shift / observed
     mean shift.

Control pairs are intergenic pairs taken in BOTH orientations (each pair enters
once as DI and once as -DI). This removes any plus-strand preference of the model
or the assay from the control, so the control centres exactly on zero.

Robustness readouts:
  - GC-stratified AUC (within quintiles of local GC defined on control pairs),
    because intronic pairs are AT-richer than intergenic ones.
  - Intronic pairs split by host-gene strand: a plus-strand preference would
    push the two halves in opposite directions.
  - Held-out chromosomes (default chr8, chr9: the Puffin test set), where the
    model never saw the measured signal, so it cannot reproduce the asymmetry
    by having memorised training loci. All quantities, including the
    calibration fit, are recomputed within this subset.

Reading frac_linear / frac_binned (the share): near 0 = the reference sequence,
as read by the model, does not account for the asymmetry; near 1 = it accounts
for all of it. The model misses part of what sequence encodes, so the share is a
conservative estimate of the sequence part, while a share from all chromosomes
can be inflated by memorisation; the held-out rows guard against that.

CIs: percentile bootstrap over pairs (cases and controls resampled separately).

Usage:
  python3 07b_calibrate.py --pred results/puffin_pairs_K562.tsv --line K562 \
      --out results/calib_K562.tsv [--heldout chr8,chr9] [--n-boot 2000] [--seed 1]
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
OBS = "DI"
METRICS = ["pred_DI", "pred_DI_point", "pred_DI_wide", "pred_DI_sum"]
N_BINS = 20
N_GC = 5
MIN_STRATUM_CTRL = 20


def auc_lower(case, ctrl_sorted):
    """P(case < ctrl) + 0.5 P(tie), averaged over case values."""
    n = len(ctrl_sorted)
    case = np.sort(case)  # sorted queries make searchsorted several times faster
    lo = np.searchsorted(ctrl_sorted, case, side="left")
    hi = np.searchsorted(ctrl_sorted, case, side="right")
    return float(np.mean((n - hi + 0.5 * (hi - lo)) / n))


def stats_once(c_obs, c_pred, c_gc, k_obs, k_pred, k_gc, gc_edges):
    """All readouts for one case set (c_*) against one control set (k_*, unsymmetrised)."""
    s_obs = np.concatenate([k_obs, -k_obs])
    s_pred = np.concatenate([k_pred, -k_pred])
    s_gc = np.concatenate([k_gc, k_gc])
    out = {}
    out["auc_obs"] = auc_lower(c_obs, np.sort(s_obs))
    out["auc_pred"] = auc_lower(c_pred, np.sort(s_pred))
    d_obs = out["auc_obs"] - 0.5
    out["frac_rank"] = (out["auc_pred"] - 0.5) / d_obs if d_obs != 0 else np.nan

    # GC-stratified AUC, weighted by the number of cases per stratum
    c_bin = np.digitize(c_gc, gc_edges)
    s_bin = np.digitize(s_gc, gc_edges)
    a_o, a_p, w = [], [], []
    for b in range(N_GC):
        mc, ms = c_bin == b, s_bin == b
        if mc.sum() == 0 or ms.sum() < MIN_STRATUM_CTRL:
            continue
        a_o.append(auc_lower(c_obs[mc], np.sort(s_obs[ms])))
        a_p.append(auc_lower(c_pred[mc], np.sort(s_pred[ms])))
        w.append(mc.sum())
    if w:
        out["auc_obs_gc"] = float(np.average(a_o, weights=w))
        out["auc_pred_gc"] = float(np.average(a_p, weights=w))
        d = out["auc_obs_gc"] - 0.5
        out["frac_rank_gc"] = (out["auc_pred_gc"] - 0.5) / d if d != 0 else np.nan
    else:
        out["auc_obs_gc"] = out["auc_pred_gc"] = out["frac_rank_gc"] = np.nan

    # calibration on symmetrised controls: means are 0 by construction, so OLS has no intercept
    sxx = float(np.sum(s_pred * s_pred))
    beta = float(np.sum(s_pred * s_obs)) / sxx if sxx > 0 else np.nan
    out["beta_ctrl"] = beta
    obs_shift = float(np.mean(c_obs))
    out["obs_mean_shift"] = obs_shift
    out["exp_shift_linear"] = beta * float(np.mean(c_pred))
    out["frac_linear"] = out["exp_shift_linear"] / obs_shift if obs_shift != 0 else np.nan

    edges = np.unique(np.quantile(s_pred, np.linspace(0, 1, N_BINS + 1)))
    if len(edges) >= 3:
        nb = len(edges) - 1
        s_b = np.clip(np.digitize(s_pred, edges[1:-1]), 0, nb - 1)
        cnt = np.bincount(s_b, minlength=nb)
        means = np.bincount(s_b, weights=s_obs, minlength=nb) / np.maximum(cnt, 1)
        c_b = np.clip(np.digitize(c_pred, edges[1:-1]), 0, len(edges) - 2)
        out["exp_shift_binned"] = float(np.mean(means[c_b]))
    else:
        out["exp_shift_binned"] = np.nan
    out["frac_binned"] = out["exp_shift_binned"] / obs_shift if obs_shift != 0 else np.nan
    return out


CI_KEYS = ["auc_obs", "auc_pred", "frac_rank", "frac_rank_gc", "frac_linear", "frac_binned"]


def analyse(case, ctrl, metric, n_boot, rng):
    c_obs, c_pred, c_gc = (case[OBS].to_numpy(float), case[metric].to_numpy(float),
                           case["gc_local"].to_numpy(float))
    k_obs, k_pred, k_gc = (ctrl[OBS].to_numpy(float), ctrl[metric].to_numpy(float),
                           ctrl["gc_local"].to_numpy(float))
    gc_edges = np.quantile(k_gc, np.linspace(0, 1, N_GC + 1))[1:-1]
    point = stats_once(c_obs, c_pred, c_gc, k_obs, k_pred, k_gc, gc_edges)
    boots = {k: [] for k in CI_KEYS}
    nc, nk = len(c_obs), len(k_obs)
    for _ in range(n_boot):
        ic = rng.integers(0, nc, nc)
        ik = rng.integers(0, nk, nk)
        b = stats_once(c_obs[ic], c_pred[ic], c_gc[ic], k_obs[ik], k_pred[ik], k_gc[ik], gc_edges)
        for k in CI_KEYS:
            boots[k].append(b[k])
    for k in CI_KEYS:
        v = np.asarray(boots[k], dtype=float)
        v = v[np.isfinite(v)]
        point[k + "_lo"] = float(np.percentile(v, 2.5)) if len(v) else np.nan
        point[k + "_hi"] = float(np.percentile(v, 97.5)) if len(v) else np.nan
    point["rho_ctrl"] = float(spearmanr(np.concatenate([k_obs, -k_obs]),
                                        np.concatenate([k_pred, -k_pred]))[0])
    return point


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pred", required=True, help="per-pair output of 07_puffin_pairs.py")
    ap.add_argument("--line", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--heldout", default="chr8,chr9", help="comma-separated chromosomes the model was not trained on")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    d = pd.read_csv(a.pred, sep="\t")
    need = ["chrom", "klasa", "host_signal", "strand_sense", "pred_status", "gc_local", OBS] + METRICS
    miss = [c for c in need if c not in d.columns]
    if miss:
        sys.exit(f"{a.pred}: missing columns {miss}")
    d = d[d["pred_status"] == "ok"].copy()
    d = d[np.isfinite(d[[OBS, "gc_local"] + METRICS].to_numpy(float)).all(axis=1)]
    heldout = set(a.heldout.split(","))

    intron = d["klasa"] == INTRON_CLASS
    active = d["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL
    case_defs = {
        "intron_all": intron,
        "intron_active": intron & active,
        "intron_silent": intron & (d["host_signal"] < ACTIVE_MIN_HOST_SIGNAL),
        "intron_active_hostplus": intron & active & (d["strand_sense"] == "+"),
        "intron_active_hostminus": intron & active & (d["strand_sense"] == "-"),
    }
    ctrl_mask = ~intron

    print(f"[{a.line}] {len(d)} scored pairs; control = {int(ctrl_mask.sum())} intergenic pairs, "
          f"both orientations; held-out set: {sorted(heldout)}")
    k = d[ctrl_mask]
    print("\n== spread in control pairs (IQR of DI, both orientations): how compressed are predictions? ==")
    for m in [OBS] + METRICS:
        v = np.concatenate([k[m].to_numpy(float), -k[m].to_numpy(float)])
        print(f"  {m:14s} IQR {np.percentile(v, 75) - np.percentile(v, 25):.4g}")
    print("  plus-strand preference in control (unsymmetrised median, plus arm = 'sense'): "
          + ", ".join(f"{m} {k[m].median():+.4g}" for m in [OBS, "pred_DI"]))

    rows = []
    for subset in ("all", "heldout"):
        sub = d if subset == "all" else d[d["chrom"].isin(heldout)]
        ctrl = sub[ctrl_mask.loc[sub.index]]
        for gname, gmask in case_defs.items():
            if subset == "heldout" and gname.endswith(("plus", "minus")):
                continue
            case = sub[gmask.loc[sub.index]]
            if len(case) < 20 or len(ctrl) < 50:
                print(f"  skip {subset}/{gname}: n_case={len(case)}, n_ctrl={len(ctrl)}")
                continue
            for m in METRICS:
                r = analyse(case, ctrl, m, a.n_boot, rng)
                rows.append({"line": a.line, "subset": subset, "group": gname, "metric": m,
                             "n_case": len(case), "n_ctrl": len(ctrl), **r})
    res = pd.DataFrame(rows)
    res["test"] = ("AUC = P(case < control) + 0.5 P(tie), controls in both orientations; "
                   "fractions = predicted / observed shift; percentile bootstrap over pairs")
    res.to_csv(a.out, sep="\t", index=False)

    show = ["subset", "group", "metric", "n_case", "auc_obs", "auc_pred",
            "frac_rank", "frac_rank_lo", "frac_rank_hi", "frac_rank_gc",
            "frac_linear", "frac_linear_lo", "frac_linear_hi", "frac_binned"]
    fmt = lambda x: f"{x:.3g}"
    print("\n== main metric (pred_DI) ==")
    print(res[res["metric"] == "pred_DI"][show].to_string(index=False, float_format=fmt))
    print("\n== intron_active, all chromosomes, every metric ==")
    sel = (res["subset"] == "all") & (res["group"] == "intron_active")
    print(res[sel][show].to_string(index=False, float_format=fmt))
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
