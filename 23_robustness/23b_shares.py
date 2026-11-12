#!/usr/bin/env python3
"""Step 23b (shares): share of the intronic asymmetry predicted from sequence, on the same pairs as
the measured asymmetry (intronic pairs in active hosts), with confidence intervals that resample
whole host genes.

Step 07b computed the share on all intronic pairs and resampled single pairs. Pairs in one host gene
are not independent (in K562, 1,863 pairs lie in about 1,100 genes), and the measured asymmetry of
Figure 1 uses intronic pairs in active hosts only. This script recomputes, for every prediction table,
with 07b's own statistics function (stats_once, imported unchanged and checked by md5):

  frac_linear  calibrated share (least squares through the origin on intergenic pairs taken in both
               orientations; expected mean shift of the intronic pairs / their observed mean shift);
  frac_rank    rank-based share, (AUC_pred - 0.5) / (AUC_obs - 0.5), which does not depend on calibration;
  frac_depth   calibrated share with the calibration done separately within three bands of read depth
               (tertiles of sens + anty among the intergenic pairs), each intronic pair calibrated in its
               own band - a check that calibration across depths does not drive the share.

Point estimates are 07b's (so frac_linear of the group intron_all reproduces 07b's tables). Confidence
intervals: percentile bootstrap, intronic pairs resampled by host gene, intergenic pairs by 1-Mb block
(chromosome, floor(mid / 1e6)), the two independently; the GC strata edges stay those of the full
control set, as in 07b.

Usage:
  python3 23b_shares.py --step07b 07_puffin_pairs/07b_calibrate.py --table LABEL:LINE=path [...] \
      --out results/23b_shares.tsv [--groups intron_active,intron_all] [--n-boot 2000] [--seed 1]
"""
import argparse
import hashlib
import importlib.util
import sys

import numpy as np
import pandas as pd

STEP07B_MD5 = "2d3d079426245a9022ca505ffd449372"
INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
N_DEPTH = 3


def load_module(name, path, md5):
    h = hashlib.md5(open(path, "rb").read()).hexdigest()
    if h != md5:
        sys.exit(f"STOP: {path} has md5 {h}, expected {md5}")
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def unit_slices(units):
    """Group row indices by unit: returns (order, offsets) so rows of unit u are order[off[u]:off[u+1]]."""
    u, inv = np.unique(np.asarray(units, str), return_inverse=True)
    order = np.argsort(inv, kind="stable")
    off = np.r_[0, np.cumsum(np.bincount(inv, minlength=len(u)))]
    return order, off, len(u)


def draw(order, off, n_units, rng):
    pick = rng.integers(0, n_units, n_units)
    return np.concatenate([order[off[p]:off[p + 1]] for p in pick])


def depth_share(c_obs, c_pred, c_dep, k_obs, k_pred, k_dep, edges):
    cb, kb = np.digitize(c_dep, edges), np.digitize(k_dep, edges)
    exp = np.zeros(len(c_obs))
    for b in range(N_DEPTH):
        sp = np.concatenate([k_pred[kb == b], -k_pred[kb == b]])
        so = np.concatenate([k_obs[kb == b], -k_obs[kb == b]])
        sxx = float(np.sum(sp * sp))
        beta = float(np.sum(sp * so)) / sxx if sxx > 0 else np.nan
        exp[cb == b] = beta * c_pred[cb == b]
    om = float(np.mean(c_obs))
    return float(np.mean(exp)) / om if om != 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step07b", required=True)
    ap.add_argument("--table", action="append", required=True, help="LABEL:LINE=path")
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups", default="intron_active,intron_all")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    s07 = load_module("s07b", a.step07b, STEP07B_MD5)
    rng = np.random.default_rng(a.seed)
    groups = [g for g in a.groups.split(",") if g]
    rows = []
    for spec in a.table:
        lab, rest = spec.split(":", 1)
        line, path = rest.split("=", 1)
        d = pd.read_csv(path, sep="\t", low_memory=False, dtype={"host": str, "klasa": str, "chrom": str})
        need = ["chrom", "mid", "klasa", "host", "host_signal", "pred_status", "gc_local", "sens", "anty", "DI", "pred_DI"]
        miss = [c for c in need if c not in d.columns]
        if miss:
            sys.exit(f"{path}: missing columns {miss}")
        d = d[d["pred_status"] == "ok"].copy()
        d = d[np.isfinite(d[["DI", "pred_DI", "gc_local"]].to_numpy(float)).all(axis=1)]
        intr = (d["klasa"] == INTRON_CLASS).to_numpy()
        act = intr & (d["host_signal"].to_numpy(float) >= ACTIVE_MIN_HOST_SIGNAL)
        k = d[~intr]
        k_obs, k_pred, k_gc = (k[c].to_numpy(float) for c in ("DI", "pred_DI", "gc_local"))
        k_dep = (k["sens"] + k["anty"]).to_numpy(float)
        gc_edges = np.quantile(k_gc, np.linspace(0, 1, s07.N_GC + 1))[1:-1]
        dep_edges = np.quantile(k_dep, np.linspace(0, 1, N_DEPTH + 1))[1:-1]
        k_unit = (k["chrom"].astype(str) + ":" + (k["mid"].astype(np.int64) // 1_000_000).astype(str)).to_numpy()
        ko, koff, kn = unit_slices(k_unit)
        for g in groups:
            m = act if g == "intron_active" else intr if g == "intron_all" else None
            if m is None:
                sys.exit(f"unknown group {g}")
            c = d[m]
            c_obs, c_pred, c_gc = (c[x].to_numpy(float) for x in ("DI", "pred_DI", "gc_local"))
            c_dep = (c["sens"] + c["anty"]).to_numpy(float)
            co, coff, cn = unit_slices(c["host"].astype(str).to_numpy())
            pt = s07.stats_once(c_obs, c_pred, c_gc, k_obs, k_pred, k_gc, gc_edges)
            pt["frac_depth"] = depth_share(c_obs, c_pred, c_dep, k_obs, k_pred, k_dep, dep_edges)
            B = {x: [] for x in ("frac_linear", "frac_rank", "frac_depth", "obs_mean_shift")}
            for _ in range(a.n_boot):
                ic = draw(co, coff, cn, rng)
                ik = draw(ko, koff, kn, rng)
                b = s07.stats_once(c_obs[ic], c_pred[ic], c_gc[ic], k_obs[ik], k_pred[ik], k_gc[ik], gc_edges)
                b["frac_depth"] = depth_share(c_obs[ic], c_pred[ic], c_dep[ic], k_obs[ik], k_pred[ik], k_dep[ik], dep_edges)
                for x in B:
                    B[x].append(b[x])
            r = {"label": lab, "line": line, "group": g, "n_pairs": len(c), "n_host_genes": cn,
                 "n_intergenic": len(k), "n_blocks": kn, "obs_mean_shift": pt["obs_mean_shift"],
                 "beta_ctrl": pt["beta_ctrl"]}
            for x in ("frac_linear", "frac_rank", "frac_depth"):
                v = np.asarray(B[x], float)
                v = v[np.isfinite(v)]
                r[x] = pt[x]
                r[x + "_lo95"] = float(np.percentile(v, 2.5)) if len(v) else np.nan
                r[x + "_hi95"] = float(np.percentile(v, 97.5)) if len(v) else np.nan
            v = np.asarray(B["obs_mean_shift"], float)
            r["obs_mean_shift_lo95"], r["obs_mean_shift_hi95"] = float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
            rows.append(r)
            print(f"  {lab:22s} {line:7s} {g:14s} n={len(c):5d} in {cn:5d} genes | linear {100 * r['frac_linear']:6.1f}% "
                  f"[{100 * r['frac_linear_lo95']:6.1f}, {100 * r['frac_linear_hi95']:6.1f}] | rank {100 * r['frac_rank']:6.1f}% "
                  f"[{100 * r['frac_rank_lo95']:6.1f}, {100 * r['frac_rank_hi95']:6.1f}] | depth bands {100 * r['frac_depth']:6.1f}% "
                  f"[{100 * r['frac_depth_lo95']:6.1f}, {100 * r['frac_depth_hi95']:6.1f}]", flush=True)
    R = pd.DataFrame(rows)
    R["ci"] = (f"percentile bootstrap, {a.n_boot} resamples (seed {a.seed}): intronic pairs by host gene, "
               "intergenic pairs by 1-Mb block")
    R.to_csv(a.out, sep="\t", index=False)
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
