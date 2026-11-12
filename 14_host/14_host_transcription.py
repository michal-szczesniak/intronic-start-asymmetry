#!/usr/bin/env python3
"""Step 14: does the part of the sense weakness that sequence does not explain depend on
the host gene's own transcription?

Steps 07b/07c: in intronic divergent PRO-cap pairs the sense arm is weaker; Puffin-D,
reading the sequence, accounts for roughly a quarter to a third of that. Candidate
explanations of the rest involve the host gene's transcription itself (Pol II passing
through the element, chromatin marked behind it such as H3K36me3, which in yeast and
mouse suppresses spurious initiation inside genes). If so, the residual sense weakness
should (1) grow with host activity and (2) follow the host's activity when the SAME
element is compared between cell lines.

Input: the per-pair tables of step 07c (results/puffin_rc_pairs_<line>.tsv; observed
DI = log2((sense+1)/(antisense+1)), predicted pred_DI, host gene, host_signal, host_tx,
k36), and the step-07b calibration tables (calib_rc_<line>.tsv) when present.

Part 1  What the host columns look like (distributions, agreement between them).
Part 2  The step-07b rows for active vs silent hosts (as already computed).
Part 3  Within each line, intronic pairs grouped by host activity:
          silent: host_signal = 0, 0 < host_signal < 5      (5 = the 07b threshold)
          active: quartiles of host_signal among active pairs (Q1..Q4)
          and, among active pairs, k36 True vs False.
        Per group: n, mean observed DI, sequence-expected DI = beta * mean pred_DI
        (beta: OLS of observed on predicted DI in intergenic pairs taken in both
        orientations, as 07b), residual = observed - expected (< 0: sense weaker than
        its sequence predicts). Contrasts: each group's residual minus the silent
        (host_signal = 0) residual. CIs: bootstrap over host genes (pairs of one gene
        move together), controls resampled by pair for beta.
Part 4  Natural experiment across cell lines. Intronic pairs found in both lines (both
        arms within --match-nt, same host gene) have the same sequence, so the
        sequence part cancels. dDI = DI_lineA - DI_lineB. Groups: host active in both /
        only A / only B / neither. Effect of host activity on DI within an element:
            (mean dDI[only A] - mean dDI[only B]) / 2
        (the overall difference between the lines cancels). < 0: the sense arm is
        weaker where the host is active. Also Spearman of dDI with the difference in
        host_signal percentile (within-line ranks), robust to scale differences.

Caveat printed with the results: pairs are found by peak calling, and a sense-strand
arm has to stand out against the host's own sense-strand signal, which is higher in
active genes; detection can therefore depend on host activity in both lines.

Usage:
  python3 14_host_transcription.py --pairs-dir ../07_puffin_pairs/results --lines K562,HCT116 \
      --out-prefix results/14_host [--match-nt 50] [--n-boot 2000] [--seed 1]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

INTRON = "intron genu"
ACTIVE = 5.0
OBS, PRED = "DI", "pred_DI"


def as_bool(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def load(pairs_dir, line):
    f = os.path.join(pairs_dir, f"puffin_rc_pairs_{line}.tsv")
    if not os.path.isfile(f):
        sys.exit(f"STOP: missing {f}")
    d = pd.read_csv(f, sep="\t", keep_default_na=False, na_values=["", "NA", "nan", "NaN"],
                    dtype={"host": str, "host_name": str, "host_strand": str, "klasa": str})
    need = ["chrom", "klasa", "host", "host_signal", "host_tx", "k36", "pos_plus", "pos_minus", "strand_sense",
            "pred_status", OBS, PRED]
    miss = [c for c in need if c not in d.columns]
    if miss:
        sys.exit(f"STOP: {f} lacks columns {miss}")
    n0 = len(d)
    fin = [c for c in [OBS, "gc_local", PRED, "pred_DI_point", "pred_DI_wide", "pred_DI_sum"] if c in d.columns]
    d = d[(d["pred_status"] == "ok") & np.isfinite(d[fin].to_numpy(float)).all(axis=1)].reset_index(drop=True)   # as 07b
    d["intron"] = d["klasa"] == INTRON
    d["k36b"] = as_bool(d["k36"])
    return d, n0


def boot_counts(n, n_boot, rng):
    idx = rng.integers(0, n, size=(n_boot, n))
    W = np.zeros((n_boot, n), dtype=np.float64)
    for b in range(n_boot):
        W[b] = np.bincount(idx[b], minlength=n)
    return W


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (np.nan, np.nan)


def pval(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return max(2 * min((v <= 0).mean(), (v >= 0).mean()), 1 / len(v)) if len(v) else np.nan


def pct_weaker(di):
    return 100 * (1 - 2 ** di)


def describe(x):
    x = pd.to_numeric(x, errors="coerce")
    q = x.quantile([0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1]).to_numpy()
    return (f"n={x.notna().sum():,} missing={x.isna().sum():,} zero={(x == 0).sum():,}; quantiles 0/10/25/50/75/90/99/100%: "
            + " / ".join(f"{v:.3g}" for v in q))


def within_line(d, line, n_boot, rng):
    """Part 3: residual sense weakness by host activity, bootstrap over host genes."""
    intr = d[d["intron"]].copy()
    ctrl = d[~d["intron"]]
    hs = pd.to_numeric(intr["host_signal"], errors="coerce")
    groups = {}
    groups["silent_0"] = (hs == 0).to_numpy()
    groups["silent_0to5"] = ((hs > 0) & (hs < ACTIVE)).to_numpy()
    act = (hs >= ACTIVE).to_numpy()
    if act.sum() >= 40:
        qs = np.quantile(hs[act], [0.25, 0.5, 0.75])
        lab = np.full(len(intr), -1)
        lab[act] = np.searchsorted(qs, hs[act].to_numpy(), side="right")
        for k in range(4):
            groups[f"active_Q{k + 1}"] = lab == k
        qtxt = ", ".join(f"{v:.3g}" for v in qs)
    else:
        qtxt = "too few active pairs for quartiles"
    groups["active_all"] = act
    groups["active_k36"] = act & intr["k36b"].to_numpy()
    groups["active_no_k36"] = act & ~intr["k36b"].to_numpy()
    groups["intron_all"] = np.ones(len(intr), dtype=bool)

    units = intr["host"].fillna("").astype(str).to_numpy()
    units = np.where(units == "", "pair" + intr.index.astype(str), units)
    uu, ui = np.unique(units, return_inverse=True)
    U = len(uu)
    o, p = intr[OBS].to_numpy(float), intr[PRED].to_numpy(float)
    ko, kp = ctrl[OBS].to_numpy(float), ctrl[PRED].to_numpy(float)
    Wg = boot_counts(U, n_boot, rng)
    Wc = boot_counts(len(ko), n_boot, rng)
    beta = float(np.sum(kp * ko) / np.sum(kp * kp))                 # symmetrised controls: sum doubles on both sides
    beta_b = (Wc @ (kp * ko)) / (Wc @ (kp * kp))
    stats, boots = {}, {}
    for g, m in groups.items():
        if m.sum() == 0:
            continue
        n_u = np.bincount(ui[m], minlength=U).astype(float)
        so = np.bincount(ui[m], weights=o[m], minlength=U)
        sp = np.bincount(ui[m], weights=p[m], minlength=U)
        nb = Wg @ n_u
        with np.errstate(invalid="ignore", divide="ignore"):
            mo_b, mp_b = (Wg @ so) / nb, (Wg @ sp) / nb
        mo, mp = so.sum() / n_u.sum(), sp.sum() / n_u.sum()
        res_b = mo_b - beta_b * mp_b
        boots[g] = res_b
        stats[g] = {"line": line, "group": g, "n_pairs": int(m.sum()), "n_genes": int((n_u > 0).sum()),
                    "host_signal_median": float(np.median(hs[m])) if m.sum() else np.nan,
                    "obs_mean_DI": mo, "obs_median_DI": float(np.median(o[m])),
                    "pct_sense_weaker_median": pct_weaker(float(np.median(o[m]))),
                    "pred_mean_DI": mp, "beta_ctrl": beta, "seq_expected_DI": beta * mp,
                    "residual_DI": mo - beta * mp, "residual_lo95": ci(res_b)[0], "residual_hi95": ci(res_b)[1],
                    "seq_share": (beta * mp) / mo if mo != 0 else np.nan}
    rows = list(stats.values())
    crows = []
    if "silent_0" in boots:
        for g in stats:
            if g in ("silent_0", "intron_all"):
                continue
            dv = boots[g] - boots["silent_0"]
            crows.append({"line": line, "contrast": f"{g} - silent_0",
                          "residual_diff": stats[g]["residual_DI"] - stats["silent_0"]["residual_DI"],
                          "lo95": ci(dv)[0], "hi95": ci(dv)[1], "p_boot": pval(dv)})
    if "active_k36" in boots and "active_no_k36" in boots:
        dv = boots["active_k36"] - boots["active_no_k36"]
        crows.append({"line": line, "contrast": "active_k36 - active_no_k36",
                      "residual_diff": stats["active_k36"]["residual_DI"] - stats["active_no_k36"]["residual_DI"],
                      "lo95": ci(dv)[0], "hi95": ci(dv)[1], "p_boot": pval(dv)})
    return pd.DataFrame(rows), pd.DataFrame(crows), qtxt


def match_lines(da, db, tol):
    """One-to-one matches of intronic pairs between two lines: both arms within tol nt, same host."""
    A = da[da["intron"]].reset_index(drop=True)
    B = db[db["intron"]].reset_index(drop=True)
    A = A.assign(ia=np.arange(len(A)))
    B = B.assign(ib=np.arange(len(B)))
    out = []
    for c, ga in A.groupby("chrom"):
        gb = B[B["chrom"] == c]
        if not len(gb):
            continue
        pa, ma = ga["pos_plus"].to_numpy(), ga["pos_minus"].to_numpy()
        pb, mb = gb["pos_plus"].to_numpy(), gb["pos_minus"].to_numpy()
        order = np.argsort(pb)
        pbs = pb[order]
        for i in range(len(ga)):
            lo, hi = np.searchsorted(pbs, pa[i] - tol, "left"), np.searchsorted(pbs, pa[i] + tol, "right")
            for j in order[lo:hi]:
                if abs(ma[i] - mb[j]) <= tol:
                    out.append((ga["ia"].iat[i], gb["ib"].iat[j], abs(pa[i] - pb[j]) + abs(ma[i] - mb[j])))
    if not out:
        return None
    M = pd.DataFrame(out, columns=["ia", "ib", "dist"]).sort_values(["dist", "ia", "ib"], kind="stable")
    M = M.drop_duplicates("ia").drop_duplicates("ib")
    M = M.merge(A, on="ia").merge(B, on="ib", suffixes=("_a", "_b"))
    M = M[M["host_a"].astype(str) == M["host_b"].astype(str)]
    hs_ok = pd.to_numeric(M["host_signal_a"], errors="coerce").notna() & pd.to_numeric(M["host_signal_b"], errors="coerce").notna()
    return M[hs_ok].reset_index(drop=True)


def cross_line(M, la, lb, da, db, n_boot, rng):
    """Part 4: within-element effect of host activity, bootstrap over host genes."""
    ha = pd.to_numeric(M["host_signal_a"], errors="coerce").to_numpy(float)
    hb = pd.to_numeric(M["host_signal_b"], errors="coerce").to_numpy(float)
    act_a, act_b = ha >= ACTIVE, hb >= ACTIVE
    dDI = (M[OBS + "_a"] - M[OBS + "_b"]).to_numpy(float)
    dpred = (M[PRED + "_a"] - M[PRED + "_b"]).to_numpy(float)
    # within-line percentile of host_signal among that line's intronic pairs
    ref_a = np.sort(pd.to_numeric(da.loc[da["intron"], "host_signal"], errors="coerce").dropna().to_numpy(float))
    ref_b = np.sort(pd.to_numeric(db.loc[db["intron"], "host_signal"], errors="coerce").dropna().to_numpy(float))
    pr_a = np.searchsorted(ref_a, ha, "right") / len(ref_a)
    pr_b = np.searchsorted(ref_b, hb, "right") / len(ref_b)
    cats = {"active_both": act_a & act_b, f"active_only_{la}": act_a & ~act_b,
            f"active_only_{lb}": ~act_a & act_b, "active_neither": ~act_a & ~act_b}
    units = M["host_a"].astype(str).to_numpy()
    uu, ui = np.unique(units, return_inverse=True)
    U = len(uu)
    W = boot_counts(U, n_boot, rng)
    rows, bm = [], {}
    for g, m in cats.items():
        if m.sum() == 0:
            continue
        n_u = np.bincount(ui[m], minlength=U).astype(float)
        s_u = np.bincount(ui[m], weights=dDI[m], minlength=U)
        with np.errstate(invalid="ignore", divide="ignore"):
            bm[g] = (W @ s_u) / (W @ n_u)
        lo, hi = ci(bm[g])
        rows.append({"quantity": f"mean dDI ({la} - {lb}), {g}", "n": int(m.sum()), "value": float(dDI[m].mean()),
                     "lo95": lo, "hi95": hi, "p_boot": np.nan})
    ka, kb = f"active_only_{la}", f"active_only_{lb}"
    if ka in bm and kb in bm:
        v = (dDI[cats[ka]].mean() - dDI[cats[kb]].mean()) / 2
        bv = (bm[ka] - bm[kb]) / 2
        rows.append({"quantity": "host-activity effect on DI within element = (only A - only B)/2",
                     "n": int(cats[ka].sum() + cats[kb].sum()), "value": float(v), "lo95": ci(bv)[0],
                     "hi95": ci(bv)[1], "p_boot": pval(bv)})
    rho, prho = spearmanr(dDI, pr_a - pr_b)
    rows.append({"quantity": f"Spearman(dDI, host_signal percentile {la} - {lb}), all matched", "n": len(M),
                 "value": float(rho), "lo95": np.nan, "hi95": np.nan, "p_boot": float(prho)})
    disc = cats[ka] | cats[kb]
    if disc.sum() >= 10:
        rho2, p2 = spearmanr(dDI[disc], (pr_a - pr_b)[disc])
        rows.append({"quantity": "Spearman as above, discordant elements only", "n": int(disc.sum()),
                     "value": float(rho2), "lo95": np.nan, "hi95": np.nan, "p_boot": float(p2)})
    rows.append({"quantity": "sanity: median |pred_DI difference| of matched elements (same sequence -> ~0)",
                 "n": len(M), "value": float(np.median(np.abs(dpred))), "lo95": np.nan, "hi95": np.nan, "p_boot": np.nan})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-dir", required=True, help="directory with puffin_rc_pairs_<line>.tsv (step 07c)")
    ap.add_argument("--lines", default="K562,HCT116")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--match-nt", type=int, default=50)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    lines = [x for x in a.lines.split(",") if x]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    rng = np.random.default_rng(a.seed)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    fmt = lambda v: f"{v:.3g}"

    data = {}
    print("== Part 1: the host columns ==")
    for L in lines:
        d, n0 = load(a.pairs_dir, L)
        data[L] = d
        intr = d[d["intron"]]
        print(f"\n[{L}] {n0:,} rows, {len(d):,} scored; classes: " +
              ", ".join(f"{k} {v:,}" for k, v in d["klasa"].value_counts().items()))
        print(f"  intronic host_signal: {describe(intr['host_signal'])}")
        print(f"  intronic host_tx:     {describe(intr['host_tx'])}")
        hs = pd.to_numeric(intr["host_signal"], errors="coerce")
        ht = pd.to_numeric(intr["host_tx"], errors="coerce")
        ok = hs.notna() & ht.notna()
        if ok.sum() > 10:
            print(f"  Spearman(host_signal, host_tx) = {spearmanr(hs[ok], ht[ok])[0]:.3f}")
        act = hs >= ACTIVE
        print(f"  k36 True: {int(intr['k36b'].sum()):,} of {len(intr):,} intronic pairs; "
              f"among active (host_signal >= {ACTIVE:g}) {100 * intr.loc[act, 'k36b'].mean():.1f}%, "
              f"among silent {100 * intr.loc[~act, 'k36b'].mean():.1f}%")
        print(f"  distinct host genes: {intr['host'].nunique():,}; pairs per gene median "
              f"{intr.groupby('host').size().median():.0f}, max {intr.groupby('host').size().max()}")

    print("\n== Part 2: step 07b/07c, active vs silent hosts (pred_DI; observed and sequence-expected mean DI) ==")
    for L in lines:
        f = os.path.join(a.pairs_dir, f"calib_rc_{L}.tsv")
        if not os.path.isfile(f):
            print(f"  [{L}] no {f}")
            continue
        c = pd.read_csv(f, sep="\t")
        c = c[c["metric"] == "pred_DI"]
        cols = ["subset", "group", "n_case", "obs_mean_shift", "exp_shift_linear", "frac_linear", "frac_linear_lo",
                "frac_linear_hi"]
        print(f"  [{L}]")
        print(c[cols].to_string(index=False, float_format=fmt))

    print("\n== Part 3: within each line, residual (observed - sequence-expected) DI by host activity ==")
    print("   DI < 0: sense weaker; residual < 0: sense weaker than its sequence predicts")
    G, C = [], []
    for L in lines:
        g, cc, qtxt = within_line(data[L], L, a.n_boot, rng)
        print(f"\n  [{L}] active quartile edges of host_signal: {qtxt}")
        show = ["group", "n_pairs", "n_genes", "host_signal_median", "obs_mean_DI", "pct_sense_weaker_median",
                "seq_expected_DI", "residual_DI", "residual_lo95", "residual_hi95", "seq_share"]
        print(g[show].to_string(index=False, float_format=fmt))
        if len(cc):
            print(cc.drop(columns=["line"]).to_string(index=False, float_format=fmt))
        G.append(g)
        C.append(cc)
    pd.concat(G).to_csv(a.out_prefix + "_groups.tsv", sep="\t", index=False)
    pd.concat(C).to_csv(a.out_prefix + "_contrasts.tsv", sep="\t", index=False)

    if len(lines) >= 2:
        la, lb = lines[0], lines[1]
        print(f"\n== Part 4: the same element in {la} and {lb} (sequence identical; matched within {a.match_nt} nt, same host) ==")
        M = match_lines(data[la], data[lb], a.match_nt)
        if M is None or len(M) < 20:
            print("  too few matched elements")
        else:
            na, nb = int(data[la]["intron"].sum()), int(data[lb]["intron"].sum())
            print(f"  matched {len(M):,} elements ({100 * len(M) / na:.0f}% of {la}, {100 * len(M) / nb:.0f}% of {lb} "
                  f"intronic pairs); {M['host_a'].nunique():,} host genes")
            X = cross_line(M, la, lb, data[la], data[lb], a.n_boot, rng)
            print(X.to_string(index=False, float_format=fmt))
            X.to_csv(a.out_prefix + "_crossline.tsv", sep="\t", index=False)
            keep = ["chrom_a", "host_a", "host_name_a", "pos_plus_a", "pos_minus_a", "pos_plus_b", "pos_minus_b",
                    "host_signal_a", "host_signal_b", "k36_a", "k36_b", OBS + "_a", OBS + "_b", PRED + "_a", PRED + "_b"]
            M[[k for k in keep if k in M.columns]].to_csv(a.out_prefix + "_crossline_pairs.tsv", sep="\t", index=False)

    print("\nCaveat: pairs come from peak calling; a sense arm must stand out against the host's own sense-strand "
          "signal, which is higher in active genes, so detection may depend on host activity.")
    print(f"CIs: bootstrap over host genes ({a.n_boot} resamples, seed {a.seed}); beta from intergenic pairs "
          f"(both orientations), resampled by pair.")
    print(f"written: {a.out_prefix}_groups.tsv, _contrasts.tsv" + (", _crossline.tsv, _crossline_pairs.tsv" if len(lines) >= 2 else ""))


if __name__ == "__main__":
    main()
