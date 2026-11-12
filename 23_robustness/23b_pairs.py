#!/usr/bin/env python3
"""Step 23b (pairs): the depth-matched asymmetry with confidence intervals that resample host genes,
a paired test, the asymmetry by host activity with a trend test, and the H3K36me3 split.

Pairs in one host gene are not independent observations (K562: 1,611 active-host pairs in 964 genes),
but step 02c resampled single pairs. This script reuses 02c's own functions (imported unchanged,
checked by md5) and:

  A. Replays 02c's main loop with the same random generator (seed, line order, number of resamples),
     so that the matched pairs are exactly those behind istart_depth_matched.tsv; with --check-matched
     the median difference and the number of matched pairs must equal that table, or the script stops.
     For the matched pairs it reports the percent by which the sense start is weaker, 100 (1 - 2^d),
     d = median DI(intronic) - median DI(matched intergenic), with a 95% CI from resampling host genes
     (each intronic pair travels with its matched intergenic pair), and a two-sided Wilcoxon
     signed-rank test on the within-pair differences. Within the five depth bands of 02c it recomputes
     the CIs resampling host genes (intronic) and 1-Mb blocks (intergenic, chromosome + floor(mid/1e6)).
  B. Matches every intronic pair, in active and in silent hosts, to an unused intergenic pair of the
     same depth (02c's match_on_depth, caliper and seed as given), and summarises the asymmetry in
     groups of host activity: silent (host_signal = 0), weak (0 < host_signal < 5) and thirds of
     active hosts (host_signal >= 5, by rank). Trend test: Spearman correlation between host_signal
     and the within-pair difference DI(intronic) - DI(matched intergenic) over all intronic pairs,
     with a 95% CI from resampling host genes; also the difference active minus silent.
  C. In the active-host matched pairs of A: pairs with and without an H3K36me3 peak within 1 kb (k36),
     each with its CI, and their difference (CI from resampling host genes).

Usage:
  python3 23b_pairs.py --step02c 02_istart_baseline/02c_depth_matched.py \
      --table K562=.../istart_directionality_K562.tsv --table HCT116=... --outdir results/23b \
      [--tag all] [--check-matched 02_istart_baseline/results/istart_depth_matched.tsv] \
      [--caliper 0.1] [--n-boot 10000] [--seed 1]
"""
import argparse
import hashlib
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

STEP02C_MD5 = "736af775e86f18a63bde0a578a81e5b0"


def load_module(name, path, md5):
    h = hashlib.md5(open(path, "rb").read()).hexdigest()
    if h != md5:
        sys.exit(f"STOP: {path} has md5 {h}, expected {md5}")
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class Clusters:
    """Resample whole units (host genes or 1-Mb blocks); returns row indices of one resample."""

    def __init__(self, units):
        u, inv = np.unique(np.asarray(units, str), return_inverse=True)
        self.n = len(u)
        self.order = np.argsort(inv, kind="stable")
        self.size = np.bincount(inv, minlength=self.n)
        self.start = np.r_[0, np.cumsum(self.size)[:-1]]

    def draw(self, rng):
        p = rng.integers(0, self.n, self.n)
        ln = self.size[p]
        tot = int(ln.sum())
        idx = np.arange(tot) - np.repeat(np.cumsum(ln) - ln, ln) + np.repeat(self.start[p], ln)
        return self.order[idx]


def pct(d):
    return 100 * (1 - 2 ** d)


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (np.nan, np.nan)


def block_units(df):
    return (df["chrom"].astype(str) + ":" + (df["mid"].astype(np.int64) // 1_000_000).astype(str)).to_numpy()


def matched_summary(x, y, host, n_boot, rng):
    """x, y: DI of intronic pairs and of their matched intergenic pairs (same order); host: unit of x."""
    d = float(np.median(x) - np.median(y))
    C = Clusters(host)
    b = np.empty(n_boot)
    for i in range(n_boot):
        k = C.draw(rng)
        b[i] = np.median(x[k]) - np.median(y[k])
    lo, hi = ci(b)
    w = wilcoxon(x - y, alternative="two-sided") if np.any(x != y) else None
    return {"n_pairs": len(x), "n_host_genes": C.n, "median_DI_intron": float(np.median(x)),
            "median_DI_intergenic": float(np.median(y)), "median_diff": d, "diff_lo95": lo, "diff_hi95": hi,
            "sense_weaker_pct": pct(d), "pct_lo95": pct(hi), "pct_hi95": pct(lo),
            "wilcoxon_signed_rank_p": float(w.pvalue) if w is not None else np.nan,
            "median_within_pair_diff": float(np.median(x - y)),
            "mannwhitney_p": float(mannwhitneyu(x, y, alternative="two-sided").pvalue)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step02c", required=True)
    ap.add_argument("--table", action="append", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default="all")
    ap.add_argument("--check-matched", default="")
    ap.add_argument("--caliper", type=float, default=0.1)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    s02 = load_module("s02c", a.step02c, STEP02C_MD5)
    os.makedirs(a.outdir, exist_ok=True)
    INTRON, ACT = s02.INTRON_CLASS, s02.ACTIVE_MIN_HOST_SIGNAL

    # ---------------- A. replay of 02c, then host-gene CIs
    rng02 = np.random.default_rng(a.seed)
    rngB = np.random.default_rng(a.seed + 1000)
    A_rows, S_rows, P_rows, B_rows, T_rows, K_rows = [], [], [], [], [], []
    data = {}
    for spec in a.table:
        line, path = spec.split("=", 1)
        d = pd.read_csv(path, sep="\t")
        if list(d.columns) != s02.EXPECTED_COLUMNS:
            sys.exit(f"{path}: unexpected columns {list(d.columns)}")
        d["depth"] = d["sens"] + d["anty"]
        gi = d[d["klasa"] == INTRON]
        act = gi[gi["host_signal"] >= ACT].copy()
        inter = d[d["klasa"] != INTRON].copy()
        data[line] = (d, gi, act, inter)
        # strata exactly as 02c (consumes rng02 in the same order)
        pooled = pd.concat([act.assign(group="intron_active"), inter.assign(group="intergenic")])
        edges = np.unique(np.quantile(pooled["depth"], np.linspace(0, 1, s02.N_STRATA + 1)))
        pooled["stratum"] = pd.cut(pooled["depth"], edges, include_lowest=True, labels=False)
        for s, sub in pooled.groupby("stratum"):
            xs = sub[sub["group"] == "intron_active"]
            ys = sub[sub["group"] == "intergenic"]
            if len(xs) < 5 or len(ys) < 5:
                continue
            s02.boot_diff_ci(xs["DI"], ys["DI"], a.n_boot, rng02)       # replay only
            Cx, Cy = Clusters(xs["host"].astype(str)), Clusters(block_units(ys))
            xv, yv = xs["DI"].to_numpy(float), ys["DI"].to_numpy(float)
            b = np.empty(a.n_boot)
            for i in range(a.n_boot):
                b[i] = np.median(xv[Cx.draw(rngB)]) - np.median(yv[Cy.draw(rngB)])
            lo, hi = ci(b)
            dd = float(np.median(xv) - np.median(yv))
            S_rows.append({"line": line, "tag": a.tag, "stratum": int(s), "depth_lo": edges[int(s)],
                           "depth_hi": edges[int(s) + 1], "n_intron_active": len(xs), "n_host_genes": Cx.n,
                           "n_intergenic": len(ys), "n_blocks": Cy.n, "median_diff": dd, "diff_lo95": lo,
                           "diff_hi95": hi, "sense_weaker_pct": pct(dd), "pct_lo95": pct(hi), "pct_hi95": pct(lo)})
        ci_, ti_ = s02.match_on_depth(act["depth"].to_numpy(float), inter["depth"].to_numpy(float), a.caliper, rng02)
        x = act["DI"].to_numpy(float)[ci_]
        y = inter["DI"].to_numpy(float)[ti_]
        s02.boot_diff_ci(x, y, a.n_boot, rng02)                          # replay only
        r = {"line": line, "tag": a.tag, "group": "intron_active", "n_intron_active": len(act), "n_pairs": len(ci_)}
        if len(ci_) >= 5:
            r.update(matched_summary(x, y, act["host"].astype(str).to_numpy()[ci_], a.n_boot, rngB))
        A_rows.append(r)
        mp = pd.DataFrame({"line": line, "chrom": act["chrom"].to_numpy()[ci_], "mid": act["mid"].to_numpy()[ci_],
                           "host": act["host"].to_numpy()[ci_], "host_name": act["host_name"].to_numpy()[ci_],
                           "host_signal": act["host_signal"].to_numpy()[ci_], "k36": act["k36"].to_numpy()[ci_],
                           "depth": act["depth"].to_numpy()[ci_], "DI": x,
                           "ctrl_chrom": inter["chrom"].to_numpy()[ti_], "ctrl_mid": inter["mid"].to_numpy()[ti_],
                           "ctrl_depth": inter["depth"].to_numpy()[ti_], "ctrl_DI": y})
        P_rows.append(mp)
        # C. H3K36me3 split within these matched pairs
        k36 = mp["k36"].map({True: True, False: False, "True": True, "False": False})
        hosts = mp["host"].astype(str).to_numpy()
        sub = {}
        for lab, m in (("k36", (k36 == True).to_numpy()), ("no_k36", (k36 == False).to_numpy())):
            if m.sum() >= 5:
                rr = {"line": line, "tag": a.tag, "set": lab}
                rr.update(matched_summary(x[m], y[m], hosts[m], a.n_boot, rngB))
                K_rows.append(rr)
                sub[lab] = m
        if len(sub) == 2:
            C = Clusters(hosts)
            b = np.empty(a.n_boot)
            m1, m0 = sub["k36"], sub["no_k36"]
            for i in range(a.n_boot):
                k = C.draw(rngB)
                k1, k0 = k[m1[k]], k[m0[k]]
                b[i] = (np.median(x[k1]) - np.median(y[k1])) - (np.median(x[k0]) - np.median(y[k0])) \
                    if len(k1) and len(k0) else np.nan
            lo, hi = ci(b)
            dd = (np.median(x[m1]) - np.median(y[m1])) - (np.median(x[m0]) - np.median(y[m0]))
            K_rows.append({"line": line, "tag": a.tag, "set": "k36 minus no_k36", "median_diff": float(dd),
                           "diff_lo95": lo, "diff_hi95": hi})

    A = pd.DataFrame(A_rows)
    if a.check_matched:
        ref = pd.read_csv(a.check_matched, sep="\t")
        for _, r in A.iterrows():
            q = ref[ref["line"] == r["line"]].iloc[0]
            ok = (int(q["n_matched"]) == int(r["n_pairs"])) and abs(float(q["median_diff"]) - r["median_diff"]) < 5e-4
            print(f"[{r['line']}] replay of 02c: n_matched {r['n_pairs']} vs {int(q['n_matched'])}, median_diff "
                  f"{r['median_diff']:.4f} vs {float(q['median_diff']):.4f}: {'identical' if ok else 'DIFFERENT'}")
            if not ok:
                sys.exit("STOP: the replay did not reproduce 02c's matched pairs")

    # ---------------- B. all intronic pairs matched, by host activity
    rngM = np.random.default_rng(a.seed)
    for line, (d, gi, act, inter) in data.items():
        ci_, ti_ = s02.match_on_depth(gi["depth"].to_numpy(float), inter["depth"].to_numpy(float), a.caliper, rngM)
        x = gi["DI"].to_numpy(float)[ci_]
        y = inter["DI"].to_numpy(float)[ti_]
        hs = gi["host_signal"].to_numpy(float)[ci_]
        hosts = gi["host"].astype(str).to_numpy()[ci_]
        grp = np.full(len(x), "", dtype=object)
        grp[hs == 0] = "silent"
        grp[(hs > 0) & (hs < ACT)] = "weak"
        am = hs >= ACT
        if am.sum() >= 3:
            third = pd.qcut(pd.Series(hs[am]).rank(method="first"), 3, labels=["active_T1", "active_T2", "active_T3"])
            grp[am] = third.astype(str).to_numpy()
        for g in ("silent", "weak", "active_T1", "active_T2", "active_T3"):
            m = grp == g
            rr = {"line": line, "tag": a.tag, "group": g, "n_intronic": len(gi), "n_matched_all": len(ci_),
                  "n_pairs": int(m.sum())}
            if m.sum() >= 5:
                rr.update({"host_signal_min": float(hs[m].min()), "host_signal_median": float(np.median(hs[m])),
                           "host_signal_max": float(hs[m].max())})
                rr.update(matched_summary(x[m], y[m], hosts[m], a.n_boot, rngB))
            B_rows.append(rr)
        delta = x - y
        if len(x) < 10:
            continue
        rho = float(spearmanr(hs, delta)[0])
        p = float(spearmanr(hs, delta)[1])
        C = Clusters(hosts)
        b = np.empty(a.n_boot)
        bd = np.empty(a.n_boot)
        sil, actm = grp == "silent", am
        for i in range(a.n_boot):
            k = C.draw(rngB)
            b[i] = spearmanr(hs[k], delta[k])[0]
            ks, ka = k[sil[k]], k[actm[k]]
            bd[i] = ((np.median(x[ka]) - np.median(y[ka])) - (np.median(x[ks]) - np.median(y[ks]))
                     if len(ks) and len(ka) else np.nan)
        lo, hi = ci(b)
        dlo, dhi = ci(bd)
        dact = float((np.median(x[actm]) - np.median(y[actm])) - (np.median(x[sil]) - np.median(y[sil])))
        T_rows.append({"line": line, "tag": a.tag, "n_pairs": len(x), "n_host_genes": C.n,
                       "spearman_rho_hostsignal_vs_withinpair_diff": rho, "rho_lo95": lo, "rho_hi95": hi,
                       "spearman_p_pairs": p, "active_minus_silent_median_diff": dact,
                       "active_minus_silent_lo95": dlo, "active_minus_silent_hi95": dhi})

    S, B, T, K = pd.DataFrame(S_rows), pd.DataFrame(B_rows), pd.DataFrame(T_rows), pd.DataFrame(K_rows)
    P = pd.concat(P_rows)
    note = f"95% CI: percentile bootstrap, {a.n_boot} resamples (seed {a.seed + 1000}), host genes (intergenic in strata: 1-Mb blocks)"
    for df, name in ((A, "matched"), (S, "strata"), (B, "host_activity"), (T, "trend"), (K, "k36")):
        df["ci"] = note
        df.to_csv(os.path.join(a.outdir, f"23b_{a.tag}_{name}.tsv"), sep="\t", index=False, float_format="%.6g")
    P.to_csv(os.path.join(a.outdir, f"23b_{a.tag}_matched_pairs.tsv"), sep="\t", index=False)
    f = lambda v: f"{v:.3g}"
    print(f"\n== A. depth-matched comparison, intronic pairs in active hosts [{a.tag}]; CIs resample host genes ==")
    A = A.reindex(columns=list(dict.fromkeys(list(A.columns) + ["median_diff", "diff_lo95", "diff_hi95", "sense_weaker_pct",
                                                                  "pct_lo95", "pct_hi95", "wilcoxon_signed_rank_p",
                                                                  "mannwhitney_p", "n_host_genes"])))
    B = B.reindex(columns=list(dict.fromkeys(list(B.columns) + ["n_host_genes", "host_signal_min", "host_signal_max",
                                                                  "sense_weaker_pct", "pct_lo95", "pct_hi95",
                                                                  "wilcoxon_signed_rank_p"])))
    print(A[["line", "n_pairs", "n_host_genes", "median_diff", "diff_lo95", "diff_hi95", "sense_weaker_pct", "pct_lo95",
             "pct_hi95", "wilcoxon_signed_rank_p", "mannwhitney_p"]].to_string(index=False, float_format=f))
    print(f"\n== A. within five depth bands (as 02c); CIs resample host genes and 1-Mb intergenic blocks [{a.tag}] ==")
    print(S[["line", "stratum", "depth_lo", "depth_hi", "n_intron_active", "n_host_genes", "n_intergenic", "sense_weaker_pct",
             "pct_lo95", "pct_hi95"]].to_string(index=False, float_format=f))
    print(f"\n== B. all intronic pairs matched on depth, by host activity [{a.tag}] ==")
    print(B[["line", "group", "n_pairs", "n_host_genes", "host_signal_min", "host_signal_max", "sense_weaker_pct", "pct_lo95",
             "pct_hi95", "wilcoxon_signed_rank_p"]].to_string(index=False, float_format=f))
    print(f"\n== B. trend: Spearman(host_signal, within-pair DI difference); active minus silent [{a.tag}] ==")
    print(T.drop(columns=["ci", "tag"], errors="ignore").to_string(index=False, float_format=f))
    print(f"\n== C. H3K36me3 within 1 kb, matched active-host pairs of A [{a.tag}] ==")
    cols = [c for c in ["line", "set", "n_pairs", "n_host_genes", "median_diff", "diff_lo95", "diff_hi95", "sense_weaker_pct",
                        "pct_lo95", "pct_hi95"] if c in K.columns]
    print(K[cols].to_string(index=False, float_format=f))
    print(f"\n{note}; tests two-sided (Wilcoxon signed-rank on within-pair differences; Mann-Whitney U as 02c).")
    print(f"written: {a.outdir}/23b_{a.tag}_*.tsv")


if __name__ == "__main__":
    main()
