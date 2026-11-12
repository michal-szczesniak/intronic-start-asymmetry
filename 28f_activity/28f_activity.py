#!/usr/bin/env python3
"""Step 28f (EXPLORATORY; specified on 10 Nov 2026 after step 28e, before running): does the near-start effect grow
with the activity of the host gene?

Near-start effect = strand difference < 5 kb from the host gene's 5' end minus that >= 5 kb (pairs: bias
100 x (1 - 2^d), d = median DI of intronic pairs - median DI of their intergenic partners; strong sites: deficit
100 x (1 - sense / antisense) in active protein-coding genes). Activity: pairs - thirds of host_signal (PRO-cap reads
at the host's annotated start, the measure used for pairs) over the host genes of the matched pairs within each line;
sites - the activity thirds A1-A3 of step 28. Per line and pooled over both lines; CIs from 2,000 bootstrap resamples
of genes drawn jointly for both lines (identifiers without version). Readout: interaction = near-start effect in the
most active third minus that in the least active third. Secondary: the same for sites and pairs > 2 kb from any CpG
island. Two-sided bootstrap P, no correction. Check: over all thirds, the near and far estimates of each line must
equal step 28d (all pairs / all sites; STOP otherwise).

Usage: python3 28f_activity.py <28d results dir> <28 results dir> <step-17 script> <out prefix>
Writes <prefix>_cells.tsv and prints a log.
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

N_BOOT = 2000
THIRDS = ["T1", "T2", "T3"]


def wmedian(v, w):
    o = np.argsort(v, kind="stable")
    v, c = v[o], np.cumsum(w[o])
    N = int(c[-1]) if len(c) else 0
    if N == 0:
        return np.nan
    if N % 2:
        return float(v[np.searchsorted(c, N // 2, side="right")])
    return 0.5 * float(v[np.searchsorted(c, N // 2 - 1, side="right")] + v[np.searchsorted(c, N // 2, side="right")])


def main():
    r28d, r28, s17p, out = sys.argv[1:5]
    spec = importlib.util.spec_from_file_location("s17", s17p)
    s17 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s17)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    rng = np.random.default_rng(1)
    lines = ["K562", "HCT116"]
    P = pd.read_csv(f"{r28d}/28d_pairs.tsv", sep="\t", dtype={"host": str})
    S = pd.read_csv(f"{r28d}/28d_sites.tsv.gz", sep="\t", dtype={"gene_id": str})
    S = S[S["metric"] == "strong"]
    for c in ("host_signal",):
        if c not in P.columns:
            sys.exit(f"STOP: 28d_pairs.tsv lacks {c}")
    act = {}
    for ln in lines:
        G = pd.read_csv(f"{r28}/28_{ln}_genes.tsv", sep="\t", dtype={"gene_id": str})
        act[ln] = dict(zip(G["gene_id"].str.split(".").str[0], G["activity"]))
    P["g"] = P["host"].str.split(".").str[0]
    S["g"] = S["gene_id"].str.split(".").str[0]
    S["third"] = [{"A1": "T1", "A2": "T2", "A3": "T3"}.get(act[ln].get(g, ""), "") for g, ln in zip(S["g"], S["line"])]
    S = S[S["third"] != ""].reset_index(drop=True)
    P["third"] = ""
    for ln in lines:
        m = (P["line"] == ln).to_numpy()
        h = P.loc[m].groupby("g")["host_signal"].first()
        e = np.quantile(h.to_numpy(float), [1 / 3, 2 / 3])
        lab = pd.Series(np.array(THIRDS)[np.searchsorted(e, h.to_numpy(float), side="right").clip(0, 2)], index=h.index)
        P.loc[m, "third"] = P.loc[m, "g"].map(lab).to_numpy()
        print(f"{ln}: host_signal thirds over {len(h):,} host genes of the pairs, edges {e[0]:.0f} and {e[1]:.0f} reads; "
              f"pairs per third " + ", ".join(f"{t} {int((P.loc[m, 'third'] == t).sum()):,}" for t in THIRDS))
    units = sorted(set(P["g"]) | {g for ln in lines for g, a in act[ln].items() if a in ("A1", "A2", "A3")})
    ui = {g: i for i, g in enumerate(units)}
    W = np.stack([np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)) for _ in range(N_BOOT)])
    for D in (P, S):
        D["u"] = D["g"].map(ui).to_numpy()
        D["near"] = (D["dist5_kb"] < 5).to_numpy()
        D["away"] = (D["cgi_dist_kb"] > 2).to_numpy()
    print(f"== 28f (EXPLORATORY): near-start effect by host activity; {len(P):,} pairs, {len(S):,} strong sites in active "
          f"genes; {len(units):,} genes resampled jointly; bootstrap {N_BOOT}")

    def pair_est(m):
        if not m.any():
            return np.nan, np.full(N_BOOT, np.nan), 0
        DI, C, u = P["DI"].to_numpy(float)[m], P["ctrl_DI"].to_numpy(float)[m], P["u"].to_numpy()[m]
        pt = 100 * (1 - 2 ** (np.median(DI) - np.median(C)))
        return pt, np.array([100 * (1 - 2 ** (wmedian(DI, w[u]) - wmedian(C, w[u]))) for w in W]), int(m.sum())

    def site_est(m):
        u, sen = S["u"].to_numpy()[m], (S["strand"] == "sense").to_numpy()[m]
        cs = np.bincount(u[sen], minlength=len(units)).astype(float)
        ca = np.bincount(u[~sen], minlength=len(units)).astype(float)
        return float(s17.deficit(cs.sum(), ca.sum())), s17.deficit(W @ cs, W @ ca), int(m.sum())

    rows = []
    st28d = pd.read_csv(f"{r28d}/28d_strata.tsv", sep="\t")
    for what, D, fn, an28d, sub28d in (("pairs", P, pair_est, "pairs", "all pairs"),
                                       ("sites", S, site_est, "sites_strong", "all sites")):
        for subset in ("all", "away from islands"):
            print(f"\n-- {what}, {subset}: near (< 5 kb) and far (>= 5 kb), %, [95% CI]; near-start effect = near - far")
            for scope in lines + ["pooled"]:
                base = np.ones(len(D), bool) if scope == "pooled" else (D["line"] == scope).to_numpy()
                if subset != "all":
                    base = base & D["away"].to_numpy()
                if subset == "all" and scope != "pooled":      # check against step 28d
                    for reg, rm in (("<5 kb", D["near"].to_numpy()), (">=5 kb", ~D["near"].to_numpy())):
                        mine = fn(base & rm)[0]
                        ref = st28d[(st28d["line"] == scope) & (st28d["analysis"] == an28d) & (st28d["subset"] == sub28d)
                                    & (st28d["region"] == reg)]["estimate"].iloc[0]
                        if not np.isclose(mine, ref):
                            sys.exit(f"STOP: {what} {scope} {reg}: {mine} differs from step 28d {ref}")
                eff = {}
                for t in THIRDS:
                    mt = base & (D["third"] == t).to_numpy()
                    a_ = fn(mt & D["near"].to_numpy())
                    b_ = fn(mt & ~D["near"].to_numpy())
                    v = a_[1] - b_[1]
                    lo, hi = s17.ci(v)
                    eff[t] = (a_[0] - b_[0], v)
                    for reg, r in (("near", a_), ("far", b_)):
                        l2, h2 = s17.ci(r[1])
                        rows.append({"what": what, "subset": subset, "scope": scope, "third": t, "cell": reg, "n": r[2],
                                     "estimate": r[0], "lo": l2, "hi": h2, "p_boot": np.nan})
                    rows.append({"what": what, "subset": subset, "scope": scope, "third": t, "cell": "near minus far",
                                 "n": a_[2] + b_[2], "estimate": a_[0] - b_[0], "lo": lo, "hi": hi, "p_boot": s17.pval(v)})
                    print(f"  {scope:7s} {t}: near {a_[0]:+6.1f} (n {a_[2]:,})  far {b_[0]:+6.1f} (n {b_[2]:,})  "
                          f"effect {a_[0] - b_[0]:+6.1f} [{lo:+.1f}, {hi:+.1f}]")
                v = eff["T3"][1] - eff["T1"][1]
                lo, hi = s17.ci(v)
                p = s17.pval(v)
                rows.append({"what": what, "subset": subset, "scope": scope, "third": "T3-T1", "cell": "interaction",
                             "n": int(base.sum()), "estimate": eff["T3"][0] - eff["T1"][0], "lo": lo, "hi": hi, "p_boot": p})
                print(f"  {scope:7s} INTERACTION (effect in T3 - effect in T1): {eff['T3'][0] - eff['T1'][0]:+.1f} "
                      f"[{lo:+.1f}, {hi:+.1f}], p = {p:.2g}")
    print("\ncheck: per-line near and far estimates over all thirds equal step 28d")
    pd.DataFrame(rows).to_csv(out + "_cells.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"written: {out}_cells.tsv")


if __name__ == "__main__":
    main()
