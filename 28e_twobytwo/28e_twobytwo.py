#!/usr/bin/env python3
"""Step 28e (EXPLORATORY; specified on 10 Nov 2026 after step 28d, before running): near the gene's start versus near
a CpG island, for the pair bias (Results 1) and the deficit of strong sites (Results 4), per cell line and pooled.

Cells: start = < 5 kb from the host gene's 5' end (else 'far'); island = <= 2 kb from a CpG island (else 'away').
Pairs: bias 100 x (1 - 2^d), d = median DI of the intronic pairs - median DI of their intergenic partners (as in
step 23b), pooled over lines by pooling the pairs. Sites: deficit 100 x (1 - sense / antisense strong sites) in active
protein-coding genes, pooled by summing over lines. CIs: 2,000 bootstrap resamples of genes drawn jointly for both
lines (gene identifiers without version; a gene present in both lines has the same weight in both); for sites, all
active genes of step 28 are resampled, also those without sites. Readouts:
  R1 / R3  (near start, away from islands) minus (far, away)     - pairs / sites
  R2 / R4  (far from start, near islands) minus (far, away)      - pairs / sites
Two-sided bootstrap P, no correction. Check: per-line estimates of the cells must equal those computed from the
step-28d strata where both exist (near start, away = '> 2 kb from an island' & '<5 kb').

Usage: python3 28e_twobytwo.py <28d results dir> <28 results dir> <step-17 script> <out prefix>
Writes <prefix>_cells.tsv (every cell and contrast) and prints a log.
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

N_BOOT = 2000
CELLS = [("near", "away"), ("near", "island"), ("far", "away"), ("far", "island")]


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
    P = pd.read_csv(f"{r28d}/28d_pairs.tsv", sep="\t", dtype={"host": str})
    S = pd.read_csv(f"{r28d}/28d_sites.tsv.gz", sep="\t", dtype={"gene_id": str})
    S = S[S["metric"] == "strong"]
    lines = ["K562", "HCT116"]
    act = {}
    for ln in lines:
        G = pd.read_csv(f"{r28}/28_{ln}_genes.tsv", sep="\t", dtype={"gene_id": str})
        act[ln] = set(G.loc[G["activity"].isin(["A1", "A2", "A3"]), "gene_id"].str.split(".").str[0])
    P["g"] = P["host"].str.split(".").str[0]
    S["g"] = S["gene_id"].str.split(".").str[0]
    S = S[[g in act[ln] for g, ln in zip(S["g"], S["line"])]]
    units = sorted(set(P["g"]) | act["K562"] | act["HCT116"])
    ui = {g: i for i, g in enumerate(units)}
    W = np.stack([np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)) for _ in range(N_BOOT)])
    for D in (P, S):
        D["u"] = D["g"].map(ui).to_numpy()
        D["start"] = np.where(D["dist5_kb"] < 5, "near", "far")
        D["isl"] = np.where(D["cgi_dist_kb"] <= 2, "island", "away")
    print(f"== 28e (EXPLORATORY): near the gene's start x near a CpG island; {len(P):,} pairs, {len(S):,} strong sites; "
          f"{len(units):,} genes resampled jointly for both lines; bootstrap {N_BOOT}")
    rows = []

    def pair_cell(m):
        DI, C, u = P["DI"].to_numpy(float)[m], P["ctrl_DI"].to_numpy(float)[m], P["u"].to_numpy()[m]
        pt = 100 * (1 - 2 ** (np.median(DI) - np.median(C)))
        bt = np.array([100 * (1 - 2 ** (wmedian(DI, w[u]) - wmedian(C, w[u]))) for w in W])
        return pt, bt

    def site_cell(m):
        u, sen = S["u"].to_numpy()[m], (S["strand"] == "sense").to_numpy()[m]
        cs = np.bincount(u[sen], minlength=len(units)).astype(float)
        ca = np.bincount(u[~sen], minlength=len(units)).astype(float)
        return float(s17.deficit(cs.sum(), ca.sum())), s17.deficit(W @ cs, W @ ca), int(cs.sum()), int(ca.sum())

    for what, D, fn in (("pairs", P, pair_cell), ("sites", S, site_cell)):
        print(f"\n-- {what}: {'bias' if what == 'pairs' else 'deficit of strong sites'}, %, [95% CI]")
        for scope in lines + ["pooled"]:
            inl = np.ones(len(D), bool) if scope == "pooled" else (D["line"] == scope).to_numpy()
            res = {}
            for st, isl in CELLS:
                m = inl & (D["start"] == st).to_numpy() & (D["isl"] == isl).to_numpy()
                r = fn(m)
                pt, bt = r[0], r[1]
                lo, hi = s17.ci(bt)
                res[(st, isl)] = (pt, bt)
                extra = f"  {r[2]:,} / {r[3]:,} sites" if what == "sites" else f"  n = {int(m.sum()):,}"
                rows.append({"what": what, "scope": scope, "cell": f"{st} start, {isl}", "n": int(m.sum()), "estimate": pt,
                             "lo": lo, "hi": hi, "p_boot": np.nan})
                print(f"  {scope:7s} {st:4s} start, {isl:6s} {pt:+6.1f} [{lo:+.1f}, {hi:+.1f}]{extra}")
            for lab, a_, b_ in ((f"R{1 if what == 'pairs' else 3}: near start & away  minus  far & away", ("near", "away"), ("far", "away")),
                                (f"R{2 if what == 'pairs' else 4}: far & island  minus  far & away", ("far", "island"), ("far", "away"))):
                v = res[a_][1] - res[b_][1]
                lo, hi = s17.ci(v)
                p = s17.pval(v)
                rows.append({"what": what, "scope": scope, "cell": lab, "n": int(inl.sum()), "estimate": res[a_][0] - res[b_][0],
                             "lo": lo, "hi": hi, "p_boot": p})
                print(f"  {scope:7s} {lab}: {res[a_][0] - res[b_][0]:+.1f} [{lo:+.1f}, {hi:+.1f}], p = {p:.2g}")

    # check against step 28d (point estimates of 'near start, away from islands' per line)
    st28d = pd.read_csv(f"{r28d}/28d_strata.tsv", sep="\t")
    R = pd.DataFrame(rows)
    for what, an in (("pairs", "pairs"), ("sites", "sites_strong")):
        for ln in lines:
            ref = st28d[(st28d["line"] == ln) & (st28d["analysis"] == an) & (st28d["subset"] == "> 2 kb from an island")
                        & (st28d["region"] == "<5 kb")]["estimate"].iloc[0]
            mine = R[(R["what"] == what) & (R["scope"] == ln) & (R["cell"] == "near start, away")]["estimate"].iloc[0]
            if not np.isclose(ref, mine):
                sys.exit(f"STOP: {what} {ln}: {mine} differs from step 28d {ref}")
    print("\ncheck: per-line 'near start, away from islands' estimates equal step 28d")
    R.to_csv(out + "_cells.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"written: {out}_cells.tsv")


if __name__ == "__main__":
    main()
