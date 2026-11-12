#!/usr/bin/env python3
"""Step 28c (EXPLORATORY, 10 Nov 2026, after step 28b): bias and deficit beyond the first 5 kb of the gene.

Step 28b placed most of the strand difference within 1-5 kb of the gene's 5' end. Step 28b printed the estimate
for the rest of the gene without a confidence interval; this script adds it, from the step-28b output tables only:
  pairs   bias 100 x (1 - 2^d) (d = median DI of intronic pairs - median DI of their intergenic partners) for pairs
          < 5 kb and >= 5 kb from the 5' end of the host gene; CIs and two-sided bootstrap P of the difference from
          2,000 resamples of host genes (each pair with its partner)
  sites   deficit of strong sites, 100 x (1 - sense / antisense), in active protein-coding genes, < 5 kb and >= 5 kb
          from the gene's 5' end and from the nearest upstream annotated TSS; 2,000 resamples of genes within
          activity thirds (all active genes of step 28, also those without sites)
Check: the point estimates for < 5 kb must equal step 28b's 1-5 kb rows (eligible positions start at about 1 kb).
Usage: python3 28c_beyond5kb.py <step-17 script> <results dir of 28b> <results dir of 28> <out.tsv>
"""
import importlib.util
import sys

import numpy as np
import pandas as pd

s17p, r28b, r28, out = sys.argv[1:5]
spec = importlib.util.spec_from_file_location("s17", s17p)
s17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s17)
N_BOOT, rows = 2000, []


def wmedian(v, w):
    o = np.argsort(v, kind="stable")
    v, c = v[o], np.cumsum(w[o])
    N = int(c[-1])
    if N == 0:
        return np.nan
    if N % 2:
        return float(v[np.searchsorted(c, N // 2, side="right")])
    return 0.5 * float(v[np.searchsorted(c, N // 2 - 1, side="right")] + v[np.searchsorted(c, N // 2, side="right")])


for line in ("K562", "HCT116"):
    rng = np.random.default_rng(1)
    S28b = pd.read_csv(f"{r28b}/28b_{line}_strata.tsv", sep="\t")
    # ---- pairs
    P = pd.read_csv(f"{r28b}/28b_{line}_pairs.tsv", sep="\t", dtype={"host": str})
    hosts, hidx = np.unique(P["host"].to_numpy(), return_inverse=True)
    W = np.stack([np.bincount(rng.integers(0, len(hosts), len(hosts)), minlength=len(hosts)) for _ in range(N_BOOT)])
    DI, C = P["DI"].to_numpy(float), P["ctrl_DI"].to_numpy(float)
    near = (P["dist5_kb"] < 5).to_numpy()
    res = {}
    for nm, m in (("<5 kb", near), (">=5 kb", ~near)):
        pt = 100 * (1 - 2 ** (np.median(DI[m]) - np.median(C[m])))
        bt = np.array([100 * (1 - 2 ** (wmedian(DI[m], w[hidx[m]]) - wmedian(C[m], w[hidx[m]]))) for w in W])
        res[nm] = (pt, bt)
        lo, hi = s17.ci(bt)
        rows.append({"line": line, "analysis": "pairs_dist5", "stratum": nm, "n": int(m.sum()),
                     "units": int(len(np.unique(hidx[m]))), "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan})
    ref = S28b[(S28b["analysis"] == "pairs_dist_bin") & (S28b["stratum"] == "1-5")]["estimate"].iloc[0]
    if not np.isclose(res["<5 kb"][0], ref):
        sys.exit(f"STOP: {line} pairs < 5 kb {res['<5 kb'][0]} differs from step 28b {ref}")
    v = res["<5 kb"][1] - res[">=5 kb"][1]
    lo, hi = s17.ci(v)
    rows.append({"line": line, "analysis": "pairs_dist5", "stratum": "<5 minus >=5 kb", "n": len(P), "units": len(hosts),
                 "estimate": res["<5 kb"][0] - res[">=5 kb"][0], "lo": lo, "hi": hi, "p_boot": s17.pval(v)})
    # ---- sites
    G = pd.read_csv(f"{r28}/28_{line}_genes.tsv", sep="\t", dtype={"gene_id": str})
    G = G[G["activity"].isin(["A1", "A2", "A3"])].reset_index(drop=True)
    gi = {g: i for i, g in enumerate(G["gene_id"])}
    parts = [np.flatnonzero(G["activity"].to_numpy() == t) for t in ("A1", "A2", "A3")]
    Ws = [s17.boot_W(len(ix), N_BOOT, rng) for ix in parts]
    S = pd.read_csv(f"{r28b}/28b_{line}_sites.tsv.gz", sep="\t", dtype={"gene_id": str})
    S = S[(S["metric"] == "strong") & S["gene_id"].isin(gi)]
    g_of = S["gene_id"].map(gi).to_numpy()
    sense = (S["strand"] == "sense").to_numpy()
    for col, ref_an in (("dist5_kb", "sites_strong_bin5"), ("dist_tss_kb", "sites_strong_bintss")):
        res = {}
        for nm, m in (("<5 kb", (S[col] < 5).to_numpy()), (">=5 kb", (S[col] >= 5).to_numpy())):
            cs = np.bincount(g_of[m & sense], minlength=len(G)).astype(float)
            ca = np.bincount(g_of[m & ~sense], minlength=len(G)).astype(float)
            bS = sum(Wt @ cs[ix] for Wt, ix in zip(Ws, parts))
            bA = sum(Wt @ ca[ix] for Wt, ix in zip(Ws, parts))
            pt, dr = float(s17.deficit(cs.sum(), ca.sum())), s17.deficit(bS, bA)
            res[nm] = (pt, dr)
            lo, hi = s17.ci(dr)
            rows.append({"line": line, "analysis": f"sites_strong_{col}", "stratum": nm, "n": int(cs.sum() + ca.sum()),
                         "units": len(G), "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan,
                         "sense": int(cs.sum()), "anti": int(ca.sum())})
        ref = S28b[(S28b["analysis"] == ref_an) & (S28b["stratum"] == "1-5")]["estimate"].iloc[0]
        if not np.isclose(res["<5 kb"][0], ref):
            sys.exit(f"STOP: {line} sites {col} < 5 kb {res['<5 kb'][0]} differs from step 28b {ref}")
        v = res["<5 kb"][1] - res[">=5 kb"][1]
        lo, hi = s17.ci(v)
        rows.append({"line": line, "analysis": f"sites_strong_{col}", "stratum": "<5 minus >=5 kb", "n": len(S),
                     "units": len(G), "estimate": res["<5 kb"][0] - res[">=5 kb"][0], "lo": lo, "hi": hi, "p_boot": s17.pval(v)})
R = pd.DataFrame(rows)
R.to_csv(out, sep="\t", index=False, float_format="%.6g")
pd.set_option("display.width", 200)
print(R.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
