"""Independent check of step 28f point estimates (direct computation). Usage: naive_check_28f.py <28d dir> <28 dir> <28f prefix>"""
import sys, numpy as np, pandas as pd
r28d, r28, pre = sys.argv[1:4]
P = pd.read_csv(f"{r28d}/28d_pairs.tsv", sep="\t"); S = pd.read_csv(f"{r28d}/28d_sites.tsv.gz", sep="\t")
S = S[S.metric == "strong"]
parts = []
for ln in ("K562", "HCT116"):
    G = pd.read_csv(f"{r28}/28_{ln}_genes.tsv", sep="\t")
    a = dict(zip(G.gene_id, G.activity))
    s = S[S.line == ln].copy(); s["third"] = s.gene_id.map(a).map({"A1": "T1", "A2": "T2", "A3": "T3"}); parts.append(s.dropna(subset=["third"]))
S = pd.concat(parts)
for ln in ("K562", "HCT116"):
    h = P[P.line == ln].groupby("host").host_signal.first()
    e1, e2 = np.quantile(h, [1/3, 2/3])
    t = h.apply(lambda v: "T1" if v < e1 else ("T2" if v < e2 else "T3"))
    P.loc[P.line == ln, "third"] = P.loc[P.line == ln, "host"].map(t)
R = pd.read_csv(pre + "_cells.tsv", sep="\t")
def est(what, d):
    if what == "pairs":
        return 100 * (1 - 2 ** (np.median(d.DI) - np.median(d.ctrl_DI)))
    return 100 * (1 - (d.strand == "sense").sum() / (d.strand == "anti").sum())
bad = n = 0
for what, D in (("pairs", P), ("sites", S)):
    for subset in ("all", "away from islands"):
        for scope in ("K562", "HCT116", "pooled"):
            d0 = D if scope == "pooled" else D[D.line == scope]
            if subset != "all":
                d0 = d0[d0.cgi_dist_kb > 2]
            for th in ("T1", "T2", "T3"):
                for cell, cond in (("near", d0.dist5_kb < 5), ("far", d0.dist5_kb >= 5)):
                    d = d0[cond & (d0.third == th)]
                    if len(d) == 0:
                        continue
                    exp = est(what, d)
                    got = R[(R.what == what) & (R.subset == subset) & (R.scope == scope) & (R.third == th) & (R.cell == cell)].estimate.iloc[0]
                    n += 1
                    if not (np.isclose(got, exp, rtol=1e-5) or (np.isnan(got) and np.isnan(exp))):
                        bad += 1; print("DIFF", what, subset, scope, th, cell, got, exp)
print(f"compared {n} values, {bad} differ")
