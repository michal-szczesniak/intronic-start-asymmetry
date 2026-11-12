"""Independent check of step 28e point estimates (direct pandas computation, no bootstrap).
Usage: python3 naive_check_28e.py <28d results dir> <28 results dir> <28e prefix>"""
import sys, numpy as np, pandas as pd
r28d, r28, pre = sys.argv[1:4]
P = pd.read_csv(f"{r28d}/28d_pairs.tsv", sep="\t"); S = pd.read_csv(f"{r28d}/28d_sites.tsv.gz", sep="\t")
S = S[S.metric == "strong"]
keep = []
for ln in ("K562", "HCT116"):
    G = pd.read_csv(f"{r28}/28_{ln}_genes.tsv", sep="\t")
    a = set(G.loc[G.activity.isin(["A1", "A2", "A3"]), "gene_id"])
    keep.append(S[(S.line == ln) & S.gene_id.isin(a)])
S = pd.concat(keep)
R = pd.read_csv(pre + "_cells.tsv", sep="\t")
bad = n = 0
for what, D in (("pairs", P), ("sites", S)):
    for scope in ("K562", "HCT116", "pooled"):
        d = D if scope == "pooled" else D[D.line == scope]
        for st in ("near", "far"):
            for isl in ("away", "island"):
                m = d[((d.dist5_kb < 5) == (st == "near")) & ((d.cgi_dist_kb <= 2) == (isl == "island"))]
                if what == "pairs":
                    exp = 100 * (1 - 2 ** (np.median(m.DI) - np.median(m.ctrl_DI)))
                else:
                    exp = 100 * (1 - (m.strand == "sense").sum() / (m.strand == "anti").sum())
                got = R[(R.what == what) & (R.scope == scope) & (R.cell == f"{st} start, {isl}")].estimate.iloc[0]
                n += 1
                if not np.isclose(got, exp, rtol=1e-5):
                    bad += 1; print("DIFF", what, scope, st, isl, got, exp)
print(f"compared {n} values, {bad} differ")
