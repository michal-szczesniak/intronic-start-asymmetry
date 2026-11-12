"""Independent loop-based check of step 28d: distances to islands, site counts and pair estimates for the subset
'> 2 kb from an island'. Usage: python3 naive_check_28d.py cgi.txt.gz <28b results> <28 results> <28d prefix> LINE"""
import gzip, sys, numpy as np, pandas as pd
cgi_p, r28b, r28, pre, line = sys.argv[1:6]
isl = {}
for l in gzip.open(cgi_p, "rt"):
    f = l.split("\t"); isl.setdefault(f[1], []).append((int(f[2]), int(f[3])))
def dist(c, p):
    best = float("inf")
    for s, e in isl.get(c, []):
        d = 0 if s <= p < e else (s - p if p < s else p - (e - 1))
        best = min(best, d)
    return best
G = pd.read_csv(f"{r28}/28_{line}_genes.tsv", sep="\t")
act = set(G.loc[G.activity.isin(["A1", "A2", "A3"]), "gene_id"])
S = pd.read_csv(f"{r28b}/28b_{line}_sites.tsv.gz", sep="\t")
S = S[(S.metric == "strong") & S.gene_id.isin(act)]
R = pd.read_csv(pre + "_strata.tsv", sep="\t")
bad = n = 0
for reg, cond in (("<5 kb", lambda d5: d5 < 5), (">=5 kb", lambda d5: d5 >= 5)):
    cs = ca = 0
    for r in S.itertuples():
        if cond(r.dist5_kb) and dist(r.chrom, r.pos) / 1e3 > 2:
            cs += r.strand == "sense"; ca += r.strand == "anti"
    row = R[(R.line == line) & (R.analysis == "sites_strong") & (R.subset == "> 2 kb from an island") & (R.region == reg)].iloc[0]
    n += 2; bad += (row.sense != cs) + (row.anti != ca)
    if (row.sense != cs) or (row.anti != ca): print("DIFF sites", reg, row.sense, row.anti, cs, ca)
P = pd.read_csv(f"{r28b}/28b_{line}_pairs.tsv", sep="\t")
for reg, cond in (("<5 kb", lambda d5: d5 < 5), (">=5 kb", lambda d5: d5 >= 5)):
    m = [cond(r.dist5_kb) and dist(r.chrom, r.mid) / 1e3 > 2 for r in P.itertuples()]
    d = P[m]
    row = R[(R.line == line) & (R.analysis == "pairs") & (R.subset == "> 2 kb from an island") & (R.region == reg)]
    if len(d) < 20:
        continue
    exp = 100 * (1 - 2 ** (np.median(d.DI) - np.median(d.ctrl_DI)))
    n += 1
    if not np.isclose(row.estimate.iloc[0], exp):
        bad += 1; print("DIFF pairs", reg, row.estimate.iloc[0], exp)
print(f"compared {n} values, {bad} differ")
