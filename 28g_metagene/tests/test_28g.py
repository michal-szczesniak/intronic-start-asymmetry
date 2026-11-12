#!/usr/bin/env python3
"""Synthetic test of 28g_metagene.py: random genes and sites; counts per bin, genes reaching each bin and the
bootstrap ratios are recomputed by explicit loops (same seed and draw order) and compared; the checks must stop
the script on a wrong expected total."""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

script, work = sys.argv[1], sys.argv[2]
os.makedirs(work, exist_ok=True)
rng = np.random.default_rng(7)
n = 300
G = pd.DataFrame({"chrom": rng.choice(["chr1", "chr2"], n), "gene_id": [f"G{i}" for i in range(n)],
                  "gene_type": rng.choice(["protein_coding", "lncRNA"], n, p=[0.85, 0.15]),
                  "activity": rng.choice(["silent", "low", "A1", "A2", "A3"], n), "start": rng.integers(0, 10**6, n)})
G["end"] = G["start"] + rng.integers(2000, 60000, n)
G.to_csv(f"{work}/genes.tsv", sep="\t", index=False)
m = 2500
gi = rng.integers(0, n, m)
dist = rng.uniform(0.3, 1.0, m) * (G["end"] - G["start"]).to_numpy()[gi] / 1e3
S = pd.DataFrame({"metric": rng.choice(["strong", "r10"], m, p=[0.8, 0.2]), "chrom": G["chrom"].to_numpy()[gi],
                  "gene_id": G["gene_id"].to_numpy()[gi], "activity": G["activity"].to_numpy()[gi],
                  "strand": rng.choice(["sense", "anti"], m), "dist5_kb": dist})
S.loc[:5, "dist5_kb"] = 0.5   # a few sites below 1 kb
S.to_csv(f"{work}/sites.tsv.gz", sep="\t", index=False, compression="gzip")

act = G[(G.gene_type == "protein_coding") & G.activity.isin(["A1", "A2", "A3"])].reset_index(drop=True)
x = S[(S.metric == "strong") & S.activity.isin(["A1", "A2", "A3"])]
x = x[x.gene_id.isin(act.gene_id)]
tot = (int((x.strand == "sense").sum()), int((x.strand == "anti").sum()))
fst = (int(((x.strand == "sense") & (x.dist5_kb < 5)).sum()), int(((x.strand == "anti") & (x.dist5_kb < 5)).sum()))
# sites of non-protein-coding active genes must not exist in the input of the script (as in step 28b)
S2 = S[~(S.activity.isin(["A1", "A2", "A3"]) & ~S.gene_id.isin(act.gene_id))]
S2.to_csv(f"{work}/sites.tsv.gz", sep="\t", index=False, compression="gzip")
x = S2[(S2.metric == "strong") & S2.activity.isin(["A1", "A2", "A3"])]
args = [sys.executable, script, "--label", "SYN", "--sites", f"{work}/sites.tsv.gz", "--genes", f"{work}/genes.tsv",
        "--expect-genes", str(len(act)), "--expect-total", f"{tot[0]},{tot[1]}", "--expect-first", f"{fst[0]},{fst[1]}",
        "--n-boot", "50", "--max-kb", "30", "--out", f"{work}/out.tsv"]
r = subprocess.run(args, capture_output=True, text=True)
print(r.stdout[-400:], r.stderr[-400:])
assert r.returncode == 0
bad = subprocess.run(args[:-6] + ["--expect-first", f"{fst[0] + 1},{fst[1]}", "--n-boot", "50", "--out", f"{work}/x.tsv"],
                     capture_output=True, text=True)
assert bad.returncode != 0 and "STOP" in (bad.stdout + bad.stderr), "the check did not stop the script"
R = pd.read_csv(f"{work}/out.tsv", sep="\t")
# independent recount
key = {g: i for i, g in enumerate(act.gene_id)}
nb = 29
CS = np.zeros((len(act), nb)); CA = np.zeros((len(act), nb))
for g, st, dd in zip(x.gene_id, x.strand, x.dist5_kb):
    if dd >= 30:
        continue
    k = 0 if dd < 1 else int(dd) - 1
    (CS if st == "sense" else CA)[key[g], k] += 1
rng2 = np.random.default_rng(1)
Ws = []
parts = [np.flatnonzero(act.activity.to_numpy() == t) for t in ("A1", "A2", "A3")]
for ix in parts:
    idx = rng2.integers(0, len(ix), size=(50, len(ix)))
    Ws.append(np.array([[np.sum(idx[b] == j) for j in range(len(ix))] for b in range(50)], float))
diff = 0
compared = 0
glen = (act.end - act.start).to_numpy() / 1e3
for k in range(nb):
    s, an = CS[:, k].sum(), CA[:, k].sum()
    reach = int(sum(1 for L in glen if L > k + 1))
    bs = np.array([sum(Ws[t][b] @ CS[parts[t], k] for t in range(3)) for b in range(50)])
    ba = np.array([sum(Ws[t][b] @ CA[parts[t], k] for t in range(3)) for b in range(50)])
    rr = bs / ba
    rr = rr[np.isfinite(rr)]
    lo_, hi_ = max(k - 1, 0), min(k + 1, nb - 1)
    s3, a3 = CS[:, lo_:hi_ + 1].sum(), CA[:, lo_:hi_ + 1].sum()
    exp = {"sense": s, "anti": an, "genes_reaching": reach, "ratio": s / an if an else np.nan,
           "lo": np.percentile(rr, 2.5), "hi": np.percentile(rr, 97.5), "sense3": s3, "anti3": a3}
    for c, v in exp.items():
        compared += 1
        if not np.isclose(R.loc[k, c], v, rtol=1e-5, equal_nan=True):
            diff += 1
            print("differs", k, c, R.loc[k, c], v)
print(f"compared {compared} values, {diff} differ")
