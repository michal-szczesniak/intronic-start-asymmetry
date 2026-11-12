#!/usr/bin/env python3
"""Synthetic test of 29b_liftover.py against an independent position-by-position mapper.

The synthetic chain has: a chain with gaps on the same strand; an inverted chain (query on the minus strand); two
chains whose target blocks overlap (positions there must be dropped as ambiguous); a chain to a sequence that is not
kept; a chain that maps a second hg19 region onto the same hg38 bases as the first (reads must be added); gaps and a
hg19 sequence without any chain (unmapped). Random single-base counts on both strands (minus as negative values).
Usage: python3 test_29b_liftover.py <29b_liftover.py> <empty work folder>
Prints 'compared N values, 0 differ' when the outputs agree."""
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pyBigWig

script, wd = sys.argv[1:3]
os.makedirs(wd, exist_ok=True)
rng = np.random.default_rng(11)
T = {"tA": 5000, "tB": 3000, "tC": 1000}                        # hg19 sequences
Q = {"qA": 6000, "qB": 4000, "qZ_alt": 2000}                    # hg38 sequences; qZ_alt is not kept
# chains: (tname, tstart, qname, qstrand, qstart, [(size, dt, dq), ..., (size,)])
CH = [("tA", 100, "qA", "+", 300, [(400, 50, 20), (300, 0, 100), (200,)]),
      ("tA", 1500, "qB", "-", 500, [(600, 30, 10), (250,)]),
      ("tA", 1700, "qA", "+", 3000, [(300,)]),                   # overlaps the inverted chain on tA 1700-2000
      ("tA", 2600, "qZ_alt", "+", 100, [(400,)]),
      ("tA", 3200, "qA", "+", 300, [(150,)]),                    # same hg38 bases as the first block of chain 1
      ("tB", 0, "qB", "+", 3000, [(500, 100, 0), (400,)])]       # tC has no chain
with gzip.open(os.path.join(wd, "syn.chain.gz"), "wt") as fh:
    for k, (tn, ts, qn, qs, qp, bl) in enumerate(CH):
        tlen = sum(b[0] + (b[1] if len(b) == 3 else 0) for b in bl)
        qlen = sum(b[0] + (b[2] if len(b) == 3 else 0) for b in bl)
        fh.write(f"chain 1000 {tn} {T[tn]} + {ts} {ts + tlen} {qn} {Q[qn]} {qs} {qp} {qp + qlen} {k + 1}\n")
        for b in bl:
            fh.write("\t".join(map(str, b)) + "\n")
        fh.write("\n")
with open(os.path.join(wd, "syn.fa.fai"), "w") as fh:
    for q, L in Q.items():
        fh.write(f"{q}\t{L}\t0\t60\t61\n")

data = {}
for strand, sign in (("+", 1), ("-", -1)):
    bw = pyBigWig.open(os.path.join(wd, f"hg19_{'plus' if strand == '+' else 'minus'}.bw"), "w")
    bw.addHeader(list(T.items()))
    for t, L in T.items():
        pos = np.sort(rng.choice(L, L // 3, replace=False))
        val = rng.integers(1, 30, len(pos))
        bw.addEntries([t] * len(pos), pos.tolist(), ends=(pos + 1).tolist(), values=(sign * val).astype(float).tolist())
        data[(t, strand)] = (pos, val)
    bw.close()

# independent mapper: every block of every chain, position by position
blocks = []
for tn, ts, qn, qs, qp, bl in CH:
    for b in bl:
        blocks.append((tn, ts, ts + b[0], qn, qs, qp))
        if len(b) == 3:
            ts += b[0] + b[1]
            qp += b[0] + b[2]
        else:
            ts += b[0]
            qp += b[0]
exp = defaultdict(int)
cnt = defaultdict(int)
for (t, strand), (pos, val) in data.items():
    for p, v in zip(pos.tolist(), val.tolist()):
        hits = [b for b in blocks if b[0] == t and b[1] <= p < b[2]]
        if not hits:
            cnt["unmapped"] += v
            continue
        if len(hits) > 1:
            cnt["ambiguous"] += v
            continue
        _, bs, _, qn, qs, qp = hits[0]
        x = qp + (p - bs)
        if qs == "-":
            x = Q[qn] - 1 - x
            s_out = "-" if strand == "+" else "+"
            cnt["flipped"] += v
        else:
            s_out = strand
        if qn == "qZ_alt":
            cnt["other"] += v
            continue
        exp[(qn, s_out, x)] += v
        cnt["kept"] += v

out = os.path.join(wd, "out", "lift")
log = subprocess.run([sys.executable, script, "--plus", os.path.join(wd, "hg19_plus.bw"), "--minus",
                      os.path.join(wd, "hg19_minus.bw"), "--chain", os.path.join(wd, "syn.chain.gz"), "--fai",
                      os.path.join(wd, "syn.fa.fai"), "--keep", "qA,qB", "--out-prefix", out],
                     capture_output=True, text=True)
print(log.stdout + log.stderr)
if log.returncode:
    sys.exit("STOP: 29b_liftover.py failed on synthetic data")
got = {}
for s_, suffix in (("+", "plus"), ("-", "minus")):
    bw = pyBigWig.open(f"{out}_{suffix}.bw")
    for q in ("qA", "qB"):
        for s0, e0, v in bw.intervals(q) or []:
            got[(q, s_, s0)] = int(round(v))
    bw.close()
keys = set(exp) | set(got)
bad = sum(exp.get(k, 0) != got.get(k, 0) for k in keys)
import csv
rows = list(csv.DictReader(open(f"{out}_liftover.tsv"), delimiter="\t"))
tab = {k: sum(int(r[k]) for r in rows) for k in ("unmapped_reads", "ambiguous_reads", "flipped_reads", "other_seq_reads",
                                                  "kept_reads")}
for k_exp, k_tab in (("unmapped", "unmapped_reads"), ("ambiguous", "ambiguous_reads"), ("flipped", "flipped_reads"),
                     ("other", "other_seq_reads"), ("kept", "kept_reads")):
    if cnt[k_exp] != tab[k_tab]:
        print(f"DIFF {k_tab}: expected {cnt[k_exp]}, table {tab[k_tab]}")
        bad += 1
print(f"independent mapper: kept {cnt['kept']}, unmapped {cnt['unmapped']}, ambiguous {cnt['ambiguous']}, "
      f"flipped {cnt['flipped']}, other sequence {cnt['other']}; positions with added reads "
      f"{sum(1 for k in exp if exp[k] > 29)}")
print(f"compared {len(keys) + 5} values, {bad} differ")
