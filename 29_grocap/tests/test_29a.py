#!/usr/bin/env python3
"""Synthetic test of 29a_inspect.py: writes small bigWig, bedGraph and chain files with known content into the
folder given, runs the inspection and compares the printed totals with the known ones.
Usage: python3 test_29a.py <inspect script> <empty work folder>"""
import gzip
import os
import re
import subprocess
import sys

import numpy as np
import pyBigWig

script, wd = sys.argv[1:3]
os.makedirs(wd, exist_ok=True)
rng = np.random.default_rng(7)
chroms = [("chr1", 1_000_000), ("chr7", 6_000_000), ("chr12", 7_000_000), ("chrX", 500_000), ("chrUn_gl000220", 50_000)]
expect = {}
for strand, sign in (("plus", 1), ("minus", -1)):
    bw = pyBigWig.open(os.path.join(wd, f"syn_{strand}.bigWig"), "w")
    bw.addHeader(chroms)
    tot_a = tot_x = tot_o = 0.0
    land = {"GAPDH": 0.0, "ACTB": 0.0}
    for c, ln in chroms:
        pos = np.sort(rng.choice(ln - 10, 300, replace=False))
        if c == "chr12" and strand == "plus":
            pos = np.unique(np.concatenate([pos, [6643584, 6643600]]))
        if c == "chr7" and strand == "minus":
            pos = np.unique(np.concatenate([pos, [5570231, 5570100]]))
        val = rng.integers(1, 20, len(pos)).astype(float)
        spans = np.where(np.arange(len(pos)) % 10 == 0, 2, 1)      # some intervals two bases wide
        ends = pos + spans
        ok = np.concatenate([[True], pos[1:] >= ends[:-1]])
        pos, ends, val = pos[ok], ends[ok], val[ok]
        bw.addEntries([c] * len(pos), pos.tolist(), ends=ends.tolist(), values=(sign * val).tolist())
        t = float(np.sum(val * (ends - pos)))
        if c in ("chr1", "chr7", "chr12"):
            tot_a += t
        elif c == "chrX":
            tot_x += t
        else:
            tot_o += t
        for nm, lc, p0 in (("GAPDH", "chr12", 6643584), ("ACTB", "chr7", 5570231)):
            if c == lc:
                m = (ends > p0 - 250) & (pos <= p0 + 250)
                ov = np.minimum(ends[m], p0 + 251) - np.maximum(pos[m], p0 - 250)
                land[nm] += float(np.sum(val[m] * ov))
    bw.close()
    expect[f"syn_{strand}.bigWig"] = (tot_a, tot_x, tot_o, land["GAPDH"], land["ACTB"])
with gzip.open(os.path.join(wd, "syn.bedGraph.gz"), "wt") as fh:
    fh.write("track type=bedGraph\nchr1\t100\t101\t3\nchr1\t200\t203\t2\nchrX\t5\t6\t4\n")
expect["syn.bedGraph.gz"] = (9.0, 4.0, 0.0, 0.0, 0.0)
with gzip.open(os.path.join(wd, "syn.over.chain.gz"), "wt") as fh:
    fh.write("chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1\n100\n\n"
             "chain 900 chr2 1000 + 0 50 chr2 1000 - 10 60 2\n50\n\n")
files = [os.path.join(wd, f) for f in ("syn_plus.bigWig", "syn_minus.bigWig", "syn.bedGraph.gz", "syn.over.chain.gz")]
out = subprocess.run([sys.executable, script] + files, capture_output=True, text=True, check=True).stdout
print(out)
bad = 0
for name, (a, x, o, g, b) in expect.items():
    blk = out.split(f"== {name}")[1].split("\n== ")[0]
    m = re.search(r"autosomes ([\d,.]+); chrX ([\d,.]+); chrY [\d,.]+; chrM [\d,.]+; other sequences ([\d,.]+)", blk)
    got = [float(v.replace(",", "")) for v in m.groups()]
    lm = re.search(r"GAPDH \(\+\) ([\d,.]+); ACTB \(-\) ([\d,.]+)", blk)
    got += [float(v.replace(",", "")) for v in lm.groups()]
    for e_, g_ in zip((a, x, o, g, b), got):
        if abs(e_ - g_) > 0.051:
            print(f"DIFF {name}: expected {e_}, printed {g_}")
            bad += 1
if "2 chains, 1 with the query on the minus strand" not in out:
    print("DIFF chain summary")
    bad += 1
if "negative values: 0" in out.split("== syn_minus.bigWig")[1].split("\n== ")[0]:
    print("DIFF minus strand negatives not counted")
    bad += 1
print(f"compared {len(expect) * 5 + 2} values, {bad} differ")
