#!/usr/bin/env python3
"""Step 29a: format check of the GRO-cap files (GEO GSM1480321 with TAP, GSM1480322 without TAP; hg19) and of
the UCSC liftover chain. No strand counts inside genes are computed here: the step only answers what the files
contain, so that the liftover (step 29b) can be written for the real format.

For every signal file (bigWig, or bedGraph, plain or gzipped) it prints: the chromosome names, the number of
intervals, how many span more than one base, the smallest non-zero absolute value and whether all values are whole
multiples of it (raw counts), how many values are negative, the total signal (value x span) on the 22 autosomes,
chrX, chrY, chrM and other sequences, and the signal on each strand within 250 nt of two well-known starts in hg19
(GAPDH, + strand, chr12:6,643,585; ACTB, - strand, chr7:5,570,232), which shows which file holds which strand.
For the chain file it prints the number of chains, the number of chains with an inverted query strand and the
assembly names in the first chain.

Usage: python3 29a_inspect.py <file> [<file> ...]
"""
import gzip
import os
import sys

import numpy as np

AUTOS = [f"chr{i}" for i in range(1, 23)]
LANDMARKS = [("GAPDH", "chr12", 6643584, "+"), ("ACTB", "chr7", 5570231, "-")]   # 0-based first base of the gene


def summarise(name, chrom_iter):
    """chrom_iter yields (chrom, starts, ends, values) as numpy arrays."""
    tot = {}
    n_iv = n_wide = n_neg = 0
    vals_all = []
    land = {lm[0]: 0.0 for lm in LANDMARKS}
    chroms = []
    for c, s, e, v in chrom_iter:
        chroms.append(c)
        if len(v) == 0:
            continue
        span = (e - s).astype(np.int64)
        n_iv += len(v)
        n_wide += int((span > 1).sum())
        n_neg += int((v < 0).sum())
        tot[c] = tot.get(c, 0.0) + float(np.sum(np.abs(v) * span))
        nz = np.abs(v[v != 0])
        if len(nz):
            vals_all.append(nz[: 200000])
        for nm, lc, pos, _ in LANDMARKS:
            if c == lc:
                m = (e > pos - 250) & (s <= pos + 250)
                ov = np.minimum(e[m], pos + 251) - np.maximum(s[m], pos - 250)
                land[nm] += float(np.sum(np.abs(v[m]) * ov))
    print(f"\n== {name}")
    print(f"  sequences: {len(chroms)} ({', '.join(chroms[:6])}{', ...' if len(chroms) > 6 else ''})")
    print(f"  intervals: {n_iv:,}; spanning more than one base: {n_wide:,}; negative values: {n_neg:,}")
    if vals_all:
        nz = np.concatenate(vals_all)
        unit = float(nz.min())
        ratio = nz / unit
        whole = float(np.mean(np.abs(ratio - np.round(ratio)) < 1e-6))
        print(f"  smallest non-zero |value| {unit:.6g}; share of values that are whole multiples of it {whole:.4f}; "
              f"largest |value| {float(nz.max()):.6g}; median {float(np.median(nz)):.6g}")
    a = sum(tot.get(c, 0.0) for c in AUTOS)
    other = sum(v for k, v in tot.items() if k not in AUTOS + ["chrX", "chrY", "chrM"])
    print(f"  total |signal| x span: autosomes {a:,.1f}; chrX {tot.get('chrX', 0):,.1f}; chrY {tot.get('chrY', 0):,.1f}; "
          f"chrM {tot.get('chrM', 0):,.1f}; other sequences {other:,.1f}")
    print("  within 250 nt of hg19 starts: " + "; ".join(f"{nm} ({st}) {land[nm]:,.1f}" for nm, _, _, st in LANDMARKS))


def bigwig_iter(path):
    import pyBigWig
    bw = pyBigWig.open(path)
    if bw is None:
        sys.exit(f"STOP: {path} cannot be opened as bigWig")
    print(f"\n-- {os.path.basename(path)}: bigWig header {bw.header()}")
    for c, ln in bw.chroms().items():
        iv = bw.intervals(c) or []
        if not iv:
            yield c, np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
            continue
        a = np.array(iv, dtype=float)
        yield c, a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2]
    bw.close()


def bedgraph_iter(path):
    op = gzip.open if path.endswith(".gz") else open
    data = {}
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(("track", "browser", "#")) or not ln.strip():
                continue
            f = ln.split()
            data.setdefault(f[0], []).append((int(f[1]), int(f[2]), float(f[3])))
    print(f"\n-- {os.path.basename(path)}: bedGraph, {sum(len(v) for v in data.values()):,} lines")
    for c, rows in data.items():
        a = np.array(rows, dtype=float)
        yield c, a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2]


def chain_info(path):
    n = n_inv = 0
    first = None
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("chain"):
                f = ln.split()
                n += 1
                n_inv += f[9] == "-"
                if first is None:
                    first = f
    print(f"\n== {os.path.basename(path)}: {n:,} chains, {n_inv:,} with the query on the minus strand; first chain: "
          f"target {first[2]} ({first[4]}), query {first[7]} ({first[9]})")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        low = p.lower()
        if low.endswith(".chain.gz"):
            chain_info(p)
        elif low.endswith((".bw", ".bigwig")):
            summarise(os.path.basename(p), bigwig_iter(p))
        elif low.endswith((".bedgraph", ".bedgraph.gz", ".bg", ".bg.gz", ".bed.gz", ".bed")):
            summarise(os.path.basename(p), bedgraph_iter(p))
        else:
            print(f"\n== {os.path.basename(p)}: format not recognised from the name; first bytes: "
                  f"{open(p, 'rb').read(64)!r}")


if __name__ == "__main__":
    main()
