#!/usr/bin/env python3
"""Step 29b: lift single-base GRO-cap counts (bigWig, one value per base, as checked in step 29a) from hg19 to
hg38 with a UCSC chain file, keeping the strand of every read.

Rules (recorded in the step-29 pre-specification)
  - every hg19 position is mapped through the ungapped blocks of the chain (target = hg19, query = hg38);
  - a position inside a block of a chain whose query strand is '-' moves to the other strand in hg38;
  - positions outside every block are dropped (counted as unmapped);
  - positions covered by blocks of more than one chain are dropped (counted as ambiguous);
  - reads of several hg19 positions that land on the same hg38 position and strand are added;
  - only hg38 sequences listed in --keep are written (default chr1-22, chrX, chrY, chrM); reads that land on other
    sequences are counted and dropped. Their lengths are taken from --fai and must equal the chain's query sizes.
Values are taken as absolute counts (the minus-strand file holds negative values); the output holds positive counts
on both strands.

Usage: python3 29b_liftover.py --plus hg19_plus.bw --minus hg19_minus.bw --chain hg19ToHg38.over.chain.gz
                               --fai hg38.fa.fai --out-prefix out/29b_GROcap_hg38
Writes <prefix>_plus.bw, <prefix>_minus.bw and <prefix>_liftover.tsv (counts per hg19 chromosome and outcome);
prints a log.
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pyBigWig


def read_chain(path):
    """Ungapped blocks per target chromosome: arrays tstart, size, qname index, qstart (query-strand coordinates),
    qstrand (+1/-1), qsize, chain id."""
    blocks = defaultdict(list)
    qnames, qsizes = [], {}
    cid = -1
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for ln in fh:
            f = ln.split()
            if not f:
                continue
            if f[0] == "chain":
                cid += 1
                tname, tpos = f[2], int(f[5])
                qname, qsize, qstrand, qpos = f[7], int(f[8]), f[9], int(f[10])
                if f[4] != "+":
                    sys.exit("STOP: chain with the target on the minus strand")
                if qname not in qsizes:
                    qnames.append(qname)
                    qsizes[qname] = qsize
                elif qsizes[qname] != qsize:
                    sys.exit(f"STOP: two sizes for {qname} in the chain")
                qi = qnames.index(qname)
                qs = 1 if qstrand == "+" else -1
                continue
            size = int(f[0])
            blocks[tname].append((tpos, size, qi, qpos, qs, cid))
            if len(f) == 3:
                tpos += size + int(f[1])
                qpos += size + int(f[2])
    out = {}
    for t, b in blocks.items():
        a = np.array(sorted(b), dtype=np.int64)
        out[t] = {"tstart": a[:, 0], "size": a[:, 1], "qi": a[:, 2], "qstart": a[:, 3], "qs": a[:, 4], "cid": a[:, 5]}
    return out, qnames, qsizes


def segments(B):
    """Cut the target into elementary segments between block boundaries; for each segment the number of blocks
    covering it and, where exactly one block covers it, that block's index."""
    s, e = B["tstart"], B["tstart"] + B["size"]
    bnd = np.unique(np.concatenate([s, e]))
    i0, i1 = np.searchsorted(bnd, s), np.searchsorted(bnd, e)
    cov = np.zeros(len(bnd) + 1, np.int64)
    idx = np.zeros(len(bnd) + 1, np.int64)
    k = np.arange(len(s), dtype=np.int64)
    np.add.at(cov, i0, 1)
    np.add.at(cov, i1, -1)
    np.add.at(idx, i0, k)
    np.add.at(idx, i1, -k)
    return bnd, np.cumsum(cov)[: len(bnd)], np.cumsum(idx)[: len(bnd)]


def lift(pos, B, seg):
    """Map 0-based hg19 positions; returns (ok, ambiguous, qi, qpos on the query strand of the chain, flip)."""
    bnd, cov, idx = seg
    j = np.searchsorted(bnd, pos, side="right") - 1
    inside = (j >= 0) & (j < len(bnd) - 1)
    jj = np.clip(j, 0, len(bnd) - 1)
    c = np.where(inside, cov[jj], 0)
    ok, amb = c == 1, c >= 2
    kk = np.where(ok, idx[jj], 0)
    off = pos - B["tstart"][kk]
    if np.any(ok & ((off < 0) | (off >= B["size"][kk]))):
        sys.exit("STOP: segment lookup outside its block")
    return ok, amb, B["qi"][kk], B["qstart"][kk] + off, B["qs"][kk] < 0


def read_bw(path):
    bw = pyBigWig.open(path)
    if bw is None:
        sys.exit(f"STOP: cannot open {path}")
    out = {}
    for c in bw.chroms():
        iv = bw.intervals(c) or []
        if not iv:
            continue
        a = np.array(iv, dtype=float)
        if np.any(a[:, 1] - a[:, 0] != 1):
            sys.exit(f"STOP: {path} {c}: intervals wider than one base (step 29a found none)")
        v = np.abs(a[:, 2])
        if np.any(np.abs(v - np.round(v)) > 1e-6):
            sys.exit(f"STOP: {path} {c}: values are not whole counts")
        out[c] = (a[:, 0].astype(np.int64), np.round(v).astype(np.int64))
    bw.close()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plus", required=True)
    ap.add_argument("--minus", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--keep", default=",".join([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]))
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    keep = [c for c in a.keep.split(",") if c]
    fai = {}
    with open(a.fai) as fh:
        for ln in fh:
            f = ln.split("\t")
            fai[f[0]] = int(f[1])
    chain, qnames, qsizes = read_chain(a.chain)
    for c in keep:
        if c not in fai:
            sys.exit(f"STOP: {c} not in {a.fai}")
        if c in qsizes and qsizes[c] != fai[c]:
            sys.exit(f"STOP: {c}: length {fai[c]} in the fai but {qsizes[c]} in the chain")
    print(f"== 29b liftover: chain {a.chain}: {sum(len(v['size']) for v in chain.values()):,} blocks on "
          f"{len(chain)} hg19 sequences; {len(qnames)} hg38 sequences")
    segs = {t: segments(B) for t, B in chain.items()}
    n_amb = sum(int(np.sum(np.diff(bn)[cv[:-1] >= 2])) for bn, cv, _ in segs.values())
    print(f"  hg19 bases covered by more than one chain block: {n_amb:,}")
    acc = {(c, s): defaultdict(int) for c in keep for s in ("+", "-")}
    rows = []
    for strand_in, path in (("+", a.plus), ("-", a.minus)):
        D = read_bw(path)
        for c, (pos, val) in D.items():
            n_all, r_all = len(pos), int(val.sum())
            if c not in chain:
                rows.append({"file_strand": strand_in, "hg19_chrom": c, "positions": n_all, "reads": r_all,
                             "unmapped_reads": r_all, "ambiguous_reads": 0, "flipped_reads": 0, "other_seq_reads": 0,
                             "kept_reads": 0})
                continue
            ok, amb, qi, qp, flip = lift(pos, chain[c], segs[c])
            qname = np.array(qnames, dtype=object)[qi]
            qsz = np.array([qsizes[q] for q in qnames], dtype=np.int64)[qi]
            hpos = np.where(flip, qsz - 1 - qp, qp)
            st_out = np.where(flip, {"+": "-", "-": "+"}[strand_in], strand_in)
            kept_mask = ok & np.isin(qname, keep)
            for q, p_, s_, v_ in zip(qname[kept_mask], hpos[kept_mask], st_out[kept_mask], val[kept_mask]):
                acc[(q, s_)][int(p_)] += int(v_)
            rows.append({"file_strand": strand_in, "hg19_chrom": c, "positions": n_all, "reads": r_all,
                         "unmapped_reads": int(val[~ok & ~amb].sum()), "ambiguous_reads": int(val[amb].sum()),
                         "flipped_reads": int(val[ok & flip].sum()),
                         "other_seq_reads": int(val[ok & ~np.isin(qname, keep)].sum()),
                         "kept_reads": int(val[kept_mask].sum())})
    tot = {k: sum(r[k] for r in rows) for k in ("reads", "unmapped_reads", "ambiguous_reads", "flipped_reads",
                                              "other_seq_reads", "kept_reads")}
    print(f"  reads in: {tot['reads']:,}; kept: {tot['kept_reads']:,} ({100 * tot['kept_reads'] / tot['reads']:.2f}%); "
          f"unmapped {tot['unmapped_reads']:,}; ambiguous {tot['ambiguous_reads']:,}; on other hg38 sequences "
          f"{tot['other_seq_reads']:,}; moved to the other strand (inverted blocks) {tot['flipped_reads']:,}")
    if tot["kept_reads"] + tot["unmapped_reads"] + tot["ambiguous_reads"] + tot["other_seq_reads"] != tot["reads"]:
        sys.exit("STOP: reads do not add up")
    for s_, suffix in (("+", "plus"), ("-", "minus")):
        bw = pyBigWig.open(f"{a.out_prefix}_{suffix}.bw", "w")
        bw.addHeader([(c, fai[c]) for c in keep])
        for c in keep:
            d = acc[(c, s_)]
            if not d:
                continue
            p = np.array(sorted(d), dtype=np.int64)
            v = np.array([d[x] for x in p], dtype=float)
            if p[-1] >= fai[c]:
                sys.exit(f"STOP: {c}: lifted position beyond the end of the sequence")
            bw.addEntries([c] * len(p), p.tolist(), ends=(p + 1).tolist(), values=v.tolist())
        bw.close()
    out_reads = {}
    for s_, suffix in (("+", "plus"), ("-", "minus")):
        bw = pyBigWig.open(f"{a.out_prefix}_{suffix}.bw")
        out_reads[s_] = sum(int(round(sum(x[2] for x in (bw.intervals(c) or [])))) for c in keep)
        bw.close()
    if out_reads["+"] + out_reads["-"] != tot["kept_reads"]:
        sys.exit("STOP: the written bigWigs do not hold the kept reads")
    auto = [f"chr{i}" for i in range(1, 23)]
    n_auto = sum(sum(acc[(c, s_)].values()) for c in auto for s_ in ("+", "-") if c in keep)
    print(f"  written: {a.out_prefix}_plus.bw ({out_reads['+']:,} reads), {a.out_prefix}_minus.bw "
          f"({out_reads['-']:,} reads); autosomal reads in hg38: {n_auto:,}")
    import csv
    with open(f"{a.out_prefix}_liftover.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  table: {a.out_prefix}_liftover.tsv")


if __name__ == "__main__":
    main()
