#!/usr/bin/env python3
"""Step 03: recover per-arm coordinates for the iSTART divergent pairs.

directionality.tsv stores only the midpoint of each pair, but sequence models of
initiation (Puffin-D, ProCapNet) need the position of each arm. The PINTS
bidirectional peak files hold them: column 5 is the plus-strand arm, column 6 the
minus-strand arm (see read_bi in iSTART/scripts/12_itss.py), and step 26 wrote
mid = (col5 + col6) // 2. We therefore rebuild mid from the peak file and join on
(chrom, mid).

Arm orientation follows step 26:
  * intronic pairs ('intron genu'): the sense arm is the one on the host strand,
    taken from the host gene's strand in the GTF;
  * intergenic pairs: sense = plus arm, antisense = minus arm (arbitrary, as in
    step 26, which calls sig(chrom, pp, '+') and sig(chrom, pm, '-')).

Rows whose midpoint matches no peak, or matches several peaks with different arm
positions, are written to the output with empty coordinates and counted in the
summary; they are not silently dropped.

Usage:
  python3 03_arm_coordinates.py --table input/istart_directionality_K562.tsv \
      --bi input/bi_K562_ENCFF819CPA.bed.gz --gtf <gencode.v47.annotation.gtf.gz> \
      --line K562 --out results/arms_K562.tsv
"""
import argparse
import gzip
import sys
from collections import defaultdict

import pandas as pd

INTRON_CLASS = "intron genu"
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]


def open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_bi(path):
    """PINTS bidirectional peaks -> {(chrom, mid): set of (plus_pos, minus_pos, confidence)}."""
    index = defaultdict(set)
    n = 0
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            plus_pos, minus_pos = int(f[4]), int(f[5])
            index[(f[0], (plus_pos + minus_pos) // 2)].add((plus_pos, minus_pos, f[3]))
            n += 1
    return index, n


def read_gene_strands(path):
    """gene_id -> strand, from the GTF gene records."""
    strands = {}
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = None
            for item in f[8].split(";"):
                item = item.strip()
                if item.startswith("gene_id "):
                    gid = item.split('"')[1]
                    break
            if gid:
                strands[gid] = f[6]
    return strands


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True)
    ap.add_argument("--bi", required=True, help="PINTS bidirectional peaks (BED, may be gzipped)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--line", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    d = pd.read_csv(a.table, sep="\t")
    if list(d.columns) != EXPECTED_COLUMNS:
        sys.exit(f"{a.table}: unexpected columns {list(d.columns)}")
    index, n_peaks = read_bi(a.bi)
    print(f"[{a.line}] {n_peaks} peaks in {a.bi}; {len(index)} distinct (chrom, mid) keys")
    strands = read_gene_strands(a.gtf)
    print(f"[{a.line}] {len(strands)} gene strands from {a.gtf}")

    rows, n_missing, n_ambiguous, n_nostrand = [], 0, 0, 0
    for r in d.itertuples(index=False):
        hits = index.get((r.chrom, int(r.mid)), set())
        positions = {(h[0], h[1]) for h in hits}
        plus_pos = minus_pos = conf = None
        if not positions:
            n_missing += 1
        elif len(positions) > 1:
            n_ambiguous += 1
        else:
            plus_pos, minus_pos = positions.pop()
            conf = sorted({h[2] for h in hits})[0]
        host_strand = ""
        if r.klasa == INTRON_CLASS:
            host_strand = strands.get(r.host, "")
            if not host_strand:
                n_nostrand += 1
        sense_strand = host_strand if host_strand else "+"
        if plus_pos is None:
            pos_sense = pos_anti = ""
        elif sense_strand == "+":
            pos_sense, pos_anti = plus_pos, minus_pos
        else:
            pos_sense, pos_anti = minus_pos, plus_pos
        rows.append({"line": a.line, "chrom": r.chrom, "mid": r.mid, "klasa": r.klasa,
                     "host": r.host, "host_name": r.host_name, "host_strand": host_strand,
                     "host_signal": r.host_signal, "pos_sense": pos_sense, "pos_antisense": pos_anti,
                     "strand_sense": sense_strand if plus_pos is not None else "",
                     "strand_antisense": ("-" if sense_strand == "+" else "+") if plus_pos is not None else "",
                     "pos_plus": "" if plus_pos is None else plus_pos,
                     "pos_minus": "" if minus_pos is None else minus_pos,
                     "arm_distance": "" if plus_pos is None else abs(plus_pos - minus_pos),
                     "confidence": conf or "", "sens": r.sens, "anty": r.anty, "DI": r.DI,
                     "host_tx": r.host_tx, "k36": r.k36,
                     "status": "ok" if plus_pos is not None else ("ambiguous" if len(hits) > 1 else "no_peak")})
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    ok = out[out["status"] == "ok"]
    print(f"[{a.line}] {len(out)} pairs: {len(ok)} resolved, {n_missing} without a matching peak, "
          f"{n_ambiguous} ambiguous, {n_nostrand} intronic rows without a host strand")
    if len(ok):
        print(f"[{a.line}] arm distance (nt): median {ok['arm_distance'].median():.0f}, "
              f"min {ok['arm_distance'].min():.0f}, max {ok['arm_distance'].max():.0f}")
        print(out.groupby(["klasa", "status"]).size().to_string())
        print(f"[{a.line}] sense arm on plus strand: "
              f"{100 * (ok['strand_sense'] == '+').mean():.1f}% of resolved pairs")
    print(f"written: {a.out}")


if __name__ == "__main__":
    main()
