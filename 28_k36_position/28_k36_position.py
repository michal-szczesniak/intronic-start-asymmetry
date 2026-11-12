#!/usr/bin/env python3
"""Step 28: H3K36me3, position along the gene and a second cell line for the gene-strand deficit of strong
measured start sites.

Background. Step 17 found fewer strong PRO-cap start sites (>= 10 reads) on the gene strand of introns in K562,
clearly so only in the more active genes (17% in active protein-coding genes); step 26b (DNA methylation) gave no
clear answer. Inside transcribed genes, SETD2 lays down H3K36me3, which rises toward the 3' end of genes and
suppresses spurious initiation (Carvalho et al. 2013; Neri et al. 2017).

Questions, specified on 10 Nov 2026 before any data were seen
  A1  Is the deficit larger inside H3K36me3 peaks than outside them, in the same active genes?
      Contrast: deficit(inside) - deficit(outside); prediction > 0.
  A2  Does the deficit grow toward the 3' end of the gene? Each gene (annotated span, 5' to 3') is cut into fifths
      P1 (5') ... P5 (3'). Contrast: deficit(P5) - deficit(P1); prediction > 0.
  B   Is the deficit of strong sites in active genes reproduced in HCT116? deficit(all positions); prediction > 0.
  Primary analyses: active protein-coding genes (host-activity thirds A1-A3 pooled, genes resampled within each
  third), strong sites; K562: A1 and A2; HCT116: A1, A2 and B. Five contrasts, two-sided, no correction for
  multiple testing. Everything else is secondary.
  Strong sites: >= max(3, round(10 * N / N_ref)) reads, N = reads of the data set on the 22 autosomes, N_ref =
  18,699,448 (the K562 data of step 17), i.e. 10 reads in K562 and the same reads per million elsewhere. Sites
  with >= 10 and >= 3 reads regardless of depth: secondary.
  Controls (no difference expected): silent protein-coding genes (inside - outside, P5 - P1); intergenic 1-Mb
  blocks of step 10, plus strand standing for the gene strand, inside and outside peaks.
  Caveat: the host gene's own transcription adds reads on the gene strand, more so in active genes and inside
  H3K36me3 peaks (step 17: an excess, not a deficit, of sites with >= 1 read). This background makes the deficit
  of strong sites smaller where transcription is high, i.e. it works against A1.

Definitions (as in step 17, whose code is loaded and whose per-gene counts are reproduced, see below)
  eligible positions  step 10's masks: protein-coding and lncRNA gene bodies minus exons +/-100 nt, annotated TSSs
                      +/-1 kb and positions of other genes; intergenic positions > 10 kb from any gene; not N
  signal              plus- and minus-strand bigWigs of RNA 5' ends; several files per strand (replicates) are
                      summed; absolute values; the read unit is the smallest non-zero value (checked: values must
                      be whole multiples of it); strands swapped if most TSS reads fall on the antisense strand
  site                eligible positions with >= T reads on one strand, merged within 50 nt (step 10's rule);
                      its position is the position with most reads (the first if tied)
  strata of a site    inside / outside a replicated H3K36me3 peak (ENCODE bed, unstranded), at the site's position;
                      fifth of the gene from the relative position along the annotated gene, 5' to 3'
  host activity       sense reads within +/- 250 nt of the gene's annotated TSSs: silent 0, low 1-4, active >= 5
                      (thirds A1-A3 among active protein-coding genes)
  deficit             100 * (1 - sense sites / antisense sites), summed over genes; CIs: 2,000 bootstrap
                      resamples of genes (blocks); contrasts use the same resamples for both strata (paired)
Validation: with --genes17 (K562), the per-gene counts of sites with >= 10 and >= 3 reads over all eligible
positions, the eligible nt and the host activity must equal step 17's table (17_obs_genes.tsv); STOP otherwise.

Outputs: <prefix>_genes.tsv (per gene: activity, eligible nt per stratum, sites per threshold, stratum and strand),
<prefix>_blocks.tsv, <prefix>_strata.tsv (every deficit and contrast; column 'primary' marks the five primary
contrasts); the log is printed.

Usage (CPU, about 15 min and 16 GB per data set):
  python3 28_k36_position.py --line K562 --plus-bw ENCFF798GNW.bigWig --minus-bw ENCFF662SHP.bigWig \
      --peaks ENCFF520LHY.bed.gz --primary A1,A2 --genes17 ../17_observed/results/17_obs_genes.tsv \
      --chroms chr1,...,chr22 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --step10 ../10_shuffle/10_shuffle_test.py --step17 ../17_observed/17_observed_sites.py --out-prefix results/28_K562
"""
import argparse
import gzip
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

N_REF = 18_699_448          # reads of the K562 data of step 17 on the 22 autosomes (17_obs.log)
STORE_MIN = 3               # reads; positions with fewer reads are not needed for any threshold
ACT = ["silent", "low", "A1", "A2", "A3"]
BINS = ["P1", "P2", "P3", "P4", "P5"]
G_STRATA = ["all", "in", "out"] + BINS
B_STRATA = ["all", "in", "out"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read_peaks(path, chroms):
    """Peak intervals (0-based, half-open) per chromosome from a BED file (first three fields)."""
    out = {c: [] for c in chroms}
    n_all = n_bad = 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.split("\t")
            n_all += 1
            s, e = int(f[1]), int(f[2])
            if e <= s:
                n_bad += 1
                continue
            if f[0] in out:
                out[f[0]].append((s, e))
    return out, n_all, n_bad


def interval_mask(iv, L):
    d = np.zeros(L + 1, dtype=np.int32)
    for s, e in iv:
        d[max(0, s)] += 1
        d[min(L, e)] -= 1
    return np.cumsum(d[:-1]) > 0


def store(acc, unit_i, sense, t, sel, base, inpk, binv, vmin):
    """Keep eligible positions with >= vmin reads: unit index, strand (1 = sense / plus), genomic position,
    reads, inside a peak, fifth of the gene (-1 for intergenic blocks)."""
    idx = np.flatnonzero(sel & (t >= vmin))
    if not len(idx):
        return
    acc["unit"].append(np.full(len(idx), unit_i, dtype=np.int64))
    acc["sense"].append(np.full(len(idx), sense, dtype=np.int8))
    acc["pos"].append(base + idx.astype(np.int64))
    acc["val"].append(t[idx].astype(np.float64))
    acc["pk"].append(inpk[idx].astype(bool))
    acc["bin"].append(binv[idx].astype(np.int8) if binv is not None else np.full(len(idx), -1, dtype=np.int8))


def concat(acc):
    keys = ("unit", "sense", "pos", "val", "pk", "bin")
    dt = (np.int64, np.int8, np.int64, np.float64, bool, np.int8)
    return {k: (np.concatenate(acc[k]) if acc[k] else np.zeros(0, dtype=d)) for k, d in zip(keys, dt)}


def site_reps(X, T, merge):
    """Indices (into X) of the representative position of every site with >= T reads: positions of one unit and
    strand with >= T reads, merged when <= merge nt apart; representative = most reads, then leftmost."""
    k = np.flatnonzero(X["val"] >= T * (1 - 1e-6))
    if not len(k):
        return k
    o = k[np.lexsort((X["pos"][k], X["sense"][k], X["unit"][k]))]
    u, s, p, v = X["unit"][o], X["sense"][o], X["pos"][o], X["val"][o]
    new = np.ones(len(o), dtype=bool)
    new[1:] = (u[1:] != u[:-1]) | (s[1:] != s[:-1]) | ((p[1:] - p[:-1]) > merge)
    cid = np.cumsum(new) - 1
    o2 = np.lexsort((p, -v, cid))
    first = np.ones(len(o2), dtype=bool)
    first[1:] = cid[o2][1:] != cid[o2][:-1]
    return o[o2[first]]


def count_sites(X, reps, n_units, strata, tag, table):
    """Add columns <tag>_<stratum>_<sense|anti> (sites per unit) to the dict-of-arrays table."""
    u, s, pk, b = X["unit"][reps], X["sense"][reps], X["pk"][reps], X["bin"][reps]
    for st in strata:
        if st == "all":
            m = np.ones(len(reps), dtype=bool)
        elif st == "in":
            m = pk
        elif st == "out":
            m = ~pk
        else:
            m = b == BINS.index(st)
        for flag, nm in ((1, "sense"), (0, "anti")):
            table[f"{tag}_{st}_{nm}"] = np.bincount(u[m & (s == flag)], minlength=n_units).astype(np.int64)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", required=True, help="label of the data set, e.g. K562 or HCT116")
    ap.add_argument("--plus-bw", required=True, help="comma-separated plus-strand bigWigs, summed")
    ap.add_argument("--minus-bw", required=True, help="comma-separated minus-strand bigWigs, summed")
    ap.add_argument("--peaks", required=True, help="H3K36me3 peaks (BED, may be gzipped)")
    ap.add_argument("--primary", default="", help="primary contrasts of this data set: any of A1,A2,B")
    ap.add_argument("--genes17", default="", help="step-17 gene table; if given, its counts must be reproduced")
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--step17", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--ref-reads", type=float, default=N_REF)
    ap.add_argument("--tss-window", type=int, default=250)
    ap.add_argument("--swap-strands", choices=("auto", "yes", "no"), default="auto")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    chroms = [c for c in a.chroms.split(",") if c]
    primary = {x for x in a.primary.split(",") if x}
    if primary - {"A1", "A2", "B"}:
        sys.exit(f"STOP: unknown primary contrast(s) {sorted(primary - {'A1', 'A2', 'B'})}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    s10 = load_module("s10", os.path.abspath(a.step10))
    s17 = load_module("s17", os.path.abspath(a.step17))
    rng = np.random.default_rng(a.seed)
    print(f"== 28 {a.line}: H3K36me3, position along the gene; bootstrap {a.n_boot}, seed {a.seed}; "
          f"primary contrasts: {', '.join(sorted(primary)) or 'none'}")

    import pyBigWig
    from pyfaidx import Fasta
    bwP = [pyBigWig.open(p) for p in a.plus_bw.split(",") if p]
    bwM = [pyBigWig.open(p) for p in a.minus_bw.split(",") if p]
    print(f"  plus strand: {a.plus_bw}\n  minus strand: {a.minus_bw}\n  peaks: {a.peaks}")
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), set(chroms))
    gtss = s17.gene_tss(os.path.abspath(a.gtf), set(chroms))
    peaks, n_pk, n_bad = read_peaks(a.peaks, chroms)
    n_pk_c = sum(len(v) for v in peaks.values())
    bp_pk = sum(e - s for v in peaks.values() for s, e in v)
    print(f"  annotation: {len(genes):,} genes; peaks: {n_pk:,} in the file ({n_bad:,} empty), {n_pk_c:,} on the "
          f"analysed chromosomes covering {bp_pk / 1e6:.1f} Mb  [{time.time() - t0:.0f} s]", flush=True)
    if n_pk_c == 0:
        sys.exit("STOP: no peaks on the analysed chromosomes (chromosome names?)")

    unit, swap = None, {"yes": True, "no": False}.get(a.swap_strands)
    grows, brows, crows = [], [], []
    GA = {k: [] for k in ("unit", "sense", "pos", "val", "pk", "bin")}
    BA = {k: [] for k in ("unit", "sense", "pos", "val", "pk", "bin")}
    print("\n== 0. signal checks ==")
    for c in chroms:
        t1 = time.time()
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        up = np.where(raw >= 97, raw - 32, raw).astype(np.uint8)
        del raw
        L = len(up)
        for bw, nm in [(b, p) for b, p in zip(bwP + bwM, a.plus_bw.split(",") + a.minus_bw.split(","))]:
            bl = bw.chroms(c)
            if bl is not None and bl != L:
                sys.exit(f"STOP: {c} is {bl:,} nt in {nm} but {L:,} nt in the genome (different assembly?)")
        pk = interval_mask(peaks[c], L)
        P = np.zeros(L, dtype=np.float64)
        M = np.zeros(L, dtype=np.float64)
        negs = 0
        for bw in bwP:
            x, n_ = s17.dense(bw, c, L)
            P += np.abs(x)
            negs += n_
            del x
        for bw in bwM:
            x, n_ = s17.dense(bw, c, L)
            M += np.abs(x)
            negs += n_
            del x
        nz = np.concatenate([P[P > 0], M[M > 0]])
        if unit is None:
            if not len(nz):
                sys.exit(f"STOP: no signal on {c}")
            unit = float(nz.min())
            mult = nz / unit
            whole = float((np.abs(mult - np.round(mult)) < 1e-3).mean())
            print(f"  {c}: non-zero positions plus {int((P > 0).sum()):,}, minus {int((M > 0).sum()):,}; negative intervals "
                  f"{negs:,} (absolute values used); smallest non-zero value {unit:.6g}; whole multiples of it {100 * whole:.2f}%")
            print("  non-zero values in reads, quantiles 50/90/99/99.9%: "
                  + ", ".join(f"{np.quantile(mult, q):.4g}" for q in (0.5, 0.9, 0.99, 0.999)) + f"; max {mult.max():.4g}")
            if whole < 0.99:
                sys.exit("STOP: the values are not whole multiples of one read unit; thresholds in reads would be meaningless")
        elif len(nz) and nz.min() < unit * (1 - 1e-6):
            print(f"  WARNING: {c}: smallest non-zero value {nz.min():.6g} is below the read unit {unit:.6g}")
        del nz

        g = genes[(genes["chrom"] == c) & genes["gene_type"].isin(("protein_coding", "lncRNA"))]
        hs = {r.gene_id: (s17.window_sum(P, M, gtss.get((c, r.gene_id), []), a.tss_window, L), r.strand)
              for r in g.itertuples(index=False)}
        hs = {k: ((sp, sm) if st == "+" else (sm, sp)) for k, ((sp, sm), st) in hs.items()}
        pcids = set(g.loc[g["gene_type"] == "protein_coding", "gene_id"])
        tS = sum(v[0] for k, v in hs.items() if k in pcids)
        tA = sum(v[1] for k, v in hs.items() if k in pcids)
        share = tS / (tS + tA) if tS + tA > 0 else np.nan
        if swap is None:
            swap = bool(share < 0.5)
            print(f"  strand check on {c}: reads within +/-{a.tss_window} nt of annotated protein-coding TSSs fall "
                  f"{100 * share:.1f}% on the gene's strand -> " + ("SWAPPING plus and minus" if swap else "orientation as labelled"))
        if swap:
            P, M = M, P
            hs = {k: (v[1], v[0]) for k, v in hs.items()}
            share = 1 - share

        ok = up != ord("N")
        del up
        glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
        vmin = STORE_MIN * unit * (1 - 1e-6)
        for r, s, e, sel in glist:
            ts, ta = (P[s:e], M[s:e]) if r.strand == "+" else (M[s:e], P[s:e])
            n = e - s
            idx = np.arange(n, dtype=np.float64)
            rel = idx / n if r.strand == "+" else (n - 1 - idx) / n
            binv = np.minimum((rel * 5).astype(np.int8), 4)
            inpk = pk[s:e]
            hS, hA = hs.get(r.gene_id, (np.nan, np.nan))
            row = {"chrom": c, "gene_id": r.gene_id, "gene_name": r.gene_name, "gene_type": r.gene_type, "strand": r.strand,
                   "start": s, "end": e, "eligible_nt": int(sel.sum()), "host_tss_sense": hS / unit, "host_tss_anti": hA / unit,
                   "elig_in": int((sel & inpk).sum()), "elig_out": int((sel & ~inpk).sum())}
            for i_, b_ in enumerate(BINS):
                row[f"elig_{b_}"] = int((sel & (binv == i_)).sum())
            i = len(grows)
            store(GA, i, 1, ts / unit, sel, s, inpk, binv, vmin / unit)
            store(GA, i, 0, ta / unit, sel, s, inpk, binv, vmin / unit)
            grows.append(row)
        for b0 in range(0, L, s10.BLOCK):
            b1 = min(L, b0 + s10.BLOCK)
            sel = ig[b0:b1]
            if sel.any():
                j = len(brows)
                store(BA, j, 1, P[b0:b1] / unit, sel, b0, pk[b0:b1], None, vmin / unit)
                store(BA, j, 0, M[b0:b1] / unit, sel, b0, pk[b0:b1], None, vmin / unit)
                brows.append({"chrom": c, "block_start": b0, "eligible_nt": int(sel.sum()),
                              "elig_in": int((sel & pk[b0:b1]).sum()), "elig_out": int((sel & ~pk[b0:b1]).sum())})
        crows.append({"chrom": c, "reads_plus": float(P.sum()) / unit, "reads_minus": float(M.sum()) / unit,
                      "tss_sense_share": share})
        print(f"  {c}: {len(glist):,} genes; reads plus {crows[-1]['reads_plus']:,.0f}, minus {crows[-1]['reads_minus']:,.0f}; "
              f"TSS reads on the gene's strand {100 * share:.1f}%; peaks {len(peaks[c]):,}  [{time.time() - t1:.0f} s]", flush=True)
        del P, M, ok, pk, ig, glist

    CR = pd.DataFrame(crows)
    N = float(CR["reads_plus"].sum() + CR["reads_minus"].sum())
    T_strong = max(3, int(round(10 * N / a.ref_reads)))
    print(f"  genome: reads plus {CR['reads_plus'].sum():,.0f}, minus {CR['reads_minus'].sum():,.0f}, total {N:,.0f} "
          f"({N / a.ref_reads:.3f} x the K562 reference); strands {'SWAPPED' if swap else 'as labelled'}")
    print(f"  strong sites: >= {T_strong} reads (max(3, round(10 x {N / a.ref_reads:.3f}))); secondary: >= 10 and >= 3 reads")

    G = pd.DataFrame(grows)
    B = pd.DataFrame(brows)
    XG, XB = concat(GA), concat(BA)
    METRICS = [("strong", T_strong), ("r10", 10), ("r3", 3)]
    gcols, bcols = {}, {}
    for tag, T in METRICS:
        count_sites(XG, site_reps(XG, T, s10.MERGE), len(G), G_STRATA, tag, gcols)
        count_sites(XB, site_reps(XB, T, s10.MERGE), len(B), B_STRATA, tag, bcols)
    G = pd.concat([G, pd.DataFrame(gcols)], axis=1)
    B = pd.concat([B, pd.DataFrame(bcols)], axis=1)
    for tag, _ in METRICS:                      # strata must add up to all sites
        for nm in ("sense", "anti"):
            assert (G[f"{tag}_in_{nm}"] + G[f"{tag}_out_{nm}"] == G[f"{tag}_all_{nm}"]).all()
            assert (sum(G[f"{tag}_{b}_{nm}"] for b in BINS) == G[f"{tag}_all_{nm}"]).all()
            assert (B[f"{tag}_in_{nm}"] + B[f"{tag}_out_{nm}"] == B[f"{tag}_all_{nm}"]).all()

    # ---- validation against step 17
    if a.genes17:
        S17 = pd.read_csv(a.genes17, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
        J = G.merge(S17, on=["chrom", "gene_id"], how="outer", suffixes=("", "_17"), indicator=True)
        only = J["_merge"].value_counts()
        both = J[J["_merge"] == "both"]
        checks = [("eligible_nt", "eligible_nt_17")] + [(f"{t}_all_{s}", f"obs_all_{t}_{s}") for t in ("r10", "r3")
                                                          for s in ("sense", "anti")]
        print(f"\n== validation against step 17 ({os.path.basename(a.genes17)}) ==\n  genes in both tables {int(only.get('both', 0)):,}; "
              f"only here {int(only.get('left_only', 0)):,}; only in step 17 {int(only.get('right_only', 0)):,}")
        bad = []
        for mine, theirs in checks:
            d = int((both[mine].to_numpy(float) != both[theirs].to_numpy(float)).sum())
            print(f"  {mine:16s} vs {theirs:20s}: {d:,} genes differ")
            if d:
                bad.append(mine)
        hdiff = np.abs(both["host_tss_sense"].to_numpy(float) - both["host_tss_sense_17"].to_numpy(float))
        nh = int((hdiff > 1e-6).sum())
        print(f"  host_tss_sense: {nh:,} genes differ")
        if bad or nh or only.get("left_only", 0) or only.get("right_only", 0):
            sys.exit("STOP: step 17 is not reproduced (code, inputs or annotation differ)")
        print("  step 17 reproduced exactly")

    # ---- host activity (protein-coding genes, as in step 17)
    pc = (G["gene_type"] == "protein_coding").to_numpy()
    h = G["host_tss_sense"].to_numpy(float)
    act = pc & (h >= s17.ACTIVE_MIN)
    lab = np.where(~pc, "-", np.where(h == 0, "silent", np.where(h < s17.ACTIVE_MIN, "low", "active"))).astype(object)
    e3 = np.quantile(h[act], [1 / 3, 2 / 3])
    lab[act] = np.array(["A1", "A2", "A3"])[np.searchsorted(e3, h[act], side="right").clip(0, 2)]
    G.insert(10, "activity", lab)
    print(f"\nprotein-coding genes {int(pc.sum()):,}: " + ", ".join(f"{g_} {int((lab == g_).sum()):,}" for g_ in ACT)
          + f"; active-third edges {e3[0]:.1f} and {e3[1]:.1f} TSS reads; lncRNA genes {int((G['gene_type'] == 'lncRNA').sum()):,}; "
          f"intergenic blocks {len(B):,}")

    print("\n== 1. eligible positions inside H3K36me3 peaks and per fifth (share of eligible nt, %) ==")
    for g_ in ACT + ["intergenic"]:
        d = B if g_ == "intergenic" else G[G["activity"] == g_]
        tot = d["eligible_nt"].sum()
        txt = f"  {g_:10s} genes/blocks {len(d):6,}; eligible {tot / 1e6:7.1f} Mb; in peaks {100 * d['elig_in'].sum() / tot:5.1f}%"
        if g_ != "intergenic":
            txt += "; fifths " + " ".join(f"{100 * d[f'elig_{b}'].sum() / tot:4.1f}" for b in BINS)
        print(txt)
    fin = {g_: G.loc[G["activity"] == g_, "elig_in"].sum() / G.loc[G["activity"] == g_, "eligible_nt"].sum() for g_ in ("silent", "A3")}
    if not fin["A3"] > fin["silent"]:
        print("  WARNING: peaks do not cover more of the most active than of silent genes; check the peak file")

    # ---- deficits
    rows = []
    f1 = lambda v: f"{v:+.1f}" if np.isfinite(v) else "nan"  # noqa: E731

    def run(subset, D, cells, strata, contrasts):
        D = D.reset_index(drop=True)
        cells = np.asarray(cells)
        parts = [np.flatnonzero(cells == cc) for cc in pd.unique(cells)]
        Ws = [s17.boot_W(len(ix), a.n_boot, rng) for ix in parts]
        for tag, T in METRICS:
            res = {}
            for st in strata:
                S, A = D[f"{tag}_{st}_sense"].to_numpy(float), D[f"{tag}_{st}_anti"].to_numpy(float)
                bS = sum(W @ S[ix] for W, ix in zip(Ws, parts))
                bA = sum(W @ A[ix] for W, ix in zip(Ws, parts))
                dr, pt = s17.deficit(bS, bA), float(s17.deficit(S.sum(), A.sum()))
                lo, hi = s17.ci(dr)
                res[st] = (pt, dr)
                lbl = "B" if (subset == "active" and tag == "strong" and st == "all") else ""
                rows.append({"line": a.line, "subset": subset, "stratum": st, "metric": tag, "min_reads": T, "n_units": len(D),
                             "sense": S.sum(), "anti": A.sum(), "deficit": pt, "lo": lo, "hi": hi, "p_boot": np.nan,
                             "label": lbl, "primary": lbl in primary and lbl != ""})
            for hi_, lo_ in contrasts:
                v = res[hi_][1] - res[lo_][1]
                lo, hi = s17.ci(v)
                c_ = f"{hi_}-{lo_}"
                lbl = {"in-out": "A1", "P5-P1": "A2"}.get(c_, "") if (subset == "active" and tag == "strong") else ""
                rows.append({"line": a.line, "subset": subset, "stratum": c_, "metric": tag, "min_reads": T, "n_units": len(D),
                             "sense": np.nan, "anti": np.nan, "deficit": res[hi_][0] - res[lo_][0], "lo": lo, "hi": hi,
                             "p_boot": s17.pval(v), "label": lbl, "primary": lbl in primary and lbl != ""})

    CON = [("in", "out"), ("P5", "P1")]
    Gp = G[G["activity"].isin(ACT)]
    Gact = Gp[Gp["activity"].isin(["A1", "A2", "A3"])]
    run("active", Gact, Gact["activity"].to_numpy(), G_STRATA, CON)
    for g_ in ACT:
        d = Gp[Gp["activity"] == g_]
        run(g_, d, np.full(len(d), g_), G_STRATA, CON)
    run("intergenic", B, np.full(len(B), "blocks"), B_STRATA, [("in", "out")])
    R = pd.DataFrame(rows)

    def show(subset, title):
        x = R[R["subset"] == subset]
        if not len(x):
            return
        tags = [t for t, _ in METRICS]
        print(f"\n-- {title} (n = {int(x['n_units'].iloc[0]):,}) --")
        print("  " + f"{'stratum':9s}" + "".join(f"{t + ' (>=' + str(T) + ')':>26s}" for t, T in METRICS) + f"{'strong sites S / A':>22s}")
        for st in dict.fromkeys(x["stratum"]):
            xs = x[x["stratum"] == st].set_index("metric")
            cells_ = [f"{f1(xs.loc[t, 'deficit'])} [{f1(xs.loc[t, 'lo'])},{f1(xs.loc[t, 'hi'])}]" for t in tags]
            sa = xs.loc["strong"]
            extra = (f"{int(sa['sense']):>11,} /{int(sa['anti']):>8,}" if np.isfinite(sa["sense"])
                     else f"{'p = ' + format(sa['p_boot'], '.2g'):>22s}")
            print("  " + f"{st:9s}" + "".join(f"{cc:>26s}" for cc in cells_) + extra)

    print("\n== 2. gene-strand deficit of measured sites (%), [95% CI]; contrasts with two-sided bootstrap p (strong sites) ==")
    show("active", "active protein-coding genes, thirds A1-A3 pooled (resampled within thirds)")
    for g_ in ACT:
        show(g_, f"protein-coding genes, activity {g_}")
    show("intergenic", "intergenic 1-Mb blocks, plus vs minus strand (control)")

    G.to_csv(a.out_prefix + "_genes.tsv", sep="\t", index=False)
    B.to_csv(a.out_prefix + "_blocks.tsv", sep="\t", index=False)
    R.to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)

    print(f"\n== headline, {a.line}: strong sites (>= {T_strong} reads), active protein-coding genes ==")
    for lbl, st, txt in (("B", "all", "deficit, all eligible positions"), ("A1", "in-out", "inside minus outside H3K36me3 peaks"),
                         ("A2", "P5-P1", "3' fifth minus 5' fifth")):
        r = R[(R["subset"] == "active") & (R["metric"] == "strong") & (R["stratum"] == st)].iloc[0]
        p = f", p = {r['p_boot']:.2g}" if np.isfinite(r["p_boot"]) else ""
        print(f"  {lbl} {'(primary) ' if lbl in primary else '(secondary)'} {txt}: {f1(r['deficit'])} [{f1(r['lo'])}, {f1(r['hi'])}]{p}")
    for subset in ("silent", "intergenic"):
        for st in ("in-out", "P5-P1"):
            x = R[(R["subset"] == subset) & (R["metric"] == "strong") & (R["stratum"] == st)]
            if len(x):
                r = x.iloc[0]
                print(f"  control, {subset}, {st}: {f1(r['deficit'])} [{f1(r['lo'])}, {f1(r['hi'])}], p = {r['p_boot']:.2g}")
    print(f"\nwritten: {a.out_prefix}_genes.tsv, _blocks.tsv, _strata.tsv  [{(time.time() - t0) / 60:.1f} min]")


if __name__ == "__main__":
    main()
