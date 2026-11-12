#!/usr/bin/env python3
"""Step 28b (EXPLORATORY, specified on 10 Nov 2026 after the step-28 results): where along the host gene lie
(1) the weaker sense start of intronic start-site pairs and (2) the gene-strand deficit of strong measured sites?

Motivation. Step 28 found the deficit of strong sites (Results 4) concentrated in the 5' fifth of active genes,
opposite to the prediction from H3K36me3 (contrast 3' fifth minus 5' fifth: K562 -30.0%, HCT116 -41.0%). A quick
look by gene length showed the 5' fifth deficit also in genes > 68 kb, where the 5' fifth spans >= 14 kb. These
analyses are post hoc: no primary contrast, two-sided bootstrap P without correction, reported as exploratory.

Q1  Pairs (Results 1). The depth-matched bias, 100 x (1 - 2^d), d = median DI of intronic pairs minus median DI of
    their intergenic partners (as in step 23b), by the position of the pair's midpoint in the host gene: fifths of
    the annotated gene span (5' -> 3') and distance from the gene's 5' end. Contrast: 5' fifth minus the rest.
    CIs: resamples of host genes, each pair travelling with its intergenic partner. Check: over all pairs the
    estimate must reproduce step 23b (--expect-bias, e.g. 33.0 in K562 and 20.6 in HCT116, to 0.05). The host gene
    is found by its gene_id in the annotation (any gene type) or, if that fails (e.g. an identifier not in GENCODE
    v47), as the only protein-coding gene that contains the pair's midpoint on the strand of the sense start
    (from the arm table, --arms); STOP if fewer than 98% of the pairs are placed.
Q2  Strong measured sites in active protein-coding genes (as in step 28; strong = --strong reads). Deficit,
    100 x (1 - sense / antisense sites), by distance of the site from the gene's 5' end and from the nearest
    annotated TSS of the same gene upstream of the site, and by fifth within thirds of gene length. Contrasts: 5'
    fifth minus the rest; first distance bin minus the rest. CIs: resamples of genes within activity thirds.
    Check: per-gene site counts must equal step 28's table (--genes28); STOP otherwise.
Distance bins (kb from the 5' end or TSS): 1-5, 5-10, 10-20, 20-40, 40-80, >= 80 (eligible positions lie >= 1 kb
from annotated TSSs, so the first bin starts at about 1 kb).

Usage (CPU, about 5 min per data set):
  python3 28b_position.py --line K562 --plus-bw P.bigWig --minus-bw M.bigWig --strong 10 --expect-bias 33.0 \
      --genes28 28_K562_genes.tsv --matched 23b_all_matched_pairs.tsv --arms arms_K562.tsv --chroms chr1,...,chr22 --genome hg38.fa \
      --gtf gencode.v47.annotation.gtf.gz --step10 10_shuffle_test.py --step17 17_observed_sites.py \
      --step28 28_k36_position.py --out-prefix results/28b_K562
Outputs: <prefix>_pairs.tsv (pairs with position), <prefix>_sites.tsv.gz (every strong and >= 10-read site of
protein-coding genes with position and distances), <prefix>_strata.tsv (all estimates and contrasts); log printed.
"""
import argparse
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

BINS_KB = [1, 5, 10, 20, 40, 80]
BIN_LAB = ["1-5", "5-10", "10-20", "20-40", "40-80", ">=80"]
FIFTHS = ["P1", "P2", "P3", "P4", "P5"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def dist_bin(kb):
    return np.clip(np.searchsorted(BINS_KB, kb, side="right") - 1, 0, len(BINS_KB) - 1)


def wmedian(v, w):
    """Median of v with integer weights w, equal to np.median(np.repeat(v, w))."""
    o = np.argsort(v, kind="stable")
    v, c = v[o], np.cumsum(w[o])
    N = int(c[-1]) if len(c) else 0
    if N == 0:
        return np.nan
    if N % 2:
        return float(v[np.searchsorted(c, N // 2, side="right")])
    return 0.5 * float(v[np.searchsorted(c, N // 2 - 1, side="right")] + v[np.searchsorted(c, N // 2, side="right")])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", required=True)
    ap.add_argument("--plus-bw", required=True)
    ap.add_argument("--minus-bw", required=True)
    ap.add_argument("--strong", type=int, required=True, help="reads for a strong site, as in step 28 for this line")
    ap.add_argument("--expect-bias", type=float, required=True, help="step-23b matched bias over all pairs, %")
    ap.add_argument("--genes28", required=True)
    ap.add_argument("--matched", required=True)
    ap.add_argument("--arms", required=True, help="iSTART arm table of this line (chrom, mid, strand_sense)")
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--step17", required=True)
    ap.add_argument("--step28", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--tss-window", type=int, default=250)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    chroms = [c for c in a.chroms.split(",") if c]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    s10 = load_module("s10", os.path.abspath(a.step10))
    s17 = load_module("s17", os.path.abspath(a.step17))
    s28 = load_module("s28", os.path.abspath(a.step28))
    rng = np.random.default_rng(a.seed)
    rows = []
    print(f"== 28b {a.line} (EXPLORATORY): position along the host gene; bootstrap {a.n_boot}, seed {a.seed}")

    G28 = pd.read_csv(a.genes28, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
    need = ["chrom", "gene_id", "gene_name", "gene_type", "strand", "start", "end", "activity",
            "strong_all_sense", "strong_all_anti", "r10_all_sense", "r10_all_anti"]
    miss = [c for c in need if c not in G28.columns]
    if miss:
        sys.exit(f"STOP: {a.genes28} lacks columns {miss}")

    # ------------------------------------------------------------------ Q1: pairs
    M = pd.read_csv(a.matched, sep="\t", dtype={"chrom": str, "ctrl_chrom": str, "host": str})
    miss = [c for c in ("line", "chrom", "mid", "host", "DI", "ctrl_DI") if c not in M.columns]
    if miss:
        sys.exit(f"STOP: {a.matched} lacks columns {miss}")
    M = M[M["line"] == a.line].reset_index(drop=True)
    # gene coordinates for the pairs' chromosomes too (the step-23b pairs may include chromosomes beyond --chroms)
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), set(chroms) | set(M["chrom"]))
    gtss = s17.gene_tss(os.path.abspath(a.gtf), set(chroms))
    pt0 = 100 * (1 - 2 ** (np.median(M["DI"].to_numpy(float)) - np.median(M["ctrl_DI"].to_numpy(float))))
    print(f"\n== Q1. pairs: {len(M):,} matched intronic pairs in {M['host'].nunique():,} host genes; bias over all pairs "
          f"{pt0:.2f}% (step 23b: {a.expect_bias:.1f}%)")
    if abs(pt0 - a.expect_bias) > 0.05 + 1e-9:
        sys.exit("STOP: the overall bias does not reproduce step 23b")
    other = ~M["chrom"].isin(chroms)
    if other.any():
        ex = M.loc[other, "chrom"].value_counts()
        print(f"  NOTE: {int(other.sum()):,} pairs lie outside the chromosomes analysed for single sites: "
              + ", ".join(f"{k} {v:,}" for k, v in ex.items()))
    ARM = pd.read_csv(a.arms, sep="\t", keep_default_na=False, na_values=[""], dtype={"chrom": str})
    miss = [c for c in ("chrom", "mid", "strand_sense") if c not in ARM.columns]
    if miss:
        sys.exit(f"STOP: {a.arms} lacks columns {miss}")
    M = M.merge(ARM[["chrom", "mid", "strand_sense"]].drop_duplicates(["chrom", "mid"]), on=["chrom", "mid"], how="left")
    if M["strand_sense"].isna().any():
        sys.exit(f"STOP: {int(M['strand_sense'].isna().sum())} pairs are missing from {a.arms}")
    by_id = {gid.split(".")[0]: r for gid, r in zip(genes["gene_id"], genes.itertuples(index=False))}
    pcg = genes[genes["gene_type"] == "protein_coding"]
    pc_by = {k: d for k, d in pcg.groupby(["chrom", "strand"])}
    in28 = {g.split(".")[0] for g in G28["gene_id"]}
    recs, how, unplaced = [], dict.fromkeys(("gene_id", "position", "none", "ambiguous"), 0), []
    not28 = 0
    for c, mid, host, st in zip(M["chrom"], M["mid"], M["host"], M["strand_sense"]):
        hid = str(host).split(".")[0]
        r = by_id.get(hid)
        if r is not None and r.chrom == c and r.start <= mid < r.end:
            how["gene_id"] += 1
            not28 += hid not in in28
            recs.append(r)
            continue
        d = pc_by.get((c, st))
        cand = [] if d is None else list(d[(d["start"] <= mid) & (mid < d["end"])].itertuples(index=False))
        if len(cand) == 1:
            how["position"] += 1
            recs.append(cand[0])
        else:
            how["ambiguous" if len(cand) > 1 else "none"] += 1
            unplaced.append(f"{host} ({c}:{mid}, {st}; {'in' if r is not None else 'not in'} the annotation)")
            recs.append(None)
    print(f"  host gene placed by its gene_id {how['gene_id']:,} (of which {not28:,} not in the step-28 table), by "
          f"position and strand of the sense start {how['position']:,}; not placed: no gene {how['none']:,}, several "
          f"genes {how['ambiguous']:,}")
    for u in unplaced[:10]:
        print(f"    not placed: {u}")
    ok = np.array([r is not None for r in recs])
    if ok.mean() < 0.98:
        sys.exit("STOP: fewer than 98% of the pairs could be placed in their host gene")
    Mp = M[ok].reset_index(drop=True)
    R_ = [r for r in recs if r is not None]
    n = np.array([r.end - r.start for r in R_], dtype=float)
    off = np.array([(m - r.start) if r.strand == "+" else (r.end - 1 - m) for m, r in zip(Mp["mid"], R_)], dtype=float)
    Mp["gene_id"] = [r.gene_id for r in R_]
    Mp["gene_len_kb"] = n / 1e3
    Mp["rel"] = off / n
    Mp["fifth"] = np.array(FIFTHS)[np.minimum((Mp["rel"].to_numpy() * 5).astype(int), 4)]
    Mp["dist5_kb"] = off / 1e3
    Mp["dist_bin"] = np.array(BIN_LAB)[dist_bin(Mp["dist5_kb"].to_numpy())]
    Mp.to_csv(a.out_prefix + "_pairs.tsv", sep="\t", index=False)

    hosts, hidx = np.unique(Mp["host"].astype(str).to_numpy(), return_inverse=True)
    W = np.stack([np.bincount(rng.integers(0, len(hosts), len(hosts)), minlength=len(hosts)) for _ in range(a.n_boot)])
    DI, CDI = Mp["DI"].to_numpy(float), Mp["ctrl_DI"].to_numpy(float)

    def bias(mask):
        pt = 100 * (1 - 2 ** (np.median(DI[mask]) - np.median(CDI[mask])))
        w = W[:, hidx[mask]]
        bt = np.array([100 * (1 - 2 ** (wmedian(DI[mask], wb) - wmedian(CDI[mask], wb))) for wb in w])
        return pt, bt

    pt_all, bt_all = bias(np.ones(len(Mp), dtype=bool))
    lo, hi = s17.ci(bt_all)
    print(f"  all placed pairs: {len(Mp):,} in {len(hosts):,} host genes; bias {pt_all:.2f}% ({lo:.1f} to {hi:.1f})")
    auto = Mp["chrom"].isin(chroms).to_numpy()
    if not auto.all():
        pa, ba = bias(auto)
        la, ha = s17.ci(ba)
        rows.append({"line": a.line, "analysis": "pairs", "stratum": "analysed_chromosomes_only", "n": int(auto.sum()),
                     "units": int(len(np.unique(hidx[auto]))), "estimate": pa, "lo": la, "hi": ha, "p_boot": np.nan})
        print(f"  pairs on the analysed chromosomes only: {int(auto.sum()):,}; bias {pa:.2f}% ({la:.1f} to {ha:.1f})")
    rows.append({"line": a.line, "analysis": "pairs", "stratum": "all", "n": len(Mp), "units": len(hosts),
                 "estimate": pt_all, "lo": lo, "hi": hi, "p_boot": np.nan})
    for col, labs in (("fifth", FIFTHS), ("dist_bin", BIN_LAB)):
        res = {}
        print(f"  by {'fifth of the gene (5 -> 3)' if col == 'fifth' else 'distance from the 5′ end (kb)'}:")
        for lab in labs:
            m = (Mp[col] == lab).to_numpy()
            if m.sum() < 10:
                print(f"    {lab:6s} n = {int(m.sum())}: too few pairs")
                continue
            pt, bt = bias(m)
            lo, hi = s17.ci(bt)
            res[lab] = (pt, bt, m)
            rows.append({"line": a.line, "analysis": f"pairs_{col}", "stratum": lab, "n": int(m.sum()),
                         "units": int(len(np.unique(hidx[m]))), "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan})
            print(f"    {lab:6s} n = {int(m.sum()):5,}  bias {pt:+6.1f}% [{lo:+.1f}, {hi:+.1f}]")
        first = labs[0]
        if first in res:
            rest = (Mp[col] != first).to_numpy()
            pr, br = bias(rest)
            v = res[first][1] - br
            lo, hi = s17.ci(v)
            p = s17.pval(v)
            rows.append({"line": a.line, "analysis": f"pairs_{col}", "stratum": f"{first}-rest", "n": len(Mp),
                         "units": len(hosts), "estimate": res[first][0] - pr, "lo": lo, "hi": hi, "p_boot": p})
            print(f"    {first} minus the rest: {res[first][0] - pr:+.1f} [{lo:+.1f}, {hi:+.1f}], p = {p:.2g} "
                  f"(rest: {pr:+.1f}%)")

    # ------------------------------------------------------------------ Q2: strong measured sites
    import pyBigWig
    from pyfaidx import Fasta
    bwP = [pyBigWig.open(p) for p in a.plus_bw.split(",") if p]
    bwM = [pyBigWig.open(p) for p in a.minus_bw.split(",") if p]
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    print(f"\n== Q2. strong measured sites (>= {a.strong} reads; >= 10 reads secondary)  [{time.time() - t0:.0f} s]")
    unit, swap = None, None
    grows = []
    GA = {k: [] for k in ("unit", "sense", "pos", "val", "pk", "bin")}
    for c in chroms:
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        up = np.where(raw >= 97, raw - 32, raw).astype(np.uint8)
        del raw
        L = len(up)
        P = np.zeros(L)
        M_ = np.zeros(L)
        for bw in bwP:
            x, _ = s17.dense(bw, c, L)
            P += np.abs(x)
        for bw in bwM:
            x, _ = s17.dense(bw, c, L)
            M_ += np.abs(x)
        del x
        if unit is None:
            nz = np.concatenate([P[P > 0], M_[M_ > 0]])
            unit = float(nz.min())
            del nz
        g = genes[(genes["chrom"] == c) & (genes["gene_type"] == "protein_coding")]
        if swap is None:
            tS = tA = 0.0
            for r in g.itertuples(index=False):
                sp, sm = s17.window_sum(P, M_, gtss.get((c, r.gene_id), []), a.tss_window, L)
                tS, tA = (tS + sp, tA + sm) if r.strand == "+" else (tS + sm, tA + sp)
            swap = bool(tS / (tS + tA) < 0.5)
            print(f"  strand check on {c}: {100 * tS / (tS + tA):.1f}% of TSS reads on the gene's strand -> "
                  + ("SWAPPING" if swap else "as labelled"))
        if swap:
            P, M_ = M_, P
        ok = up != ord("N")
        del up
        glist, _ = s10.build_masks(c, L, ok, genes, exons, tss)
        vmin = min(10, a.strong) * (1 - 1e-6)
        dummy_pk = None
        for r, s, e, sel in glist:
            if r.gene_type != "protein_coding":
                continue
            ts, ta = (P[s:e], M_[s:e]) if r.strand == "+" else (M_[s:e], P[s:e])
            if dummy_pk is None or len(dummy_pk) < e - s:
                dummy_pk = np.zeros(max(e - s, 1), dtype=bool)
            i = len(grows)
            s28.store(GA, i, 1, ts / unit, sel, s, dummy_pk[:e - s], None, vmin)
            s28.store(GA, i, 0, ta / unit, sel, s, dummy_pk[:e - s], None, vmin)
            grows.append({"chrom": c, "gene_id": r.gene_id, "strand": r.strand, "start": s, "end": e,
                          "tss": sorted(gtss.get((c, r.gene_id), []))})
        print(f"  {c}: {len(grows):,} protein-coding genes so far  [{time.time() - t0:.0f} s]", flush=True)
        del P, M_, ok, glist

    X = s28.concat(GA)
    GR = pd.DataFrame(grows)
    J = GR.merge(G28[need], on=["chrom", "gene_id"], how="left", suffixes=("", "_28"))
    if J["activity"].isna().any():
        sys.exit("STOP: genes missing from the step-28 table")
    site_rows = []
    for tag, T in (("strong", a.strong), ("r10", 10)):
        reps = s28.site_reps(X, T, s10.MERGE)
        u, sflag, pos, val = X["unit"][reps], X["sense"][reps], X["pos"][reps], X["val"][reps]
        cnt_s = np.bincount(u[sflag == 1], minlength=len(GR))
        cnt_a = np.bincount(u[sflag == 0], minlength=len(GR))
        d_s = int((cnt_s != J[f"{tag}_all_sense"].to_numpy()).sum())
        d_a = int((cnt_a != J[f"{tag}_all_anti"].to_numpy()).sum())
        print(f"  check against step 28, {tag}: genes with a different count, sense {d_s:,}, antisense {d_a:,}")
        if d_s or d_a:
            sys.exit("STOP: step 28 is not reproduced")
        st_ = J["start"].to_numpy()[u]
        en_ = J["end"].to_numpy()[u]
        plus = J["strand"].to_numpy()[u] == "+"
        off = np.where(plus, pos - st_, en_ - 1 - pos).astype(float)
        dt = np.empty(len(u))
        tss_l = J["tss"].to_numpy()
        for k in range(len(u)):
            ts_ = np.asarray(tss_l[u[k]], dtype=np.int64)
            up_ = ts_[ts_ <= pos[k]] if plus[k] else ts_[ts_ >= pos[k]]
            dt[k] = (pos[k] - up_.max() if plus[k] else up_.min() - pos[k]) if len(up_) else off[k]
        site_rows.append(pd.DataFrame({"metric": tag, "min_reads": T, "chrom": J["chrom"].to_numpy()[u],
                                       "gene_id": J["gene_id"].to_numpy()[u], "activity": J["activity"].to_numpy()[u],
                                       "gene_len_kb": (en_ - st_) / 1e3, "strand": np.where(sflag == 1, "sense", "anti"),
                                       "pos": pos, "reads": val, "rel": off / (en_ - st_), "dist5_kb": off / 1e3,
                                       "dist_tss_kb": dt / 1e3, "gene_row": u}))
    S = pd.concat(site_rows, ignore_index=True)
    S["fifth"] = np.array(FIFTHS)[np.minimum((S["rel"].to_numpy() * 5).astype(int), 4)]
    S["bin5"] = np.array(BIN_LAB)[dist_bin(S["dist5_kb"].to_numpy())]
    S["bintss"] = np.array(BIN_LAB)[dist_bin(S["dist_tss_kb"].to_numpy())]
    S.drop(columns=["gene_row"]).to_csv(a.out_prefix + "_sites.tsv.gz", sep="\t", index=False, compression="gzip")

    act = J["activity"].isin(["A1", "A2", "A3"]).to_numpy()
    gidx = np.flatnonzero(act)
    pos_in = -np.ones(len(J), dtype=np.int64)
    pos_in[gidx] = np.arange(len(gidx))
    third = J["activity"].to_numpy()[gidx]
    glen = (J["end"] - J["start"]).to_numpy()[gidx] / 1e3
    e3 = np.quantile(glen, [1 / 3, 2 / 3])
    lthird = np.searchsorted(e3, glen, side="right")
    parts = [np.flatnonzero(third == t) for t in ("A1", "A2", "A3")]
    Ws = [s17.boot_W(len(ix), a.n_boot, rng) for ix in parts]
    print(f"  active protein-coding genes {len(gidx):,}; gene-length thirds: < {e3[0]:.0f} kb, {e3[0]:.0f}-{e3[1]:.0f} kb, "
          f">= {e3[1]:.0f} kb")

    for tag in ("strong", "r10"):
        Sx = S[(S["metric"] == tag) & act[S["gene_row"].to_numpy()]]
        g_of = pos_in[Sx["gene_row"].to_numpy()]
        sense = (Sx["strand"] == "sense").to_numpy()

        def est(m, gmask=None):
            cs = np.bincount(g_of[m & sense], minlength=len(gidx)).astype(float)
            ca = np.bincount(g_of[m & ~sense], minlength=len(gidx)).astype(float)
            if gmask is not None:
                cs, ca = cs * gmask, ca * gmask
            bS = sum(Wt @ cs[ix] for Wt, ix in zip(Ws, parts))
            bA = sum(Wt @ ca[ix] for Wt, ix in zip(Ws, parts))
            return float(s17.deficit(cs.sum(), ca.sum())), s17.deficit(bS, bA), int(cs.sum()), int(ca.sum())

        print(f"\n  -- {tag} sites (>= {a.strong if tag == 'strong' else 10} reads), active genes: deficit %, [95% CI], "
              f"sense / antisense sites")
        groups = [("fifth", "fifth of the gene", FIFTHS, None), ("bin5", "kb from the gene's 5' end", BIN_LAB, None),
                  ("bintss", "kb from the nearest upstream TSS of the gene", BIN_LAB, None)]
        for k, nm in enumerate(("short", "middle", "long")):
            groups.append(("fifth", f"fifth of the gene, {nm} third of gene length", FIFTHS, (lthird == k).astype(float)))
        for col, title, labs, gmask in groups:
            print(f"    by {title}:")
            res = {}
            for lab in labs:
                m = (Sx[col] == lab).to_numpy()
                pt, dr, ns, na = est(m, gmask)
                lo, hi = s17.ci(dr)
                res[lab] = (pt, dr)
                rows.append({"line": a.line, "analysis": f"sites_{tag}_{col}" + ("" if gmask is None else f"_{title.split(', ')[1].split()[0]}len"),
                             "stratum": lab, "n": ns + na, "units": int(gmask.sum()) if gmask is not None else len(gidx),
                             "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan, "sense": ns, "anti": na})
                print(f"      {lab:6s} {pt:+6.1f} [{lo:+.1f}, {hi:+.1f}]  {ns:,} / {na:,}")
            m_rest = (Sx[col] != labs[0]).to_numpy()
            pr, drr, _, _ = est(m_rest, gmask)
            v = res[labs[0]][1] - drr
            lo, hi = s17.ci(v)
            p = s17.pval(v)
            rows.append({"line": a.line, "analysis": f"sites_{tag}_{col}" + ("" if gmask is None else f"_{title.split(', ')[1].split()[0]}len"),
                         "stratum": f"{labs[0]}-rest", "n": len(Sx), "units": len(gidx), "estimate": res[labs[0]][0] - pr,
                         "lo": lo, "hi": hi, "p_boot": p, "sense": np.nan, "anti": np.nan})
            print(f"      {labs[0]} minus the rest: {res[labs[0]][0] - pr:+.1f} [{lo:+.1f}, {hi:+.1f}], p = {p:.2g} (rest {pr:+.1f})")

    pd.DataFrame(rows).to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)
    print(f"\nwritten: {a.out_prefix}_pairs.tsv, _sites.tsv.gz, _strata.tsv  [{(time.time() - t0) / 60:.1f} min]")


if __name__ == "__main__":
    main()
