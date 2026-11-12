#!/usr/bin/env python3
"""Step 17: observed K562 PRO-cap initiation in introns, sense vs antisense, as the arbiter
between Puffin-D and ProCapNet.

Steps 13 (Puffin-D) and 16 (ProCapNet) mapped potential start sites over the same eligible
intronic positions and disagree: Puffin-D predicts a sense deficit that grows with the
gene's keto skew (T/G-rich sense strand), ProCapNet a small sense excess that grows with
it. Here the same positions are read in the observed K562 PRO-cap (ENCODE plus- and
minus-strand bigWigs of RNA 5' ends; the data ProCapNet was trained on):

  eligible positions  step 10's own masks (protein-coding / lncRNA gene bodies minus exons
                      +/-100 nt, TSSs +/-1 kb and other genes; intergenic >10 kb from
                      genes), usable = not N: the same positions as ProCapNet (16d)
  observed sites      positions with >= r reads (--min-reads, default 1,3,10), merged
                      within 50 nt per strand (step 10's n_sites); also the read sums
  host activity       observed sense reads within +/- --tss-window nt of the gene's
                      annotated TSSs (outside the counted positions, which exclude TSSs
                      +/-1 kb): silent = 0 reads, low = 1-4, active = 5 or more (tertiles
                      A1-A3). In silent genes the host's own transcription cannot add
                      sense signal
  composition         gene GC and sense-strand keto skew, as in 13d

Deficit = 100 * (1 - sense / antisense), summed over genes (> 0: fewer on the sense strand);
CIs: percentile bootstrap over genes. Printed:
  0. signal checks: values, read unit, strand orientation at annotated TSSs
  1. overall: protein-coding (all / silent / active), lncRNA, intergenic blocks (plus vs
     minus); all and unique (not soft-masked) positions
  2-4. protein-coding genes present in all three tables (observed; Puffin-D, 13b tables;
     ProCapNet, 16d tables), observed next to both models' weak / strong / summed signal:
     2. by host activity
     3. by keto-skew quintile within silent, active and all genes (+ Q5 - Q1)
     4. by GC quintile within silent genes (+ Q5 - Q1)
  5. does a model's per-gene strand prediction (log2 of its summed sense / antisense
     signal) predict the observed per-gene strand bias? Observed deficit by quintile of
     each model's prediction (+ Q5 - Q1; < 0 if the model is informative)
  6. per chromosome: observed deficit vs chromosome GC / keto and vs each model
Caveat: observed initiation needs open chromatin, so this tests the models on realised
start sites, not on all potential ones.

Usage (CPU):
  python3 17_observed_sites.py --plus-bw ENCFF798GNW.bigWig --minus-bw ENCFF662SHP.bigWig \
      --chroms chr1,...,chr22 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --step10 ../10_shuffle/10_shuffle_test.py --puffin-dir ../13_genomewide/results/shuffle \
      --procapnet-dir ../16_procapnet/results/shuffle_pc --out-prefix results/17_obs
"""
import argparse
import gzip
import importlib.util
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HELDOUT = ("chr8", "chr9", "chr10")
ACTIVE_MIN = 5
MODEL_KEYS = (("top0.001", "weak"), ("top0.0001", "strong"), ("sum", "sum"))
RE_ID = re.compile(r'gene_id "([^"]+)"')
LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def counts(codes, lo, hi):
    x = codes[max(0, lo):max(0, hi)]
    return np.bincount(x[x >= 0], minlength=4).astype(float)


def boot_W(n, n_boot, rng):
    idx = rng.integers(0, n, size=(n_boot, n))
    W = np.zeros((n_boot, n), dtype=np.float64)
    for b in range(n_boot):
        W[b] = np.bincount(idx[b], minlength=n)
    return W


def deficit(S, A):
    S, A = np.asarray(S, float), np.asarray(A, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(A > 0, 100 * (1 - S / A), np.nan)


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (np.nan, np.nan)


def pval(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return max(2 * min((v <= 0).mean(), (v >= 0).mean()), 1 / len(v)) if len(v) else np.nan


def gene_tss(path, chroms):
    """(chrom, gene_id) -> annotated TSSs of the gene's transcripts (0-based)."""
    out = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "transcript" or f[0] not in chroms:
                continue
            out.setdefault((f[0], RE_ID.search(f[8]).group(1)), []).append(int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1)
    return out


def dense(bw, c, L):
    """bigWig chromosome as a dense array (0 where no data); also the number of negative intervals."""
    arr = np.zeros(L, dtype=np.float64)
    if c not in bw.chroms():
        return arr, 0
    iv = bw.intervals(c)
    if not iv:
        return arr, 0
    x = np.array(iv, dtype=np.float64)
    s, e, v = x[:, 0].astype(np.int64), x[:, 1].astype(np.int64), x[:, 2]
    n = e - s
    if (n == 1).all():
        arr[s] = v
    else:
        arr[np.repeat(s - np.cumsum(n) + n, n) + np.arange(int(n.sum()))] = np.repeat(v, n)
    return arr, int((v < 0).sum())


def window_sum(P, M, tl, W, L):
    """Reads on each strand over the union of +/- W windows around the TSSs tl."""
    sp = sm = 0.0
    hi_prev = -1
    for t in sorted(tl):
        lo, hi = max(0, t - W, hi_prev), min(L, t + W + 1)
        if hi > lo:
            sp += float(P[lo:hi].sum())
            sm += float(M[lo:hi].sum())
            hi_prev = hi
    return sp, sm


def count_obs(row, n_sites, ts, ta, sel, rep, thr, unit):
    for cls, cm in (("unique", sel & ~rep), ("all", sel)):
        for k, v in thr:
            row[f"obs_{cls}_{k}_sense"] = n_sites(cm & (ts >= v))
            row[f"obs_{cls}_{k}_anti"] = n_sites(cm & (ta >= v))
        row[f"obs_{cls}_sum_sense"] = float(ts[cm].sum()) / unit
        row[f"obs_{cls}_sum_anti"] = float(ta[cm].sum()) / unit


def load_model(d, chroms, pre):
    parts = []
    want = [f"real_all_{k}_{s}" for k, _ in MODEL_KEYS for s in ("sense", "anti")]
    for c in chroms:
        f = os.path.join(d, f"shuffle_{c}_genes.tsv")
        if not os.path.isfile(f):
            sys.exit(f"STOP: missing {f}")
        x = pd.read_csv(f, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
        parts.append(x[["chrom", "gene_id", "eligible_nt"] + want])
    X = pd.concat(parts, ignore_index=True)
    return X.rename(columns={"eligible_nt": f"{pre}_eligible_nt", **{w: f"{pre}_{w[len('real_all_'):]}" for w in want}})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plus-bw", required=True)
    ap.add_argument("--minus-bw", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--puffin-dir", default="")
    ap.add_argument("--procapnet-dir", default="")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-reads", default="1,3,10")
    ap.add_argument("--tss-window", type=int, default=250)
    ap.add_argument("--swap-strands", choices=("auto", "yes", "no"), default="auto",
                    help="auto: swap if most TSS-proximal reads fall on the antisense strand (checked on the first chromosome)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    chroms = [c for c in a.chroms.split(",") if c]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    s10 = load_module("s10", os.path.abspath(a.step10))
    rk = [int(x) for x in a.min_reads.split(",")]
    okeys = [f"r{r}" for r in rk]
    rng = np.random.default_rng(a.seed)
    f1 = lambda v: f"{v:+.1f}" if np.isfinite(v) else "nan"

    import pyBigWig
    bwP, bwM = pyBigWig.open(a.plus_bw), pyBigWig.open(a.minus_bw)
    from pyfaidx import Fasta
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), set(chroms))
    gtss = gene_tss(os.path.abspath(a.gtf), set(chroms))
    print(f"annotation: {len(genes):,} genes, {sum(len(v) for v in gtss.values()):,} transcripts on {len(chroms)} chromosomes  "
          f"[{time.time() - t0:.0f} s]", flush=True)

    unit, swap = None, {"yes": True, "no": False}.get(a.swap_strands)
    grows, brows, crows = [], [], []
    for c in chroms:
        t1 = time.time()
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        is_rep = raw >= 97
        up = np.where(is_rep, raw - 32, raw).astype(np.uint8)
        del raw
        L = len(up)
        for bw, nm in ((bwP, a.plus_bw), (bwM, a.minus_bw)):
            bl = bw.chroms(c)
            if bl is not None and bl != L:
                sys.exit(f"STOP: {c} is {bl:,} nt in {nm} but {L:,} nt in the genome (different assembly?)")
        P, negP = dense(bwP, c, L)
        M, negM = dense(bwM, c, L)
        P, M = np.abs(P), np.abs(M)
        nz = np.concatenate([P[P > 0], M[M > 0]])
        if unit is None:
            if not len(nz):
                sys.exit(f"STOP: no signal on {c} in either bigWig")
            unit = float(nz.min())
            mult = nz / unit
            print(f"\n== 0. signal checks ==\n  {c}: non-zero positions plus {int((P > 0).sum()):,}, minus {int((M > 0).sum()):,}; "
                  f"negative intervals in the bigWigs: plus {negP:,}, minus {negM:,} (absolute values used)")
            print(f"  smallest non-zero value {unit:.6g}; values that are whole numbers: {100 * float((np.abs(nz - np.round(nz)) < 1e-6).mean()):.2f}%; "
                  f"whole multiples of the smallest: {100 * float((np.abs(mult - np.round(mult)) < 1e-3).mean()):.2f}%")
            print(f"  non-zero values in reads (= multiples of {unit:.4g}), quantiles 50/90/99/99.9%: "
                  + ", ".join(f"{np.quantile(mult, q):.4g}" for q in (0.5, 0.9, 0.99, 0.999)) + f"; max {mult.max():.4g}")
            print(f"  read unit = {unit:.6g}; site thresholds: " + ", ".join(f">= {r} reads" for r in rk))
        elif len(nz) and nz.min() < unit * (1 - 1e-6):
            print(f"  WARNING: {c}: smallest non-zero value {nz.min():.6g} is below the read unit {unit:.6g}")
        thr = [(k, r * unit * (1 - 1e-6)) for k, r in zip(okeys, rk)]

        g = genes[(genes["chrom"] == c) & genes["gene_type"].isin(("protein_coding", "lncRNA"))]
        hs = {r.gene_id: (window_sum(P, M, gtss.get((c, r.gene_id), []), a.tss_window, L), r.strand)
              for r in g.itertuples(index=False)}
        hs = {k: ((sp, sm) if st == "+" else (sm, sp)) for k, ((sp, sm), st) in hs.items()}
        pcids = set(g.loc[g["gene_type"] == "protein_coding", "gene_id"])
        tS = sum(v[0] for k, v in hs.items() if k in pcids)
        tA = sum(v[1] for k, v in hs.items() if k in pcids)
        share = tS / (tS + tA) if tS + tA > 0 else np.nan
        if swap is None:
            swap = bool(share < 0.5)
            print(f"  strand check on {c}: reads within +/-{a.tss_window} nt of annotated protein-coding TSSs fall "
                  f"{100 * share:.1f}% on the gene's strand -> " + ("SWAPPING the plus and minus bigWigs" if swap else "orientation as labelled"))
        if swap:
            P, M = M, P
            hs = {k: (v[1], v[0]) for k, v in hs.items()}
            share = 1 - share
        if share < 0.5:
            print(f"  WARNING: {c}: only {100 * share:.1f}% of TSS-proximal reads on the gene's strand")

        ok = up != ord("N")
        glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
        codes = LUT[up]
        for r, s, e, sel in glist:
            ts, ta = (P[s:e], M[s:e]) if r.strand == "+" else (M[s:e], P[s:e])
            hS, hA = hs.get(r.gene_id, (np.nan, np.nan))
            m = counts(codes, s, e)
            A_, C_, G_, T_ = m if r.strand == "+" else m[::-1]
            tot = m.sum()
            row = {"chrom": c, "gene_id": r.gene_id, "gene_name": r.gene_name, "gene_type": r.gene_type,
                   "strand": r.strand, "start": s, "end": e, "eligible_nt": int(sel.sum()),
                   "host_tss_sense": hS / unit, "host_tss_anti": hA / unit,
                   "GC": (G_ + C_) / tot if tot else np.nan, "keto": (G_ + T_ - A_ - C_) / tot if tot else np.nan}
            count_obs(row, s10.n_sites, ts, ta, sel, is_rep[s:e], thr, unit)
            grows.append(row)
        for b0 in range(0, L, s10.BLOCK):
            b1 = min(L, b0 + s10.BLOCK)
            if ig[b0:b1].any():
                row = {"chrom": c, "block_start": b0, "eligible_nt": int(ig[b0:b1].sum())}
                count_obs(row, s10.n_sites, P[b0:b1], M[b0:b1], ig[b0:b1], is_rep[b0:b1], thr, unit)
                brows.append(row)
        crows.append({"chrom": c, "reads_plus": float(P.sum()) / unit, "reads_minus": float(M.sum()) / unit,
                      "tss_sense_share": share, "genes": len(glist)})
        print(f"  {c}: {len(glist):,} genes; reads plus {crows[-1]['reads_plus']:,.0f}, minus {crows[-1]['reads_minus']:,.0f}; "
              f"TSS reads on the gene's strand {100 * share:.1f}%  [{time.time() - t1:.0f} s]", flush=True)
        del P, M, up, is_rep, ok, codes
    G, B, CR = pd.DataFrame(grows), pd.DataFrame(brows), pd.DataFrame(crows)
    G.to_csv(a.out_prefix + "_genes.tsv", sep="\t", index=False)
    B.to_csv(a.out_prefix + "_blocks.tsv", sep="\t", index=False)
    print(f"  genome: reads plus {CR['reads_plus'].sum():,.0f}, minus {CR['reads_minus'].sum():,.0f}; plus/minus strands "
          f"{'SWAPPED' if swap else 'as labelled'}  [{(time.time() - t0) / 60:.1f} min]", flush=True)

    rows = []
    OBS = [(f"obs {k}", f"obs_all_{k}") for k in okeys] + [("obs sum", "obs_all_sum")]
    OBS_U = [(f"obs {k}", f"obs_unique_{k}") for k in okeys] + [("obs sum", "obs_unique_sum")]

    def run_bins(df, lab, order, metrics, analysis, subset, med=None, contrast=None):
        res = {}
        for b in order:
            d = df[lab == b]
            if len(d) < 2:
                continue
            W = boot_W(len(d), a.n_boot, rng)
            for ml, col in metrics:
                S, A = d[f"{col}_sense"].to_numpy(float), d[f"{col}_anti"].to_numpy(float)
                dr, pt = deficit(W @ S, W @ A), float(deficit(S.sum(), A.sum()))
                res[(b, ml)] = (pt, dr)
                lo, hi = ci(dr)
                rows.append({"analysis": analysis, "subset": subset, "bin": b, "metric": ml, "n_units": len(d),
                             "sense": S.sum(), "anti": A.sum(), "deficit": pt, "lo": lo, "hi": hi,
                             "median": float(d[med].median()) if med else np.nan})
        if contrast:
            bh, bl = contrast
            for ml, _ in metrics:
                if (bh, ml) in res and (bl, ml) in res:
                    v = res[(bh, ml)][1] - res[(bl, ml)][1]
                    lo, hi = ci(v)
                    rows.append({"analysis": analysis, "subset": subset, "bin": f"{bh}-{bl}", "metric": ml,
                                 "deficit": res[(bh, ml)][0] - res[(bl, ml)][0], "lo": lo, "hi": hi, "p_boot": pval(v)})
        return res

    def table(analysis, subset, labels, title, med=None):
        R = pd.DataFrame(rows)
        x = R[(R["analysis"] == analysis) & (R["subset"] == subset)]
        if not len(x):
            return
        wd = 23
        print(f"\n-- {title} --")
        print("  " + f"{'bin':14s}{'n':>7s}" + (f"{med:>9s}" if med else "") + "".join(f"{lb:>{wd}s}" for lb in labels))
        for b in dict.fromkeys(x["bin"]):
            xb = x[x["bin"] == b].set_index("metric")
            n = xb["n_units"].iloc[0]
            cells = [f"{f1(xb.loc[lb, 'deficit'])}[{f1(xb.loc[lb, 'lo'])},{f1(xb.loc[lb, 'hi'])}]" if lb in xb.index else "-"
                     for lb in labels]
            mtxt = (f"{xb['median'].iloc[0]:9.3f}" if np.isfinite(xb["median"].iloc[0]) else " " * 9) if med else ""
            print("  " + f"{b:14s}" + (f"{int(n):7,d}" if np.isfinite(n) else " " * 7) + mtxt + "".join(f"{cc:>{wd}s}" for cc in cells))
            if not np.isfinite(n):
                print("  " + " " * (21 + (9 if med else 0)) + "".join(
                    f"{'p=' + format(xb.loc[lb, 'p_boot'], '.2g'):>{wd}s}" if lb in xb.index else " " * wd for lb in labels))

    # ---- 1. overall
    P_all = G[(G["gene_type"] == "protein_coding") & G["GC"].notna()].reset_index(drop=True)
    host = P_all["host_tss_sense"].to_numpy(float)
    print(f"\nprotein-coding genes {len(P_all):,}: silent (0 TSS reads) {int((host == 0).sum()):,}, low (<{ACTIVE_MIN}) "
          f"{int(((host > 0) & (host < ACTIVE_MIN)).sum()):,}, active {int((host >= ACTIVE_MIN).sum()):,}; "
          f"lncRNA genes {int((G['gene_type'] == 'lncRNA').sum()):,}; intergenic blocks {len(B):,}")
    sets = (("pc_all", P_all), ("pc_silent", P_all[host == 0]), ("pc_active", P_all[host >= ACTIVE_MIN]),
            ("lncRNA", G[G["gene_type"] == "lncRNA"]), ("intergenic", B))
    for cls, mets in (("all", OBS), ("unique", OBS_U)):
        for nm, df in sets:
            run_bins(df, np.full(len(df), nm), [nm], mets, "overall", cls)
    print("\n== 1. observed sense deficit (%), all sets (intergenic: plus vs minus) ==")
    for cls in ("all", "unique"):
        table("overall", cls, [ml for ml, _ in OBS], f"{cls} positions")
    R = pd.DataFrame(rows)
    x = R[(R["analysis"] == "overall") & (R["subset"] == "all") & (R["metric"] == f"obs {okeys[0]}")]
    print("  sites (sense / antisense), " + okeys[0] + ": " + "; ".join(
        f"{b} {int(s):,} / {int(an):,}" for b, s, an in zip(x["bin"], x["sense"], x["anti"])))

    # ---- joined with the models
    J = P_all.copy()
    MOD = []
    print("\njoining with the model tables (protein-coding genes):")
    for pre, name, d in (("pf", "Puffin-D", a.puffin_dir), ("pc", "ProCapNet", a.procapnet_dir)):
        if not d:
            continue
        X = load_model(d, chroms, pre)
        n0 = len(J)
        J = J.merge(X, on=["chrom", "gene_id"], how="inner")
        print(f"  {name}: {len(X):,} genes in its tables; joined {len(J):,} of {n0:,} protein-coding genes; "
              f"median eligible nt {name} / observed {float((J[f'{pre}_eligible_nt'] / J['eligible_nt']).median()):.3f}")
        MOD += [(f"{name} {lb}", f"{pre}_{k}") for k, lb in MODEL_KEYS]
    J = J.reset_index(drop=True)
    jh = J["host_tss_sense"].to_numpy(float)
    OL = [ml for ml, _ in OBS]
    ML = [ml for ml, _ in MOD]
    ALL = OBS + MOD

    def both_tables(analysis, subset, title, med=None):
        table(analysis, subset, OL, title + ": observed", med)
        if ML:
            table(analysis, subset, ML, title + ": models", med)

    # ---- 2. host activity
    lab = np.where(np.isnan(jh), "-", np.where(jh == 0, "silent", np.where(jh < ACTIVE_MIN, "low", "active")))
    act = jh >= ACTIVE_MIN
    if act.sum() >= 3:
        e = np.quantile(jh[act], [1 / 3, 2 / 3])
        lab = np.where(act, np.array(["A1", "A2", "A3"])[np.searchsorted(e, jh, side="right").clip(0, 2)], lab)
    order = ["silent", "low", "A1", "A2", "A3"]
    run_bins(J, lab, order, ALL, "host", "pc", "host_tss_sense", ("A3", "silent"))
    print(f"\n== 2. protein-coding genes by host activity (sense reads within +/-{a.tss_window} nt of TSSs; "
          f"median shown); sense deficit (%) ==")
    both_tables("host", "pc", "by host activity", "host")
    kmed = pd.Series(J["keto"].to_numpy(float)).groupby(lab).median()
    print("  median keto skew: " + ", ".join(f"{b} {kmed[b]:+.4f}" for b in order if b in kmed.index))

    # ---- 3. keto quintiles, 4. GC quintiles within silent
    def quint(col, sub_mask):
        e = np.quantile(J[col].to_numpy(float), [0.2, 0.4, 0.6, 0.8])
        q = np.searchsorted(e, J[col].to_numpy(float), side="right")
        return np.where(sub_mask, np.array(["Q1", "Q2", "Q3", "Q4", "Q5"])[q], "-")

    Q = ["Q1", "Q2", "Q3", "Q4", "Q5"]
    print("\n== 3. protein-coding genes by sense-strand keto-skew quintile (edges from all joined genes); sense deficit (%) ==")
    for sub, msk in (("silent", jh == 0), ("active", act), ("all", np.ones(len(J), bool))):
        run_bins(J, quint("keto", msk), Q, ALL, "keto", sub, "keto", ("Q5", "Q1"))
        both_tables("keto", sub, f"keto quintiles, {sub} genes", "keto")
    print("\n== 4. silent protein-coding genes by GC quintile (edges from all joined genes); sense deficit (%) ==")
    run_bins(J, quint("GC", jh == 0), Q, ALL, "GC", "silent", "GC", ("Q5", "Q1"))
    both_tables("GC", "silent", "GC quintiles, silent genes", "GC")

    # ---- 5. model per-gene strand prediction vs observed
    if MOD:
        print("\n== 5. observed deficit by quintile of each model's per-gene strand prediction "
              "(log2 summed sense / antisense; Q5 = most sense-favoured) ==")
        for pre, name in (("pf", "Puffin-D"), ("pc", "ProCapNet")):
            if f"{pre}_sum_sense" not in J:
                continue
            J[f"{pre}_bias"] = np.log2((J[f"{pre}_sum_sense"] + 1e-9) / (J[f"{pre}_sum_anti"] + 1e-9))
            for sub, msk in (("silent", jh == 0), ("all", np.ones(len(J), bool))):
                d = J[msk]
                e = np.quantile(d[f"{pre}_bias"], [0.2, 0.4, 0.6, 0.8])
                q = np.where(msk, np.array(Q)[np.searchsorted(e, J[f"{pre}_bias"].to_numpy(float), side="right")], "-")
                mets = OBS + [(ml, col) for ml, col in MOD if col.startswith(pre + "_")]
                run_bins(J, q, Q, mets, f"by_{pre}_bias", sub, f"{pre}_bias", ("Q5", "Q1"))
                table(f"by_{pre}_bias", sub, [ml for ml, _ in mets], f"{name} prediction quintiles, {sub} genes", "bias")
                for k in (okeys[0], "sum"):
                    s_, a_ = d[f"obs_all_{k}_sense"].to_numpy(float), d[f"obs_all_{k}_anti"].to_numpy(float)
                    has = (s_ + a_) >= 1
                    if has.sum() >= 10:
                        rho, p = spearmanr(d.loc[has, f"{pre}_bias"], np.log2((s_[has] + 0.5) / (a_[has] + 0.5)))
                        print(f"  per gene, {sub}: Spearman({name} log2 sense/antisense, observed log2 (sense+0.5)/(antisense+0.5) "
                              f"of {k}) = {rho:+.3f} (p = {p:.2g}, {int(has.sum()):,} genes with signal)")

    # ---- 6. per chromosome
    crow = []
    for c in chroms:
        d = J[J["chrom"] == c]
        if not len(d):
            continue
        rr = {"chrom": c, "genes": len(d), "silent_genes": int((d["host_tss_sense"] == 0).sum()),
              "GC_mean": float(d["GC"].mean()), "keto_mean": float(d["keto"].mean())}
        for ml, col in ALL:
            if ml.endswith(" sum") and ml != "obs sum":
                continue
            rr[ml] = float(deficit(d[f"{col}_sense"].sum(), d[f"{col}_anti"].sum()))
        crow.append(rr)
    C = pd.DataFrame(crow).merge(CR[["chrom", "reads_plus", "reads_minus", "tss_sense_share"]], on="chrom", how="left")
    C.to_csv(a.out_prefix + "_chroms.tsv", sep="\t", index=False)
    print("\n== 6. per chromosome, joined protein-coding genes; sense deficit (%) ==")
    print(C.drop(columns=["reads_plus", "reads_minus"]).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    ref = f"obs {okeys[min(1, len(okeys) - 1)]}"            # r3: little host background (section 2)
    if len(C) >= 5:
        for col in ["GC_mean", "keto_mean"] + [ml for ml in C.columns if ml.endswith(" weak") or ml.endswith(" strong")]:
            ok_ = C[[ref, col]].dropna()
            rho, p = spearmanr(ok_[ref], ok_[col])
            print(f"  Spearman across chromosomes, {ref} deficit vs {col}: {rho:+.3f} (p = {p:.2g})")
        for col in [ref] + [ml for ml in C.columns if ml.endswith(" weak")]:
            rk_ = C[col].rank(ascending=False)
            hold = ", ".join(f"{h} {int(rk_[C['chrom'] == h].iloc[0])}" for h in HELDOUT if (C["chrom"] == h).any())
            print(f"  rank (1 = largest deficit of {len(C)}) by {col}: {hold}")

    R = pd.DataFrame(rows)
    R.to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)

    # ---- headline
    print("\n== headline: sense deficit (%) [95% CI] ==")

    def pick(an, sub, b, ml):
        x = R[(R["analysis"] == an) & (R["subset"] == sub) & (R["bin"] == b) & (R["metric"] == ml)]
        if not len(x):
            return "-"
        r = x.iloc[0]
        return f"{f1(r['deficit'])} [{f1(r['lo'])},{f1(r['hi'])}]" + (f" p={r['p_boot']:.2g}" if "-" in b else "")

    for title, an, sub, b in (("silent genes", "host", "pc", "silent"), ("active genes, top host tertile", "host", "pc", "A3"),
                              ("keto Q5-Q1, silent genes", "keto", "silent", "Q5-Q1"),
                              ("keto Q5-Q1, active genes", "keto", "active", "Q5-Q1"),
                              ("keto Q5-Q1, all genes", "keto", "all", "Q5-Q1")):
        print(f"  {title}:")
        print("    " + ";  ".join(f"{ml} {pick(an, sub, b, ml)}" for ml in OL))
        if ML:
            print("    " + ";  ".join(f"{ml} {pick(an, sub, b, ml)}" for ml in ML))
    print(f"\nCIs: percentile bootstrap over genes (blocks), {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out_prefix}_genes.tsv, _blocks.tsv, _strata.tsv, _chroms.tsv  [{(time.time() - t0) / 60:.1f} min]")


if __name__ == "__main__":
    main()
