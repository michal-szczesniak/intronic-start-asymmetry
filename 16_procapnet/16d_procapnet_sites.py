#!/usr/bin/env python3
"""Step 16d: potential start sites in the ProCapNet tracks of 16c (real and shuffled), counted
exactly as in step 10, written as step-10 tables so that 13c and 13d run on them unchanged.

  thresholds: top 1e-3 / 3e-4 / 1e-4 of all real predictions on --threshold-chroms
              (both strands; default chr8, chr9 - as for Puffin-D), unless --thresholds;
  usable positions: real prediction present on both strands and not N (as step 10);
  masks, site merging (within 50 nt) and counting: step 10's own functions;
  output: <outdir>/shuffle_<chrom>_{genes,blocks,summary}.tsv with columns real_* and
          shuf1_*, as 13b wrote for Puffin-D.
Also printed: the top-1e-3 quantile of the real predictions per chromosome. Every chromosome
comes from a different held-out model; a model whose outputs run higher or lower shifts how
many sites pass the fixed thresholds (it does not change sense vs antisense within a model).

Usage (CPU):
  python3 16d_procapnet_sites.py --tracks-dir $SCRATCH/16c_tracks --chroms chr1,...,chr22 \
      --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz --step10 ../10_shuffle/10_shuffle_test.py \
      --outdir results/shuffle_pc
"""
import argparse
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_tracks(d, c, tag):
    return tuple(np.load(os.path.join(d, f"{c}.pc_{tag}_{s}.f16.npy")).astype(np.float32) for s in ("plus", "minus"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks-dir", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threshold-chroms", default="chr8,chr9")
    ap.add_argument("--thresholds", default="")
    ap.add_argument("--n-boot", type=int, default=2000)
    a = ap.parse_args()
    s10 = load_module("s10", os.path.abspath(a.step10))
    chroms = [c for c in a.chroms.split(",") if c]
    missing = [c for c in chroms if not os.path.isfile(os.path.join(a.tracks_dir, f"{c}.pc_done.tsv"))]
    if missing:
        sys.exit(f"STOP: 16c not finished for {','.join(missing)}")
    os.makedirs(a.outdir, exist_ok=True)
    tkeys = [f"top{t:g}" for t in s10.TAILS]
    t0 = time.time()

    if a.thresholds:
        thr = dict(zip(s10.TAILS, (float(x) for x in a.thresholds.split(","))))
        src = "given"
    else:
        vals = []
        for c in a.threshold_chroms.split(","):
            for T in load_tracks(a.tracks_dir, c, "real"):
                vals.append(T[np.isfinite(T)])
        allv = np.concatenate(vals)
        thr = {t: float(np.quantile(allv, 1 - t)) for t in s10.TAILS}
        src = f"quantiles of {len(allv):,} real predictions on {a.threshold_chroms}"
        del vals, allv
    print("thresholds (" + src + "): " + ", ".join(f"top {t:g} -> {v:.4g}" for t, v in thr.items()), flush=True)

    from pyfaidx import Fasta
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), set(chroms))
    diag = []
    for c in chroms:
        t1 = time.time()
        P, M = load_tracks(a.tracks_dir, c, "real")
        Ps, Ms = load_tracks(a.tracks_dir, c, "shuf1")
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        if len(raw) != len(P):
            sys.exit(f"STOP: {c} track length {len(P)} != genome length {len(raw)}")
        is_rep = raw >= 97
        up = np.where(is_rep, raw - 32, raw).astype(np.uint8)
        del raw
        L = len(up)
        ok = np.isfinite(P) & np.isfinite(M) & (up != ord("N"))
        glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
        grows = [{"chrom": c, "gene_id": r.gene_id, "gene_name": r.gene_name, "gene_type": r.gene_type,
                  "strand": r.strand, "start": s, "end": e, "eligible_nt": int(sel.sum())} for r, s, e, sel in glist]
        bspans = [(b0, min(L, b0 + s10.BLOCK)) for b0 in range(0, L, s10.BLOCK) if ig[b0:min(L, b0 + s10.BLOCK)].any()]
        brows = [{"chrom": c, "block_start": b0, "eligible_nt": int(ig[b0:b1].sum())} for b0, b1 in bspans]
        n_shuf_missing = 0
        for tag, (TP, TM) in (("real", (P, M)), ("shuf1", (Ps, Ms))):
            for row, (r, s, e, sel) in zip(grows, glist):
                ts, ta = (TP[s:e], TM[s:e]) if r.strand == "+" else (TM[s:e], TP[s:e])
                if tag == "shuf1":
                    n_shuf_missing += int((sel & ~(np.isfinite(ts) & np.isfinite(ta))).sum())
                s10.count_into(row, tag, np.nan_to_num(ts, nan=-np.inf), np.nan_to_num(ta, nan=-np.inf),
                               sel, is_rep[s:e], thr, tkeys)
            for row, (b0, b1) in zip(brows, bspans):
                s10.count_into(row, tag, np.nan_to_num(TP[b0:b1], nan=-np.inf), np.nan_to_num(TM[b0:b1], nan=-np.inf),
                               ig[b0:b1], is_rep[b0:b1], thr, tkeys)
        if n_shuf_missing:
            sys.exit(f"STOP: {c}: {n_shuf_missing} eligible positions lack a shuffled prediction")
        G, B = pd.DataFrame(grows), pd.DataFrame(brows)
        G.to_csv(os.path.join(a.outdir, f"shuffle_{c}_genes.tsv"), sep="\t", index=False)
        B.to_csv(os.path.join(a.outdir, f"shuffle_{c}_blocks.tsv"), sep="\t", index=False)
        rng = np.random.default_rng(1)
        srows = []
        for label, df in (("intron_protein_coding", G[G["gene_type"] == "protein_coding"]),
                          ("intron_lncRNA", G[G["gene_type"] == "lncRNA"]), ("intergenic_plus_vs_minus", B)):
            for k in tkeys:
                So, Ao = df[f"real_all_{k}_sense"], df[f"real_all_{k}_anti"]
                Ss, As = df[f"shuf1_all_{k}_sense"], df[f"shuf1_all_{k}_anti"]
                if len(df) < 2 or min(Ao.sum(), As.sum()) == 0:
                    continue
                (rlo, rhi), (slo, shi), (dlo, dhi) = s10.paired_boot(So, Ao, Ss, As, a.n_boot, rng)
                srows.append({"set": label, "metric": k, "n_units": len(df), "real_sense": So.sum(), "real_anti": Ao.sum(),
                              "real_deficit_pct": 100 * (1 - So.sum() / Ao.sum()), "real_lo95": rlo, "real_hi95": rhi,
                              "shuf_sense": Ss.sum(), "shuf_anti": As.sum(),
                              "shuf_deficit_pct": 100 * (1 - Ss.sum() / As.sum()), "shuf_lo95": slo, "shuf_hi95": shi})
        S = pd.DataFrame(srows)
        S.to_csv(os.path.join(a.outdir, f"shuffle_{c}_summary.tsv"), sep="\t", index=False)
        fin = np.concatenate([P[np.isfinite(P)], M[np.isfinite(M)]])
        q = float(np.quantile(fin, 1 - 1e-3))
        x = S[(S["set"] == "intron_protein_coding") & (S["metric"] == tkeys[0])]
        diag.append({"chrom": c, "top1e-3_quantile": q, "eligible_genes": len(G), "blocks": len(B),
                     "pc_real_deficit_top0.001": float(x["real_deficit_pct"].iloc[0]) if len(x) else np.nan,
                     "pc_shuf_deficit_top0.001": float(x["shuf_deficit_pct"].iloc[0]) if len(x) else np.nan})
        print(f"  {c}: {len(G):,} genes, {len(B):,} blocks; real top-1e-3 quantile {q:.4g}; pc weak-site deficit "
              f"real {diag[-1]['pc_real_deficit_top0.001']:+.1f}%, shuffled {diag[-1]['pc_shuf_deficit_top0.001']:+.1f}%  "
              f"[{time.time() - t1:.0f} s]", flush=True)
        del P, M, Ps, Ms, up, is_rep, ok
    D = pd.DataFrame(diag)
    D.to_csv(os.path.join(a.outdir, "16d_per_chrom.tsv"), sep="\t", index=False)
    print(f"\ntop-1e-3 quantile of real predictions across chromosomes: min {D['top1e-3_quantile'].min():.4g}, "
          f"median {D['top1e-3_quantile'].median():.4g}, max {D['top1e-3_quantile'].max():.4g} "
          f"(threshold for weak sites {thr[s10.TAILS[0]]:.4g})")
    print(f"written: {a.outdir}/shuffle_<chrom>_{{genes,blocks,summary}}.tsv, 16d_per_chrom.tsv  "
          f"[{(time.time() - t0) / 60:.1f} min]")


if __name__ == "__main__":
    main()
