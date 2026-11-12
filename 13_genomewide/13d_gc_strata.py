#!/usr/bin/env python3
"""Step 13d: does the sense deficit of potential start sites depend on regional composition?

Step 13c: genome-wide, medium and strong sites show a sense deficit on every
chromosome, but weak sites (top 1e-3) vary from -12% to +13% between chromosomes, the
negative ones being the GC-rich, gene-dense chromosomes; chr8, 9, 10 (the held-out
group) happen to be among the highest. This step asks whether gene composition explains
that, using the same per-gene tables (results/shuffle/shuffle_<chrom>_genes.tsv) and
the genome:

  per gene (protein-coding and lncRNA; whole gene body, read on the gene's sense strand):
    GC      fraction G+C
    keto    keto skew ((G+T) - (A+C)) / (A+C+G+T) - the T- and G-richness that
            transcription-coupled mutation leaves on the sense strand
    GC_reg  fraction G+C of the gene +/- --flank (isochore-level)
  per intergenic 1 Mb block: GC and plus-strand keto skew (control).

Readouts, protein-coding introns, as step 13c (real deficit, shuffled deficit, drop =
real - shuffled; deficit = 100 * (1 - total sense / total antisense); paired bootstrap
over genes):
  1. by GC quintile (and Q5 - Q1 contrasts), for each site strength
  2. by keto-skew quintile
  3. GC tertile x keto tertile, weak sites
  4. held-out (chr8, 9, 10) vs training chromosomes WITHIN GC quintiles, weak sites:
     if composition explains the 13c group difference, it should vanish here
  5. per chromosome: mean GC and keto skew of genes vs the chromosome's weak-site
     deficit (Spearman across chromosomes)
  6. intergenic blocks by GC quintile (should stay ~0)
  7. lncRNA introns by GC quintile, weak sites
The gene body includes exons, but the counted positions are intronic (exons +/-100 nt
and TSSs +/-1 kb excluded), so gene-body composition stands in for the counted introns.

Usage:
  python3 13d_gc_strata.py --dir results/shuffle --chroms chr1,...,chr22 --genome hg38.fa \
      --out-prefix results/13d_gc [--flank 50000] [--n-boot 2000] [--seed 1]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

TKEYS = ("top0.001", "top0.0003", "top0.0001")
HELDOUT = ("chr8", "chr9", "chr10")
LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i


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


def group_stats(d, m, W):
    """Real / shuffled deficit and drop for rows d (one group), metric m, with bootstrap draws."""
    So, Ao = d[f"real_all_{m}_sense"].to_numpy(float), d[f"real_all_{m}_anti"].to_numpy(float)
    Ss, As = d[f"shuf1_all_{m}_sense"].to_numpy(float), d[f"shuf1_all_{m}_anti"].to_numpy(float)
    br, bs = deficit(W @ So, W @ Ao), deficit(W @ Ss, W @ As)
    r, s = float(deficit(So.sum(), Ao.sum())), float(deficit(Ss.sum(), As.sum()))
    return {"n_units": len(d), "sites": int(So.sum() + Ao.sum()), "real": r, "real_lo": ci(br)[0], "real_hi": ci(br)[1],
            "shuf": s, "shuf_lo": ci(bs)[0], "shuf_hi": ci(bs)[1], "drop": r - s, "drop_lo": ci(br - bs)[0],
            "drop_hi": ci(br - bs)[1]}, (br, bs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    chroms = [c for c in a.chroms.split(",") if c]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    rng = np.random.default_rng(a.seed)
    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)

    G, B = [], []
    for c in chroms:
        fg, fb = (os.path.join(a.dir, f"shuffle_{c}_{k}.tsv") for k in ("genes", "blocks"))
        if not (os.path.isfile(fg) and os.path.isfile(fb)):
            sys.exit(f"STOP: missing tables for {c}")
        g = pd.read_csv(fg, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
        b = pd.read_csv(fb, sep="\t")
        codes = LUT[np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)]
        L = len(codes)
        feats = []
        for r in g.itertuples(index=False):
            m = counts(codes, int(r.start), int(r.end))
            rg = counts(codes, int(r.start) - a.flank, min(L, int(r.end) + a.flank))
            A_, C_, G_, T_ = m if r.strand == "+" else m[::-1]          # sense strand: complement swaps A<->T, C<->G
            tot = m.sum()
            feats.append((((G_ + C_) / tot) if tot else np.nan, ((G_ + T_ - A_ - C_) / tot) if tot else np.nan,
                          (rg[1] + rg[2]) / rg.sum() if rg.sum() else np.nan))
        g[["GC", "keto", "GC_reg"]] = np.array(feats, dtype=float).reshape(-1, 3)
        bf = []
        for r in b.itertuples(index=False):
            m = counts(codes, int(r.block_start), min(L, int(r.block_start) + 1_000_000))
            tot = m.sum()
            bf.append((((m[1] + m[2]) / tot) if tot else np.nan, ((m[2] + m[3] - m[0] - m[1]) / tot) if tot else np.nan))
        b[["GC", "keto"]] = np.array(bf, dtype=float).reshape(-1, 2)
        G.append(g)
        B.append(b)
        print(f"  {c}: {len(g):,} genes, {len(b):,} blocks", flush=True)
    G = pd.concat(G, ignore_index=True)
    B = pd.concat(B, ignore_index=True)
    P = G[(G["gene_type"] == "protein_coding") & G["GC"].notna()].reset_index(drop=True)
    N = G[(G["gene_type"] == "lncRNA") & G["GC"].notna()].reset_index(drop=True)
    B = B[B["GC"].notna()].reset_index(drop=True)
    print(f"\nprotein-coding genes {len(P):,}, lncRNA genes {len(N):,}, intergenic blocks {len(B):,}")
    print(f"protein-coding genes: Spearman(GC, keto skew) = {spearmanr(P['GC'], P['keto'])[0]:.3f}; "
          f"Spearman(GC, GC +/- {a.flank // 1000} kb) = {spearmanr(P['GC'], P['GC_reg'])[0]:.3f}; "
          f"median keto skew {P['keto'].median():+.4f} (intergenic blocks, plus strand: {B['keto'].median():+.4f})")

    rows = []

    def binned(df, col, k, label, metrics, unit):
        edges = np.quantile(df[col], np.linspace(0, 1, k + 1)[1:-1])
        q = np.searchsorted(edges, df[col].to_numpy(), side="right")
        draws, stats = {}, {}
        for qi in range(k):
            d = df[q == qi]
            W = boot_W(len(d), a.n_boot, rng)
            for m in metrics:
                st, dr = group_stats(d, m, W)
                draws[(qi, m)], stats[(qi, m)] = dr, st
                rows.append({"analysis": label, "unit": unit, "bin": f"Q{qi + 1}", "metric": m,
                             f"{col}_median": float(d[col].median()), **st})
        for m in metrics:
            (r1, s1), (rk, sk) = draws[(0, m)], draws[(k - 1, m)]
            for qn, v1, vk in (("real", r1, rk), ("shuf", s1, sk), ("drop", r1 - s1, rk - sk)):
                val = stats[(k - 1, m)][qn] - stats[(0, m)][qn]
                rows.append({"analysis": label, "unit": unit, "bin": f"Q{k}-Q1 {qn}", "metric": m,
                             "real": val, "real_lo": ci(vk - v1)[0], "real_hi": ci(vk - v1)[1], "p_boot": pval(vk - v1)})
        return edges

    e_gc = binned(P, "GC", 5, "pc_by_GC", TKEYS, "gene")
    binned(P, "keto", 5, "pc_by_keto", TKEYS, "gene")
    binned(P, "GC_reg", 5, "pc_by_GC_regional", TKEYS[:1], "gene")
    binned(B, "GC", 5, "intergenic_by_GC", TKEYS[:1], "1Mb_block")
    binned(N, "GC", 5, "lnc_by_GC", TKEYS[:1], "gene")

    # GC tertile x keto tertile, weak sites
    t_gc = np.searchsorted(np.quantile(P["GC"], [1 / 3, 2 / 3]), P["GC"].to_numpy(), side="right")
    t_kt = np.searchsorted(np.quantile(P["keto"], [1 / 3, 2 / 3]), P["keto"].to_numpy(), side="right")
    for i in range(3):
        for j in range(3):
            d = P[(t_gc == i) & (t_kt == j)]
            if len(d) < 30:
                continue
            st, _ = group_stats(d, TKEYS[0], boot_W(len(d), a.n_boot, rng))
            rows.append({"analysis": "pc_GCtertile_x_ketotertile", "unit": "gene", "bin": f"GC{i + 1}_keto{j + 1}",
                         "metric": TKEYS[0], **st})

    # held-out vs training within GC quintiles, weak sites
    q = np.searchsorted(e_gc, P["GC"].to_numpy(), side="right")
    held = P["chrom"].isin(HELDOUT).to_numpy()
    for qi in range(5):
        res = {}
        for name, msk in (("heldout", held), ("training", ~held)):
            d = P[(q == qi) & msk]
            st, dr = group_stats(d, TKEYS[0], boot_W(len(d), a.n_boot, rng))
            res[name] = (st, dr)
            rows.append({"analysis": "pc_heldout_vs_training_by_GC", "unit": "gene", "bin": f"Q{qi + 1} {name}",
                         "metric": TKEYS[0], **st})
        dv = res["training"][1][0] - res["heldout"][1][0]
        rows.append({"analysis": "pc_heldout_vs_training_by_GC", "unit": "gene", "bin": f"Q{qi + 1} training-heldout real",
                     "metric": TKEYS[0], "real": res["training"][0]["real"] - res["heldout"][0]["real"],
                     "real_lo": ci(dv)[0], "real_hi": ci(dv)[1], "p_boot": pval(dv)})

    # per chromosome
    crow = []
    for c in chroms:
        d = P[P["chrom"] == c]
        if not len(d):
            continue
        So, Ao = d["real_all_top0.001_sense"].sum(), d["real_all_top0.001_anti"].sum()
        Ss, As = d["shuf1_all_top0.001_sense"].sum(), d["shuf1_all_top0.001_anti"].sum()
        w = d["real_all_top0.001_sense"] + d["real_all_top0.001_anti"]
        crow.append({"chrom": c, "genes": len(d), "GC_mean": float(d["GC"].mean()),
                     "GC_site_weighted": float(np.average(d["GC"], weights=w)) if w.sum() else np.nan,
                     "keto_mean": float(d["keto"].mean()), "real_deficit": float(deficit(So, Ao)),
                     "shuf_deficit": float(deficit(Ss, As))})
    C = pd.DataFrame(crow)
    R = pd.DataFrame(rows)
    R.to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)
    C.to_csv(a.out_prefix + "_chroms.tsv", sep="\t", index=False)

    f1 = lambda v: f"{v:+.1f}" if v == v else "nan"
    show = lambda r: (f"real {f1(r['real'])} [{f1(r['real_lo'])}, {f1(r['real_hi'])}]   shuffled {f1(r['shuf'])} "
                      f"[{f1(r['shuf_lo'])}, {f1(r['shuf_hi'])}]   drop {f1(r['drop'])} [{f1(r['drop_lo'])}, {f1(r['drop_hi'])}]")
    for label, col, title, metrics in (("pc_by_GC", "GC", "protein-coding introns by gene GC quintile", TKEYS),
                                       ("pc_by_keto", "keto", "protein-coding introns by sense-strand keto skew quintile", TKEYS),
                                       ("pc_by_GC_regional", "GC_reg", f"protein-coding introns by regional GC (gene +/- {a.flank // 1000} kb)", TKEYS[:1]),
                                       ("intergenic_by_GC", "GC", "intergenic control (plus vs minus) by block GC quintile", TKEYS[:1]),
                                       ("lnc_by_GC", "GC", "lncRNA introns by gene GC quintile", TKEYS[:1])):
        print(f"\n== {title}; sense deficit (%) ==")
        for m in metrics:
            print(f"  {m}")
            x = R[(R["analysis"] == label) & (R["metric"] == m)]
            for r in x.to_dict("records"):
                if r["bin"].startswith("Q") and " " not in r["bin"]:
                    print(f"    {r['bin']} {col} {r[col + '_median']:.3f}  n={int(r['n_units']):6,d}  sites {int(r['sites']):7,d}  {show(r)}")
                else:
                    print(f"    {r['bin']:14s} {f1(r['real'])} [{f1(r['real_lo'])}, {f1(r['real_hi'])}] p={r['p_boot']:.3g}")
    print("\n== GC tertile x keto tertile, protein-coding introns, top0.001: real deficit [CI] (shuffled) ==")
    x = R[R["analysis"] == "pc_GCtertile_x_ketotertile"].set_index("bin")
    print("            keto low                       keto mid                       keto high")
    for i in range(3):
        cells = []
        for j in range(3):
            k = f"GC{i + 1}_keto{j + 1}"
            cells.append(f"{f1(x.loc[k, 'real'])} [{f1(x.loc[k, 'real_lo'])},{f1(x.loc[k, 'real_hi'])}] ({f1(x.loc[k, 'shuf'])})"
                         if k in x.index else "-")
        print(f"  GC {['low ', 'mid ', 'high'][i]}  " + "   ".join(f"{c:30s}" for c in cells))
    print("\n== held-out (chr8,9,10) vs training chromosomes within gene-GC quintiles, top0.001 ==")
    for r in R[R["analysis"] == "pc_heldout_vs_training_by_GC"].to_dict("records"):
        if r["bin"].endswith("real"):
            print(f"    {r['bin']:30s} {f1(r['real'])} [{f1(r['real_lo'])}, {f1(r['real_hi'])}] p={r['p_boot']:.3g}")
        else:
            print(f"    {r['bin']:30s} n={int(r['n_units']):6,d}  real {f1(r['real'])} [{f1(r['real_lo'])}, {f1(r['real_hi'])}]  "
                  f"shuffled {f1(r['shuf'])}")
    print("\n== per chromosome, protein-coding genes, top0.001 ==")
    print(C.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    for col in ("GC_mean", "GC_site_weighted", "keto_mean"):
        print(f"  Spearman across chromosomes, real deficit vs {col}: {spearmanr(C['real_deficit'], C[col])[0]:+.3f} "
              f"(p = {spearmanr(C['real_deficit'], C[col])[1]:.2g}); shuffled deficit vs {col}: "
              f"{spearmanr(C['shuf_deficit'], C[col])[0]:+.3f}")
    print(f"\nCIs: paired percentile bootstrap over genes (blocks), {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out_prefix}_strata.tsv, {a.out_prefix}_chroms.tsv")


if __name__ == "__main__":
    main()
