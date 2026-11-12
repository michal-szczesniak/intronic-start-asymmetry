#!/usr/bin/env python3
"""Step 13c: genome-wide H1 map and shuffle test, from the per-chromosome step-10 tables of 13b.

Each 13b task wrote, for one chromosome, results/shuffle/shuffle_<chrom>_{genes,blocks}.tsv:
per gene (eligible intronic positions of protein-coding and lncRNA genes) and per 1 Mb
intergenic block, the number of potential start sites on the sense and antisense
strand (plus / minus for intergenic) at three strengths (top 1e-3 / 3e-4 / 1e-4 of the
chr8+chr9 predictions, fixed thresholds 0.0439 / 0.1432 / 0.3293), for the real
sequence and for the composition-preserving (dinucleotide) shuffle, plus the sum of
the predictions. This script pools them.

For every group of chromosomes and every set (protein-coding introns, lncRNA introns,
intergenic) it reports, as step 10:
  real deficit  = 100 * (1 - total sense / total antisense)   (> 0: fewer sites on the sense strand)
  shuffled deficit, and drop = real - shuffled (negative: the deficit GROWS after
  shuffling, i.e. composition alone predicts more than is seen)
each with a paired percentile bootstrap over units (genes / 1 Mb blocks).

Groups: genome (all chromosomes given), each chromosome, and
  heldout_test  = chr8, chr9   (Puffin-D test chromosomes)
  heldout_valid = chr10        (validation)
  heldout       = chr8, chr9, chr10
  training      = all others.
Puffin-D was trained on the observed PRO-cap of the training chromosomes. Its
predictions for their REAL sequence may therefore partly reproduce the observed
signal, including its non-sequence (chromatin) part, while shuffled sequences are new
for every chromosome. The real deficit and the drop can thus be inflated on training
chromosomes. The held-out group is the clean estimate; the training - heldout
difference (independent bootstraps of the two groups) shows whether this matters.

Consistency across chromosomes: for each strength, the number of chromosomes with a
real deficit > 0 and with a drop < 0, with a two-sided sign test.

Optional check (--check-genes): the real-sequence counts of chr8/chr9 must equal those
of the original step-10 run (same tracks, thresholds and masks), gene by gene. The
shuffled counts are not expected to match exactly (13b shuffles each chromosome with
its own random stream, and GPU predictions can differ in the last digits).

Usage:
  python3 13c_shuffle_combine.py --dir results/shuffle --chroms chr1,...,chr22 \
      --out-prefix results/13c_shuffle_gw [--n-boot 2000] [--seed 1] [--skip-missing] \
      [--check-genes ../10_shuffle/results/shuffle_chr8_chr9_genes.tsv]
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

TKEYS = ("top0.001", "top0.0003", "top0.0001")
METRICS = TKEYS + ("sum",)
SEQS = ("unique", "all")
HELDOUT = {"heldout_test": ("chr8", "chr9"), "heldout_valid": ("chr10",), "heldout": ("chr8", "chr9", "chr10")}


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
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (np.nan, np.nan)


def p_two_sided(v):
    v = v[np.isfinite(v)]
    return max(2 * min((v <= 0).mean(), (v >= 0).mean()), 1 / len(v)) if len(v) else np.nan


def sign_test(k, n):
    """Two-sided exact binomial test of k successes in n with p = 0.5."""
    if n == 0:
        return np.nan
    pk = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
    return min(1.0, sum(p for p in pk if p <= pk[k] * (1 + 1e-9)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory with shuffle_<chrom>_{genes,blocks}.tsv")
    ap.add_argument("--chroms", required=True, help="comma-separated")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--skip-missing", action="store_true", help="pool the chromosomes that are ready, list the others")
    ap.add_argument("--check-genes", default="", help="step-10 genes table of chr8+chr9 to compare real counts with")
    a = ap.parse_args()
    chroms = [c for c in a.chroms.split(",") if c]
    if len(set(chroms)) != len(chroms):
        sys.exit("STOP: a chromosome is listed twice")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)

    G, B, missing = [], [], []
    for c in chroms:
        fg = os.path.join(a.dir, f"shuffle_{c}_genes.tsv")
        fb = os.path.join(a.dir, f"shuffle_{c}_blocks.tsv")
        fs = os.path.join(a.dir, f"shuffle_{c}_summary.tsv")      # written last by step 10: the done marker
        if not (os.path.isfile(fg) and os.path.isfile(fb) and os.path.isfile(fs)):
            missing.append(c)
            continue
        g = pd.read_csv(fg, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
        b = pd.read_csv(fb, sep="\t")
        for df in (g, b):
            if not (df["chrom"].astype(str) == c).all():
                sys.exit(f"STOP: {c} tables contain other chromosomes")
        G.append(g)
        B.append(b)
    if missing and not a.skip_missing:
        sys.exit(f"STOP: not ready: {','.join(missing)} (use --skip-missing to pool the others)")
    chroms = [c for c in chroms if c not in missing]
    if not chroms:
        sys.exit("STOP: no chromosome ready")
    G = pd.concat(G, ignore_index=True)
    B = pd.concat(B, ignore_index=True)
    seeds = sorted({int(col.split("_")[0][4:]) for col in G.columns if col.startswith("shuf") and col.endswith("_sense")})
    for df, nm in ((G, "genes"), (B, "blocks")):
        need = [f"{t}_{s}_{m}_{d}" for t in ["real"] + [f"shuf{x}" for x in seeds] for s in SEQS for m in METRICS
                for d in ("sense", "anti")]
        lack = [x for x in need if x not in df.columns]
        if lack:
            sys.exit(f"STOP: {nm} tables lack columns, e.g. {lack[:3]}")
        nan_rows = int(df[need].isna().any(axis=1).sum())
        if nan_rows:
            sys.exit(f"STOP: {nan_rows} {nm} rows have empty values (a chromosome table is incomplete?)")
    if G.duplicated(["chrom", "gene_id"]).any():
        sys.exit("STOP: duplicated genes")
    print(f"pooled {len(chroms)} chromosomes: {','.join(chroms)}")
    if missing:
        print(f"NOT READY (left out): {','.join(missing)}")
    cnt = G.groupby("chrom")["gene_type"].value_counts().unstack(fill_value=0)
    print(f"genes: {len(G):,} ({int((G['gene_type'] == 'protein_coding').sum()):,} protein-coding, "
          f"{int((G['gene_type'] == 'lncRNA').sum()):,} lncRNA); intergenic 1 Mb blocks: {len(B):,}; shuffle seeds: {seeds}")
    print("per chromosome (protein-coding / lncRNA / blocks): " +
          "; ".join(f"{c} {int(cnt.loc[c].get('protein_coding', 0))}/{int(cnt.loc[c].get('lncRNA', 0))}/"
                    f"{int((B['chrom'] == c).sum())}" for c in chroms))

    if a.check_genes:
        ref = pd.read_csv(a.check_genes, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str})
        both = [c for c in ref["chrom"].astype(str).unique() if c in chroms]
        rc = [x for x in ref.columns if x.startswith("real_")]
        A_ = G[G["chrom"].isin(both)].merge(ref[ref["chrom"].isin(both)], on=["chrom", "gene_id"], suffixes=("", "_ref"),
                                            how="outer", indicator=True)
        only = int((A_["_merge"] != "both").sum())
        A_ = A_[A_["_merge"] == "both"]
        diff = sum(int((A_[x] != A_[x + "_ref"]).sum()) for x in rc if not x.endswith(("sum_sense", "sum_anti")))
        sdiff = max(float(np.abs(A_[x] - A_[x + "_ref"]).max()) for x in rc if x.endswith(("sum_sense", "sum_anti")))
        print(f"\ncheck vs {a.check_genes} ({','.join(both)}): {len(A_):,} genes in both, {only} in only one; "
              f"real site counts differing: {diff} cells (expected 0); largest difference in real sums: {sdiff:.3g}")
        sc = [x for x in ref.columns if x.startswith("shuf1_") and not x.endswith(("sum_sense", "sum_anti"))]
        if sc and "shuf1_all_top0.001_sense" in G.columns:
            for c in both:
                Ac = A_[A_["chrom"] == c]
                same = np.mean([(Ac[x] == Ac[x + "_ref"]).mean() for x in sc])
                print(f"  {c}: shuffled (seed 1) site counts equal in {100 * same:.1f}% of cells "
                      f"(13b uses its own random stream per chromosome; equality not required)")

    groups = [("genome", tuple(chroms))]
    for gname, gch in HELDOUT.items():
        sel = tuple(c for c in gch if c in chroms)
        if sel:
            groups.append((gname, sel))
    train = tuple(c for c in chroms if c not in HELDOUT["heldout"])
    if train:
        groups.append(("training", train))
    groups += [(c, (c,)) for c in chroms]

    sets = [("intron_protein_coding", G[G["gene_type"] == "protein_coding"], "gene"),
            ("intron_lncRNA", G[G["gene_type"] == "lncRNA"], "gene"),
            ("intergenic_plus_vs_minus", B, "1Mb_block")]
    rng = np.random.default_rng(a.seed)
    rows, boots = [], {}
    for gname, gch in groups:
        for label, df, unit in sets:
            d = df[df["chrom"].isin(gch)]
            if not len(d):
                continue
            W = boot_W(len(d), a.n_boot, rng)
            for seq in SEQS:
                for m in METRICS:
                    So, Ao = d[f"real_{seq}_{m}_sense"].to_numpy(float), d[f"real_{seq}_{m}_anti"].to_numpy(float)
                    for sd in seeds:
                        Ss, As = d[f"shuf{sd}_{seq}_{m}_sense"].to_numpy(float), d[f"shuf{sd}_{seq}_{m}_anti"].to_numpy(float)
                        if min(Ao.sum(), As.sum()) == 0:
                            continue
                        dr, dsh = float(deficit(So.sum(), Ao.sum())), float(deficit(Ss.sum(), As.sum()))
                        br, bs = deficit(W @ So, W @ Ao), deficit(W @ Ss, W @ As)
                        boots[(gname, label, seq, m, sd)] = (br, bs)
                        (rlo, rhi), (slo, shi), (dlo, dhi) = ci(br), ci(bs), ci(br - bs)
                        rows.append({"group": gname, "chroms": ",".join(gch), "set": label, "unit": unit,
                                     "n_units": len(d), "seq": seq, "metric": m, "seed": sd,
                                     "real_sense": So.sum(), "real_anti": Ao.sum(),
                                     "real_deficit_pct": dr, "real_lo95": rlo, "real_hi95": rhi,
                                     "shuf_sense": Ss.sum(), "shuf_anti": As.sum(),
                                     "shuf_deficit_pct": dsh, "shuf_lo95": slo, "shuf_hi95": shi,
                                     "drop_pct_points": dr - dsh, "drop_lo95": dlo, "drop_hi95": dhi,
                                     "p_drop": p_two_sided(br - bs),
                                     "remaining_fraction": dsh / dr if dr != 0 else np.nan,
                                     "sites_shuf_vs_real": (Ss.sum() + As.sum()) / (So.sum() + Ao.sum())})
    R = pd.DataFrame(rows)
    R.to_csv(a.out_prefix + "_summary.tsv", sep="\t", index=False)

    # training - heldout, independent bootstraps of the two groups
    drows = []
    if "training" in dict(groups) and "heldout" in dict(groups):
        for (gname, label, seq, m, sd), (br_t, bs_t) in boots.items():
            if gname != "training" or (("heldout", label, seq, m, sd) not in boots):
                continue
            br_h, bs_h = boots[("heldout", label, seq, m, sd)]
            rt = R[(R.group == "training") & (R.set == label) & (R.seq == seq) & (R.metric == m) & (R.seed == sd)].iloc[0]
            rh = R[(R.group == "heldout") & (R.set == label) & (R.seq == seq) & (R.metric == m) & (R.seed == sd)].iloc[0]
            for q, vt, vh, bt, bh in (("real_deficit", rt.real_deficit_pct, rh.real_deficit_pct, br_t, br_h),
                                      ("shuf_deficit", rt.shuf_deficit_pct, rh.shuf_deficit_pct, bs_t, bs_h),
                                      ("drop", rt.drop_pct_points, rh.drop_pct_points, br_t - bs_t, br_h - bs_h)):
                lo, hi = ci(bt - bh)
                drows.append({"set": label, "seq": seq, "metric": m, "seed": sd, "quantity": q, "training": vt,
                              "heldout": vh, "training_minus_heldout": vt - vh, "lo95": lo, "hi95": hi,
                              "p_boot": p_two_sided(bt - bh)})
    Dg = pd.DataFrame(drows)
    Dg.to_csv(a.out_prefix + "_training_vs_heldout.tsv", sep="\t", index=False)

    # consistency across chromosomes
    crow = []
    per = R[R["group"].isin(chroms)]
    for (label, seq, m, sd), x in per.groupby(["set", "seq", "metric", "seed"], sort=False):
        n = len(x)
        kr = int((x["real_deficit_pct"] > 0).sum())
        kd = int((x["drop_pct_points"] < 0).sum())
        crow.append({"set": label, "seq": seq, "metric": m, "seed": sd, "n_chroms": n,
                     "chroms_real_deficit_gt0": kr, "p_sign_real": sign_test(kr, n),
                     "chroms_drop_lt0": kd, "p_sign_drop": sign_test(kd, n),
                     "real_deficit_min": x["real_deficit_pct"].min(), "real_deficit_max": x["real_deficit_pct"].max(),
                     "shuf_deficit_min": x["shuf_deficit_pct"].min(), "shuf_deficit_max": x["shuf_deficit_pct"].max()})
    Cn = pd.DataFrame(crow)
    Cn.to_csv(a.out_prefix + "_consistency.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    f1 = lambda v: f"{v:+.1f}"
    show = lambda r: (f"real {f1(r.real_deficit_pct)} [{f1(r.real_lo95)}, {f1(r.real_hi95)}]   "
                      f"shuffled {f1(r.shuf_deficit_pct)} [{f1(r.shuf_lo95)}, {f1(r.shuf_hi95)}]   "
                      f"drop {f1(r.drop_pct_points)} [{f1(r.drop_lo95)}, {f1(r.drop_hi95)}]")
    sd0 = seeds[0]
    big = ["genome", "heldout_test", "heldout_valid", "heldout", "training"]
    for label, title in (("intron_protein_coding", "protein-coding introns: sense deficit (%) of potential start sites"),
                         ("intergenic_plus_vs_minus", "intergenic control: plus vs minus (should be ~0)"),
                         ("intron_lncRNA", "lncRNA introns")):
        print(f"\n== {title}; seq = all, shuffle seed {sd0} ==")
        for m in (TKEYS if label == "intron_protein_coding" else TKEYS[:1]):
            print(f"  {m}")
            for gname in big:
                x = R[(R.group == gname) & (R.set == label) & (R.seq == "all") & (R.metric == m) & (R.seed == sd0)]
                if len(x):
                    r = x.iloc[0]
                    print(f"    {gname:14s} n={r.n_units:6,d}  sites {int(r.real_sense + r.real_anti):8,d}  {show(r)}")
    print(f"\n== protein-coding introns, per chromosome (seq = all, top0.001) ==")
    for c in chroms:
        x = R[(R.group == c) & (R.set == "intron_protein_coding") & (R.seq == "all") & (R.metric == "top0.001") & (R.seed == sd0)]
        if len(x):
            r = x.iloc[0]
            print(f"    {c:6s} n={r.n_units:5,d}  sites {int(r.real_sense + r.real_anti):7,d}  {show(r)}")
    print("\n== consistency across chromosomes (seq = all) ==")
    if len(Cn):
        print(Cn[(Cn.seq == "all") & (Cn.seed == sd0)].drop(columns=["seq", "seed"]).to_string(
            index=False, float_format=lambda v: f"{v:.3g}"))
    if len(Dg):
        print("\n== training - heldout (seq = all); real deficit / drop larger on training chromosomes would point "
              "to memorised signal ==")
        print(Dg[(Dg.seq == "all") & (Dg.seed == sd0) & (Dg.set != "intron_lncRNA")].drop(columns=["seq", "seed"]).to_string(
            index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n== background: sum of predictions, protein-coding introns (seq = all) ==")
    for gname in big:
        x = R[(R.group == gname) & (R.set == "intron_protein_coding") & (R.seq == "all") & (R.metric == "sum") & (R.seed == sd0)]
        if len(x):
            r = x.iloc[0]
            print(f"    {gname:14s} {show(r)}")
    print(f"\nCIs: paired percentile bootstrap over units, {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out_prefix}_summary.tsv, _training_vs_heldout.tsv, _consistency.tsv")


if __name__ == "__main__":
    main()
