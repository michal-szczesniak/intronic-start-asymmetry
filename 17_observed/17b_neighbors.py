#!/usr/bin/env python3
"""Step 17b: is the observed sense deficit in silent genes real initiation, or transcription
from neighbouring genes?

Step 17: in introns of protein-coding genes silent in K562 (no PRO-cap reads at their
TSSs) observed start sites are fewer on the sense strand (+28% at >=1 read), but the
deficit sits in genes where Puffin-D (100 kb input) predicts antisense-strand signal. In
active genes PRO-cap carries a large sense background proportional to host transcription,
so transcription by NEIGHBOURING genes could put background on either strand of a silent
gene. Expected strand of that background, relative to the silent gene:
  antisense  readthrough of an active opposite-strand gene (convergent, past its 3' end),
             or the upstream-antisense (PROMPT) transcription of an active same-strand gene
  sense      readthrough of an active same-strand gene upstream (tandem), or the PROMPT of
             an active opposite-strand gene (divergent)
Coverage of each strand by neighbour transcription is built from active genes (protein-
coding and lncRNA, >= --active-min TSS reads, from the step-17 table): gene body plus
--readthrough nt past its 3' end on its own strand, --prompt nt upstream of its TSS on the
other strand. Each silent gene is classed by which strands of its body are touched:
  isolated / antisense-only / sense-only / both
Settings: (readthrough, prompt) = (20 kb, 5 kb) and (50 kb, 10 kb); and a strict one where
any gene with >= 1 TSS read counts as a neighbour (50 kb, 10 kb).
Readouts: observed deficit (sites >= 1, 3, 10 reads; read sum) per class, isolated vs other
classes, the isolated genes by keto quintile, and the two models in the same classes.
If the silent-gene deficit is neighbour background, it should be large in antisense-only,
negative in sense-only and ~0 in isolated genes. Uses only the step-17 genes table and the
model tables (no bigWigs).

Usage:
  python3 17b_neighbors.py --genes results/17_obs_genes.tsv --step17 17_observed_sites.py \
      --puffin-dir ../13_genomewide/results/shuffle --procapnet-dir ../16_procapnet/results/shuffle_pc \
      --out-prefix results/17b_neighbors
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def strand_cover(act, L, rt, pr):
    """Per-strand coverage by neighbour transcription: body + readthrough on the gene's strand,
    PROMPT upstream of the TSS on the other strand."""
    d = {"+": np.zeros(L + 1, np.int32), "-": np.zeros(L + 1, np.int32)}

    def add(st, lo, hi):
        lo, hi = np.clip(lo, 0, L), np.clip(hi, 0, L)
        np.add.at(d[st], lo, 1)
        np.add.at(d[st], hi, -1)

    p, m = act[act["strand"] == "+"], act[act["strand"] == "-"]
    add("+", p["start"].to_numpy(), p["end"].to_numpy() + rt)
    add("-", m["start"].to_numpy() - rt, m["end"].to_numpy())
    add("-", p["start"].to_numpy() - pr, p["start"].to_numpy())
    add("+", m["end"].to_numpy(), m["end"].to_numpy() + pr)
    return {st: np.cumsum(v[:-1]) > 0 for st, v in d.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--step17", required=True)
    ap.add_argument("--puffin-dir", default="")
    ap.add_argument("--procapnet-dir", default="")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--active-min", type=float, default=5)
    ap.add_argument("--min-reads", default="1,3,10")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    s17 = load_module("s17", os.path.abspath(a.step17))
    rng = np.random.default_rng(a.seed)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    f1 = lambda v: f"{v:+.1f}" if np.isfinite(v) else "nan"
    okeys = [f"r{int(r)}" for r in a.min_reads.split(",")]

    G = pd.read_csv(a.genes, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
    chroms = list(dict.fromkeys(G["chrom"]))
    S = G[(G["gene_type"] == "protein_coding") & G["GC"].notna() & (G["host_tss_sense"] == 0)].copy()
    print(f"genes in the step-17 table {len(G):,}; silent protein-coding genes {len(S):,}; active neighbours "
          f"(>= {a.active_min:g} TSS reads): {int((G['host_tss_sense'] >= a.active_min).sum()):,}, "
          f"any reads: {int((G['host_tss_sense'] > 0).sum()):,}")

    MOD = []
    for pre, name, d in (("pf", "Puffin-D", a.puffin_dir), ("pc", "ProCapNet", a.procapnet_dir)):
        if d:
            n0 = len(S)
            S = S.merge(s17.load_model(d, chroms, pre), on=["chrom", "gene_id"], how="inner")
            print(f"  {name}: joined {len(S):,} of {n0:,} silent genes")
            MOD += [(f"{name} {lb}", f"{pre}_{k}") for k, lb in s17.MODEL_KEYS]
    S = S.reset_index(drop=True)
    OBS = [(f"obs {k}", f"obs_all_{k}") for k in okeys] + [("obs sum", "obs_all_sum")]
    OBS_U = [(f"obs {k} uniq", f"obs_unique_{k}") for k in okeys] + [("obs sum uniq", "obs_unique_sum")]

    settings = (("rt20k_pr5k", 20_000, 5_000, a.active_min), ("rt50k_pr10k", 50_000, 10_000, a.active_min),
                ("strict_rt50k_pr10k", 50_000, 10_000, 1e-9))
    rows = []

    def stats(d, metrics, W):
        out = {}
        for ml, col in metrics:
            Sv, Av = d[f"{col}_sense"].to_numpy(float), d[f"{col}_anti"].to_numpy(float)
            out[ml] = (float(s17.deficit(Sv.sum(), Av.sum())), s17.deficit(W @ Sv, W @ Av), Sv.sum(), Av.sum())
        return out

    def run(label, sub, lab, order, metrics, contrasts):
        res = {}
        for b in order:
            d = S[lab == b]
            if len(d) < 2:
                continue
            res[b] = stats(d, metrics, s17.boot_W(len(d), a.n_boot, rng))
            for ml, (pt, dr, sv, av) in res[b].items():
                lo, hi = s17.ci(dr)
                rows.append({"analysis": label, "subset": sub, "bin": b, "metric": ml, "n_units": len(d),
                             "sense": sv, "anti": av, "deficit": pt, "lo": lo, "hi": hi})
        for bh, bl in contrasts:
            if bh in res and bl in res:
                for ml in res[bh]:
                    v = res[bh][ml][1] - res[bl][ml][1]
                    lo, hi = s17.ci(v)
                    rows.append({"analysis": label, "subset": sub, "bin": f"{bh}-{bl}", "metric": ml,
                                 "deficit": res[bh][ml][0] - res[bl][ml][0], "lo": lo, "hi": hi, "p_boot": s17.pval(v)})

    def table(label, sub, labels, title):
        R = pd.DataFrame(rows)
        x = R[(R["analysis"] == label) & (R["subset"] == sub)]
        if not len(x):
            return
        wd = 23
        print(f"\n-- {title} --")
        print("  " + f"{'class':18s}{'n':>7s}" + "".join(f"{lb:>{wd}s}" for lb in labels))
        for b in dict.fromkeys(x["bin"]):
            xb = x[x["bin"] == b].set_index("metric")
            n = xb["n_units"].iloc[0]
            cells = [f"{f1(xb.loc[lb, 'deficit'])}[{f1(xb.loc[lb, 'lo'])},{f1(xb.loc[lb, 'hi'])}]" if lb in xb.index else "-"
                     for lb in labels]
            print("  " + f"{b:18s}" + (f"{int(n):7,d}" if np.isfinite(n) else " " * 7) + "".join(f"{c:>{wd}s}" for c in cells))
            if not np.isfinite(n):
                print("  " + " " * 25 + "".join(f"{'p=' + format(xb.loc[lb, 'p_boot'], '.2g'):>{wd}s}" if lb in xb.index
                                                else " " * wd for lb in labels))

    order = ["isolated", "antisense_only", "sense_only", "both"]
    contr = [("antisense_only", "isolated"), ("sense_only", "isolated"), ("both", "isolated")]
    for name, rt, pr, amin in settings:
        fs, fa = np.zeros(len(S)), np.zeros(len(S))
        act_all = G[G["host_tss_sense"] >= amin]
        for c in chroms:
            idx = np.flatnonzero((S["chrom"] == c).to_numpy())
            if not len(idx):
                continue
            act = act_all[act_all["chrom"] == c]
            L = int(max(G.loc[G["chrom"] == c, "end"].max(), act["end"].max() if len(act) else 0)) + rt + pr + 1
            cov = strand_cover(act, L, rt, pr)
            for i in idx:
                r = S.iloc[i]
                s, e = int(r["start"]), int(r["end"])
                sense_cov, anti_cov = (cov["+"], cov["-"]) if r["strand"] == "+" else (cov["-"], cov["+"])
                fs[i], fa[i] = sense_cov[s:e].mean(), anti_cov[s:e].mean()
        cls = np.where((fs == 0) & (fa == 0), "isolated", np.where(fs == 0, "antisense_only",
                                                                    np.where(fa == 0, "sense_only", "both")))
        S[f"class_{name}"], S[f"fsense_{name}"], S[f"fanti_{name}"] = cls, fs, fa
        print(f"\n== silent genes, neighbour setting {name} (readthrough {rt // 1000} kb, PROMPT {pr // 1000} kb, "
              f"neighbours with >= {amin:.3g} TSS reads) ==")
        print("  " + ", ".join(f"{b} {int((cls == b).sum()):,}" for b in order) +
              f"; mean fraction of the body covered: sense {fs.mean():.3f}, antisense {fa.mean():.3f}")
        run(name, "obs", cls, order, OBS, contr)
        run(name, "obs_unique", cls, order, OBS_U, contr)
        table(name, "obs", [m for m, _ in OBS], "observed, all positions")
        table(name, "obs_unique", [m for m, _ in OBS_U], "observed, unique positions")
        if MOD:
            run(name, "models", cls, order, MOD, contr)
            table(name, "models", [m for m, _ in MOD], "models")

    # isolated silent genes (strict setting) by keto quintile
    iso = S["class_strict_rt50k_pr10k"].to_numpy() == "isolated"
    e = np.quantile(S["keto"].to_numpy(float), [0.2, 0.4, 0.6, 0.8])
    q = np.where(iso, np.array(["Q1", "Q2", "Q3", "Q4", "Q5"])[np.searchsorted(e, S["keto"].to_numpy(float), side="right")], "-")
    print(f"\n== strictly isolated silent genes ({int(iso.sum()):,}) by keto quintile (edges from all silent genes) ==")
    run("iso_keto", "obs", q, ["Q1", "Q2", "Q3", "Q4", "Q5"], OBS, [("Q5", "Q1")])
    table("iso_keto", "obs", [m for m, _ in OBS], "observed")
    if MOD:
        run("iso_keto", "models", q, ["Q1", "Q2", "Q3", "Q4", "Q5"], MOD, [("Q5", "Q1")])
        table("iso_keto", "models", [m for m, _ in MOD], "models")

    R = pd.DataFrame(rows)
    R.to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)
    S.drop(columns=[c for c in S.columns if c.startswith(("pf_", "pc_"))]).to_csv(a.out_prefix + "_silent_genes.tsv",
                                                                                   sep="\t", index=False)
    print(f"\nCIs: percentile bootstrap over genes, {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {a.out_prefix}_strata.tsv, {a.out_prefix}_silent_genes.tsv")


if __name__ == "__main__":
    main()
