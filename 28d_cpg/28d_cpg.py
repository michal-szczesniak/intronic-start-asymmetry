#!/usr/bin/env python3
"""Step 28d (EXPLORATORY; specified on 10 Nov 2026 after step 28c, before running): are CpG islands behind the
strand difference in the first 5 kb of genes?

Background. Steps 28b and 28c placed the deficit of strong sense sites (Results 4) and most of the weaker sense start
of pairs (Results 1) within 1-5 kb of the gene's 5' end. Promoter CpG islands often reach into this region, and
antisense initiation is known to occur near their downstream edge; such initiation would add antisense sites there
without any suppression of sense starts.

Question (specified before running). After removing sites (or pairs) within 2 kb of any CpG island (UCSC
cpgIslandExt, hg38), (1) is the deficit of strong sites within 1-5 kb of the gene's 5' end still above zero, and
(2) is the bias of pairs within 1-5 kb still above zero? Reading rule: if the 95% CI stays above zero, CpG islands
do not explain the concentration near the gene's start; if the estimate falls to about the level beyond 5 kb, the
concentration is tied to CpG islands. Margins of 1 and 5 kb, sites and pairs within 2 kb of an island, and the
difference < 5 kb minus >= 5 kb are reported as well. Two-sided bootstrap P, no correction; all exploratory.

Inputs: step-28b tables (sites with positions, pairs with midpoints), step-28 gene tables (active genes, activity
thirds; all of them enter the resampling), the CpG-island table, step-17 code (bootstrap helpers).
Distance of a site or pair midpoint to the nearest island: 0 inside an island, otherwise to the nearer island end.
Checks: with no exclusion, the < 5 kb estimates must equal step 28b (STOP otherwise); the island table must hold
>= 20,000 islands on the analysed chromosomes (STOP otherwise).

Usage: python3 28d_cpg.py --cgi cpgIslandExt.txt.gz --r28b <28b results> --r28 <28 results> \
           --step17 17_observed_sites.py --lines K562,HCT116 --out-prefix results/28d
Outputs: <prefix>_strata.tsv (every estimate), <prefix>_sites.tsv.gz and <prefix>_pairs.tsv (with distance to the
nearest island); log printed.
"""
import argparse
import gzip
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

MARGINS = (1, 2, 5)          # kb; 2 kb is the margin of the question
N_BOOT = 2000


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def wmedian(v, w):
    o = np.argsort(v, kind="stable")
    v, c = v[o], np.cumsum(w[o])
    N = int(c[-1]) if len(c) else 0
    if N == 0:
        return np.nan
    if N % 2:
        return float(v[np.searchsorted(c, N // 2, side="right")])
    return 0.5 * float(v[np.searchsorted(c, N // 2 - 1, side="right")] + v[np.searchsorted(c, N // 2, side="right")])


def read_cgi(path):
    """UCSC cpgIslandExt: bin, chrom, chromStart, chromEnd, name, length, cpgNum, gcNum, perCpg, perGc, obsExp."""
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 11:
                sys.exit(f"STOP: {path}: expected 11 columns (UCSC cpgIslandExt), found {len(f)}")
            rows.append((f[1], int(f[2]), int(f[3])))
    d = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    return {c: (g["start"].to_numpy(), g["end"].to_numpy()) for c, g in d.sort_values(["chrom", "start"]).groupby("chrom")}


def dist_to_cgi(cgi, chrom, pos):
    """Distance (nt) from each 0-based position to the nearest island [start, end) on the same chromosome: 0 inside,
    otherwise to the nearest island base; inf if the chromosome has no island."""
    out = np.full(len(pos), np.inf)
    for c in np.unique(chrom):
        if c not in cgi:
            continue
        m = chrom == c
        s, e = cgi[c]
        emax = np.maximum.accumulate(e)              # furthest end among islands starting at or before each one
        p = pos[m].astype(np.int64)
        k = np.searchsorted(s, p, side="right") - 1  # last island starting at or before p (-1 if none)
        kk = np.clip(k, 0, len(s) - 1)
        left = np.where(k >= 0, np.where(p < emax[kk], 0, p - (emax[kk] - 1)), np.inf)
        nxt = np.clip(k + 1, 0, len(s) - 1)
        right = np.where(k + 1 < len(s), s[nxt] - p, np.inf)
        out[m] = np.minimum(left, right)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cgi", required=True)
    ap.add_argument("--r28b", required=True)
    ap.add_argument("--r28", required=True)
    ap.add_argument("--step17", required=True)
    ap.add_argument("--lines", default="K562,HCT116")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-islands", type=int, default=20000, help="STOP below this number of islands (tests: 0)")
    a = ap.parse_args()
    s17 = load_module("s17", a.step17)
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    cgi = read_cgi(a.cgi)
    rows, all_sites, all_pairs = [], [], []
    print(f"== 28d (EXPLORATORY): CpG islands and the strand difference near the gene's start; bootstrap {N_BOOT}")
    for line in [x for x in a.lines.split(",") if x]:
        rng = np.random.default_rng(a.seed)
        S28b = pd.read_csv(f"{a.r28b}/28b_{line}_strata.tsv", sep="\t")
        G = pd.read_csv(f"{a.r28}/28_{line}_genes.tsv", sep="\t", dtype={"gene_id": str})
        chroms = sorted(G["chrom"].unique())
        n_isl = sum(len(cgi[c][0]) for c in chroms if c in cgi)
        print(f"\n### {line}: CpG islands on the {len(chroms)} chromosomes of the gene table: {n_isl:,}")
        if n_isl < a.min_islands:
            sys.exit(f"STOP: fewer than {a.min_islands:,} CpG islands on the analysed chromosomes")

        # ---------------- sites
        G = G[G["activity"].isin(["A1", "A2", "A3"])].reset_index(drop=True)
        gi = {g: i for i, g in enumerate(G["gene_id"])}
        parts = [np.flatnonzero(G["activity"].to_numpy() == t) for t in ("A1", "A2", "A3")]
        Ws = [s17.boot_W(len(ix), N_BOOT, rng) for ix in parts]
        S = pd.read_csv(f"{a.r28b}/28b_{line}_sites.tsv.gz", sep="\t", dtype={"gene_id": str})
        S = S[(S["metric"] == "strong") & S["gene_id"].isin(gi)].reset_index(drop=True)
        S["cgi_dist_kb"] = dist_to_cgi(cgi, S["chrom"].to_numpy(), S["pos"].to_numpy()) / 1e3
        all_sites.append(S.assign(line=line))
        g_of = S["gene_id"].map(gi).to_numpy()
        sense = (S["strand"] == "sense").to_numpy()
        near5 = (S["dist5_kb"] < 5).to_numpy()
        cd = S["cgi_dist_kb"].to_numpy()

        def site_est(m):
            cs = np.bincount(g_of[m & sense], minlength=len(G)).astype(float)
            ca = np.bincount(g_of[m & ~sense], minlength=len(G)).astype(float)
            bS = sum(W @ cs[ix] for W, ix in zip(Ws, parts))
            bA = sum(W @ ca[ix] for W, ix in zip(Ws, parts))
            return float(s17.deficit(cs.sum(), ca.sum())), s17.deficit(bS, bA), int(cs.sum()), int(ca.sum())

        print(f"  strong sites in {len(G):,} active protein-coding genes: within 2 kb of an island, < 5 kb from the 5' end: "
              f"sense {int((near5 & sense & (cd <= 2)).sum()):,} of {int((near5 & sense).sum()):,}, antisense "
              f"{int((near5 & ~sense & (cd <= 2)).sum()):,} of {int((near5 & ~sense).sum()):,}; >= 5 kb: sense "
              f"{int((~near5 & sense & (cd <= 2)).sum()):,} of {int((~near5 & sense).sum()):,}, antisense "
              f"{int((~near5 & ~sense & (cd <= 2)).sum()):,} of {int((~near5 & ~sense).sum()):,}")
        print("  deficit of strong sites, %, [95% CI], sense / antisense sites:")
        subsets = [("all sites", np.ones(len(S), bool))] + [(f"> {m} kb from an island", cd > m) for m in MARGINS] + \
                  [("<= 2 kb from an island", cd <= 2)]
        for nm, keep in subsets:
            r = {}
            for reg, rm in (("<5 kb", near5), (">=5 kb", ~near5)):
                pt, dr, ns, na = site_est(keep & rm)
                lo, hi = s17.ci(dr)
                r[reg] = (pt, dr)
                rows.append({"line": line, "analysis": "sites_strong", "subset": nm, "region": reg, "n": ns + na,
                             "sense": ns, "anti": na, "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan})
                print(f"    {nm:24s} {reg:7s} {pt:+6.1f} [{lo:+.1f}, {hi:+.1f}]  {ns:,} / {na:,}")
            v = r["<5 kb"][1] - r[">=5 kb"][1]
            lo, hi = s17.ci(v)
            p = s17.pval(v)
            rows.append({"line": line, "analysis": "sites_strong", "subset": nm, "region": "<5 minus >=5 kb", "n": int(keep.sum()),
                         "sense": np.nan, "anti": np.nan, "estimate": r["<5 kb"][0] - r[">=5 kb"][0], "lo": lo, "hi": hi, "p_boot": p})
            print(f"    {nm:24s} <5 minus >=5 kb {r['<5 kb'][0] - r['>=5 kb'][0]:+6.1f} [{lo:+.1f}, {hi:+.1f}], p = {p:.2g}")
            if nm == "all sites":
                ref = S28b[(S28b["analysis"] == "sites_strong_bin5") & (S28b["stratum"] == "1-5")]["estimate"].iloc[0]
                if not np.isclose(r["<5 kb"][0], ref):
                    sys.exit(f"STOP: {line}: sites < 5 kb {r['<5 kb'][0]} differ from step 28b {ref}")

        # ---------------- pairs
        P = pd.read_csv(f"{a.r28b}/28b_{line}_pairs.tsv", sep="\t", dtype={"host": str})
        P["cgi_dist_kb"] = dist_to_cgi(cgi, P["chrom"].to_numpy(), P["mid"].to_numpy()) / 1e3
        all_pairs.append(P.assign(line=line))
        hosts, hidx = np.unique(P["host"].to_numpy(), return_inverse=True)
        Wp = np.stack([np.bincount(rng.integers(0, len(hosts), len(hosts)), minlength=len(hosts)) for _ in range(N_BOOT)])
        DI, C = P["DI"].to_numpy(float), P["ctrl_DI"].to_numpy(float)
        pn5 = (P["dist5_kb"] < 5).to_numpy()
        pcd = P["cgi_dist_kb"].to_numpy()

        def pair_est(m):
            pt = 100 * (1 - 2 ** (np.median(DI[m]) - np.median(C[m])))
            bt = np.array([100 * (1 - 2 ** (wmedian(DI[m], w[hidx[m]]) - wmedian(C[m], w[hidx[m]]))) for w in Wp])
            return pt, bt

        print(f"  pairs: {len(P):,}; within 2 kb of an island: {int((pn5 & (pcd <= 2)).sum()):,} of {int(pn5.sum()):,} "
              f"(< 5 kb from the 5' end) and {int((~pn5 & (pcd <= 2)).sum()):,} of {int((~pn5).sum()):,} (>= 5 kb)")
        print("  bias of pairs, %, [95% CI]:")
        subsets = [("all pairs", np.ones(len(P), bool))] + [(f"> {m} kb from an island", pcd > m) for m in MARGINS] + \
                  [("<= 2 kb from an island", pcd <= 2)]
        for nm, keep in subsets:
            r = {}
            for reg, rm in (("<5 kb", pn5), (">=5 kb", ~pn5)):
                m = keep & rm
                if m.sum() < 20:
                    print(f"    {nm:24s} {reg:7s} n = {int(m.sum())}: too few pairs")
                    continue
                pt, bt = pair_est(m)
                lo, hi = s17.ci(bt)
                r[reg] = (pt, bt)
                rows.append({"line": line, "analysis": "pairs", "subset": nm, "region": reg, "n": int(m.sum()),
                             "sense": np.nan, "anti": np.nan, "estimate": pt, "lo": lo, "hi": hi, "p_boot": np.nan})
                print(f"    {nm:24s} {reg:7s} {pt:+6.1f} [{lo:+.1f}, {hi:+.1f}]  n = {int(m.sum()):,}")
            if len(r) == 2:
                v = r["<5 kb"][1] - r[">=5 kb"][1]
                lo, hi = s17.ci(v)
                p = s17.pval(v)
                rows.append({"line": line, "analysis": "pairs", "subset": nm, "region": "<5 minus >=5 kb", "n": int(keep.sum()),
                             "sense": np.nan, "anti": np.nan, "estimate": r["<5 kb"][0] - r[">=5 kb"][0], "lo": lo, "hi": hi,
                             "p_boot": p})
                print(f"    {nm:24s} <5 minus >=5 kb {r['<5 kb'][0] - r['>=5 kb'][0]:+6.1f} [{lo:+.1f}, {hi:+.1f}], p = {p:.2g}")
            if nm == "all pairs":
                ref = S28b[(S28b["analysis"] == "pairs_dist_bin") & (S28b["stratum"] == "1-5")]["estimate"].iloc[0]
                if not np.isclose(r["<5 kb"][0], ref):
                    sys.exit(f"STOP: {line}: pairs < 5 kb {r['<5 kb'][0]} differ from step 28b {ref}")

        for an, lab in (("sites_strong", "deficit of strong sites"), ("pairs", "bias of pairs")):
            x = [rw for rw in rows if rw["line"] == line and rw["analysis"] == an and rw["subset"] == "> 2 kb from an island"
                 and rw["region"] == "<5 kb"]
            if x:
                rw = x[0]
                print(f"  ANSWER {line}, {lab} within 1-5 kb, > 2 kb from any island: {rw['estimate']:+.1f} "
                      f"[{rw['lo']:+.1f}, {rw['hi']:+.1f}] (n = {rw['n']:,})")

    pd.DataFrame(rows).to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False, float_format="%.6g")
    pd.concat(all_sites).to_csv(a.out_prefix + "_sites.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.concat(all_pairs).to_csv(a.out_prefix + "_pairs.tsv", sep="\t", index=False)
    print(f"\nwritten: {a.out_prefix}_strata.tsv, _sites.tsv.gz, _pairs.tsv")


if __name__ == "__main__":
    main()
