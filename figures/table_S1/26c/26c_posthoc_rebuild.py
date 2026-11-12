#!/usr/bin/env python3
"""Step 26c (rebuild, 11 Nov 2026): post hoc weighting of the DNA-methylation comparison of step 26b.

Why a rebuild. The original step-26c script (26c_weights_posthoc.py, md5 8340fc22abf81c1cf21919ef1a648e63)
and its result tables were not kept; only its primary log
survived (26c_posthoc_primary.log as sent on 9 Nov 2026, md5 3a082807aeaa6052660f5be40d629757). This script
recomputes the same quantities from the step-26b gene tables and writes a log in the same layout; the run is
accepted only if every number of that surviving log is reproduced exactly (--check-log).

Question (chosen after seeing step 26b, therefore post hoc). The planned comparison of step 26b pools the
strong sites of the three activity thirds, so each activity third counts in proportion to its number of
sites. Within the highly methylated third most opposite-strand sites come from the least active genes, which
show no deficit. Does the high-minus-low methylation contrast change when each activity third gets the same
weight?

Definitions (as in step 26b)
  genes      active protein-coding genes with >= 20 CpGs (activity thirds A1-A3 and methylation thirds M1-M3
             formed within each activity third, columns 'activity' and 'meth_third' of the 26b gene table)
  deficit    100 x (1 - sense / antisense), summed over the genes of a group
  pooled     contrast of the deficits of the methylation thirds pooled over A1-A3 (the step-26b comparison)
  equal      mean over A1-A3 of the within-third contrasts (each activity third the same weight)
  bootstrap  2,000 resamples, seed 1; genes resampled with replacement within each of the nine cells
             (A x M), drawn once per cell in the order A1-M1, A1-M2, A1-M3, A2-M1, ... and shared by all metrics;
             95% percentile intervals; two-sided p = max(1/B, 2 x min(share <= 0, share >= 0))
  metrics    measured sites with >= 10, >= 3 and >= 1 reads (obs_all_r10/r3/r1), strong potential sites of
             Puffin-D (pf_top0.0001) and of ProCapNet (pc_top0.0001)
  cells      per cell: genes, median host activity (host_tss_sense), median eligible kb, median methylation,
             strong sites (>= 10 reads) on each strand, opposite-strand sites per 100 kb of eligible sequence,
             deficit, and weight = share of the opposite-strand strong sites of the methylation third that the
             cell contributes (the weight it has in the pooled comparison)
  rho        Spearman correlation, inside each activity third, of methylation with host activity and with
             eligible nt

Usage: python3 26c_posthoc_rebuild.py --genes 26b_meth_<set>_genes.tsv --set primary|alternative
           --out-prefix results/26c_posthoc [--expect-md5 <md5 of the gene table>] [--check-log <old log>]
Writes <prefix>_<set>_contrasts.tsv, <prefix>_<set>_cells.tsv, <prefix>_<set>_rho.tsv and <prefix>_<set>.log.
"""
import argparse
import hashlib
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ACTS = ["A1", "A2", "A3"]
MS = ["M1", "M2", "M3"]
METRICS = [("obs r10", "obs_all_r10"), ("obs r3", "obs_all_r3"), ("obs r1", "obs_all_r1"),
           ("Puffin-D strong", "pf_top0.0001"), ("ProCapNet strong", "pc_top0.0001")]
B, SEED = 2000, 1


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def deficit(s, a):
    return 100.0 * (1.0 - s / a)


def pval(v):
    return max(1.0 / B, 2.0 * min(float(np.mean(v <= 0)), float(np.mean(v >= 0))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--set", required=True, choices=["primary", "alternative"])
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--expect-md5")
    ap.add_argument("--check-log")
    a = ap.parse_args()
    m5 = md5(a.genes)
    if a.expect_md5 and m5 != a.expect_md5:
        sys.exit(f"STOP: md5 of {a.genes} is {m5}, expected {a.expect_md5}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)) or ".", exist_ok=True)
    G = pd.read_csv(a.genes, sep="\t")
    A = G[G["activity"].isin(ACTS)].reset_index(drop=True)
    if set(A["meth_third"]) != set(MS):
        sys.exit("STOP: methylation thirds missing in the active genes")
    log = []
    say = log.append
    say(f"== 26c post hoc: {os.path.basename(a.genes)}; active protein-coding genes {len(A):,}; bootstrap {B}, seed {SEED}")

    rng = np.random.default_rng(SEED)
    cells = {}
    for ac in ACTS:
        for m in MS:
            ix = np.flatnonzero(((A["activity"] == ac) & (A["meth_third"] == m)).to_numpy())
            cells[(ac, m)] = (ix, rng.integers(0, len(ix), size=(B, len(ix))))

    # 1. cells, strong sites
    S10, N10 = A["obs_all_r10_sense"].to_numpy(float), A["obs_all_r10_anti"].to_numpy(float)
    crow = []
    say("")
    say("== 1. cells, strong sites (>= 10 reads) ==")
    say("  cell     genes    host  elig kb   meth   sense    anti anti/100kb  deficit  weight")
    for m in MS:
        tot_anti = sum(N10[cells[(ac, m)][0]].sum() for ac in ACTS)
        for ac in ACTS:
            ix = cells[(ac, m)][0]
            d = A.iloc[ix]
            s, n = S10[ix].sum(), N10[ix].sum()
            r = {"set": a.set, "cell": f"{ac}-{m}", "activity_third": ac, "meth_third": m, "genes": len(ix),
                 "host_activity_median": float(d["host_tss_sense"].median()),
                 "eligible_kb_median": float(d["eligible_nt"].median()) / 1000.0, "meth_median": float(d["meth"].median()),
                 "sense_r10": s, "anti_r10": n, "anti_per_100kb": n / d["eligible_nt"].sum() * 1e5,
                 "deficit_r10": deficit(s, n), "weight_in_pooled": n / tot_anti}
            crow.append(r)
            say(f"  {r['cell']:7s} {r['genes']:6,} {r['host_activity_median']:7.0f} {r['eligible_kb_median']:8.1f} "
                f"{r['meth_median']:6.3f} {s:7.0f} {n:7.0f} {r['anti_per_100kb']:10.3f} {r['deficit_r10']:+8.1f} "
                f"{r['weight_in_pooled']:7.2f}")

    # 2. pooled deficits
    say("")
    say("== 2. pooled deficits recomputed (as 26b) ==")
    for m in MS:
        s = sum(S10[cells[(ac, m)][0]].sum() for ac in ACTS)
        n = sum(N10[cells[(ac, m)][0]].sum() for ac in ACTS)
        say(f"  {m}: {deficit(s, n):+.1f}  (sense {s:.0f} / antisense {n:.0f})")

    # 3. contrasts
    say("")
    say("== 3. POST HOC: equal-weight contrast over A1-A3 (mean of the within-third contrasts) ==")
    say("  metric             contrast                pooled (26b)                equal weight       p   A1 / A2 / A3")
    rows = []
    for lbl, col in METRICS:
        S, N = A[f"{col}_sense"].to_numpy(float), A[f"{col}_anti"].to_numpy(float)
        bs = {k: S[ix][I].sum(1) for k, (ix, I) in cells.items()}
        bn = {k: N[ix][I].sum(1) for k, (ix, I) in cells.items()}
        pt = {k: (S[ix].sum(), N[ix].sum()) for k, (ix, I) in cells.items()}
        for hi in ("M3", "M2"):
            pool_pt = (deficit(sum(pt[(ac, hi)][0] for ac in ACTS), sum(pt[(ac, hi)][1] for ac in ACTS))
                       - deficit(sum(pt[(ac, "M1")][0] for ac in ACTS), sum(pt[(ac, "M1")][1] for ac in ACTS)))
            pool_b = (deficit(sum(bs[(ac, hi)] for ac in ACTS), sum(bn[(ac, hi)] for ac in ACTS))
                      - deficit(sum(bs[(ac, "M1")] for ac in ACTS), sum(bn[(ac, "M1")] for ac in ACTS)))
            within = [deficit(*pt[(ac, hi)]) - deficit(*pt[(ac, "M1")]) for ac in ACTS]
            eq_b = np.mean([deficit(bs[(ac, hi)], bn[(ac, hi)]) - deficit(bs[(ac, "M1")], bn[(ac, "M1")]) for ac in ACTS],
                           axis=0)
            plo, phi = np.percentile(pool_b, [2.5, 97.5])
            elo, ehi = np.percentile(eq_b, [2.5, 97.5])
            r = {"set": a.set, "metric": lbl, "contrast": f"{hi}-M1", "pooled": pool_pt, "pooled_lo95": plo,
                 "pooled_hi95": phi, "equal_weight": float(np.mean(within)), "equal_lo95": elo, "equal_hi95": ehi,
                 "p_equal": pval(eq_b), "within_A1": within[0], "within_A2": within[1], "within_A3": within[2],
                 "n_genes": len(A), "n_boot": B, "seed": SEED}
            rows.append(r)
            ps = f"{pool_pt:+.1f} ({plo:+.1f}, {phi:+.1f})"
            es = f"{r['equal_weight']:+.1f} ({elo:+.1f}, {ehi:+.1f})"
            say(f"  {lbl:17s} {r['contrast']:>9s} {ps:>27s} {es:>27s} {r['p_equal']:7.2g}   "
                + " / ".join(f"{x:+.1f}" for x in within))

    # 4. rho inside activity thirds
    say("")
    say("== 4. inside each activity third: Spearman rho of methylation with host activity and with eligible nt ==")
    rr = []
    for ac in ACTS:
        d = A[A["activity"] == ac]
        r1 = spearmanr(d["meth"], d["host_tss_sense"])
        r2 = spearmanr(d["meth"], d["eligible_nt"])
        rr.append({"set": a.set, "activity_third": ac, "n_genes": len(d), "rho_meth_host_activity": r1.statistic,
                   "p_host_activity": r1.pvalue, "rho_meth_eligible_nt": r2.statistic, "p_eligible_nt": r2.pvalue})
        say(f"  {ac}: n {len(d):,}; host activity rho {r1.statistic:+.3f} (p {r1.pvalue:.2g}); "
            f"eligible nt rho {r2.statistic:+.3f} (p {r2.pvalue:.2g})")

    pre = f"{a.out_prefix}_{a.set}"
    for name, df in (("contrasts", pd.DataFrame(rows)), ("cells", pd.DataFrame(crow)), ("rho", pd.DataFrame(rr))):
        df.insert(1, "source_genes_md5", m5)
        df.to_csv(f"{pre}_{name}.tsv", sep="\t", index=False, float_format="%.10g")
    text = "\n".join(log) + "\n"
    with open(f"{pre}.log", "w") as fh:
        fh.write(text)
    print(text, end="")

    if a.check_log:
        # compare the layout-identical log with the surviving original, line by line, ignoring white space runs
        old = open(a.check_log).read().split("\n")
        new = text.split("\n")
        norm = lambda s: re.sub(r"\s+", " ", s.strip())
        diffs = [(o, n) for o, n in zip(old, new) if norm(o) != norm(n)]
        if len(old) != len(new) and not (len(old) == len(new) + 1 or len(new) == len(old) + 1):
            diffs.append((f"{len(old)} lines", f"{len(new)} lines"))
        # the first line names the input file; it may differ only in the path
        diffs = [d for d in diffs if not d[0].startswith("== 26c post hoc:")]
        if diffs:
            for o, n in diffs[:10]:
                print("OLD:", o, "\nNEW:", n)
            sys.exit(f"STOP: {len(diffs)} lines differ from {a.check_log}")
        print(f"check: every line of {a.check_log} reproduced ({len(old)} lines)")


if __name__ == "__main__":
    main()
