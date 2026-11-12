#!/usr/bin/env python3
"""Step 23c: do transposable elements explain the weaker sense start of intronic pairs?

Transposable elements (TEs) carry their own promoters, and in introns they are kept more often in
one orientation than the other. If many divergent pairs lie in TEs, and a TE's own promoter favours
one direction, the orientation bias of TEs in introns could by itself make the sense start of
intronic pairs weaker. This step asks:

  1. How many pairs overlap a TE (RepeatMasker classes SINE, LINE, LTR, DNA, Retroposon, RC;
     a trailing '?' in the class is ignored), by group (intergenic, intron_active, intron_silent)
     and class? The element of a pair is the region from 50 nt before its first start to 50 nt
     after its second; the TE assigned to a pair is the one with the largest overlap.
  2. Orientation: for intronic pairs, is the TE on the host gene's strand ('same') or on the
     other strand ('opposite')? Background: the same fraction for all TEs of that class whose
     midpoint lies in the bodies of the same host genes (GENCODE gene span).
  3. Does a TE impose a direction? For pairs in TEs, the asymmetry relative to the TE,
     DI_TE = log2((reads at the start on the TE's strand + 1) / (reads at the other start + 1)),
     summarised by class in intergenic pairs (no gene) and intronic pairs.
  4. Tables for step 02c without the pairs that overlap a TE ('noTE') and without the pairs that
     overlap any RepeatMasker element ('noRepeat'), in exactly the 02c input format (same columns,
     same row order), so that 02c can be run on them unchanged. Pairs whose arm coordinates could
     not be recovered (status != 'ok' in the step-04 table) are removed from both filtered tables
     and counted.

The step-04 table (arms_<line>.qc.tsv) and the 02c input (istart_directionality_<line>.tsv) hold
the same pairs in the same order; this is checked on chrom, mid, klasa, sens and anty (as 22f).
Arm positions are 0-based (PINTS BED); RepeatMasker coordinates are 0-based, end-exclusive (UCSC);
GTF coordinates are 1-based and converted.

Column names of the RepeatMasker table are read from its schema file (rmsk.sql), not assumed.

Usage:
  python3 23c_transposons.py --table K562=.../istart_directionality_K562.tsv --arms K562=.../arms_K562.qc.tsv \
      [--table HCT116=... --arms HCT116=...] --rmsk rmsk.txt.gz --rmsk-sql rmsk.sql \
      --gtf gencode.v47.annotation.gtf.gz --outdir results/23c [--n-boot 2000] [--seed 1]
"""
import argparse
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
EXPECTED_COLUMNS = ["chrom", "mid", "klasa", "host", "host_name", "host_signal",
                    "sens", "anty", "DI", "host_tx", "k36"]
KEYS = ["chrom", "mid", "klasa", "sens", "anty"]
ARM_COLS = ["host", "host_strand", "pos_sense", "pos_antisense", "strand_sense", "status"]
TE_CLASSES = ("SINE", "LINE", "LTR", "DNA", "Retroposon", "RC")
PAD = 50


def pairs(specs):
    out = {}
    for s in specs:
        line, path = s.split("=", 1)
        out[line] = path
    return out


def rmsk_columns(sql_path):
    sql = open(sql_path).read()
    i = sql.find("CREATE TABLE")
    if i < 0:
        sys.exit(f"{sql_path}: no CREATE TABLE statement")
    names = re.findall(r"^\s*`(\w+)`\s", sql[i:sql.find(";", i)], flags=re.M)
    need = ["genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily"]
    miss = [c for c in need if c not in names]
    if miss:
        sys.exit(f"{sql_path}: columns {miss} not in schema {names}")
    return names


def read_rmsk(path, names, chroms):
    r = pd.read_csv(path, sep="\t", header=None, names=names,
                    usecols=["genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily"],
                    dtype={"genoName": str, "strand": str, "repName": str, "repClass": str, "repFamily": str})
    r = r[r["genoName"].isin(chroms)].copy()
    r["cls"] = r["repClass"].str.rstrip("?")
    r["is_te"] = r["cls"].isin(TE_CLASSES)
    bad = ~r["strand"].isin(["+", "-"])
    if bad.any():
        sys.exit(f"{path}: {int(bad.sum())} rows with strand not + or -")
    out = {}
    for c, g in r.groupby("genoName", sort=False):
        g = g.sort_values("genoStart")
        out[c] = {k: g[k].to_numpy() for k in ("genoStart", "genoEnd", "strand", "repName", "cls", "repFamily", "is_te")}
        out[c]["mid"] = (out[c]["genoStart"] + out[c]["genoEnd"]) // 2
    return out, r


def gene_spans(gtf):
    """gene_id -> (chrom, start0, end0) from GTF gene records (1-based -> 0-based half-open)."""
    spans = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if m:
                spans[m.group(1)] = (f[0], int(f[3]) - 1, int(f[4]))
    return spans


def best_overlap(R, chrom, s, e, te_only):
    """Index of the element with the largest overlap with [s, e), or -1."""
    if chrom not in R:
        return -1, 0
    a = R[chrom]
    m = (a["genoStart"] < e) & (a["genoEnd"] > s)
    if te_only:
        m &= a["is_te"]
    idx = np.flatnonzero(m)
    if not len(idx):
        return -1, 0
    ov = np.minimum(a["genoEnd"][idx], e) - np.maximum(a["genoStart"][idx], s)
    j = int(np.argmax(ov))
    return int(idx[j]), int(ov[j])


def in_any(R, chrom, s, e, te_only):
    j, _ = best_overlap(R, chrom, s, e, te_only)
    return j >= 0


def boot_median_ci(x, n_boot, rng):
    x = np.asarray(x, float)
    if len(x) < 5:
        return np.nan, np.nan
    b = np.median(x[rng.integers(0, len(x), size=(n_boot, len(x)))], axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def boot_frac_ci_clustered(flag, unit, n_boot, rng):
    """95% CI of a proportion, resampling units (e.g. host genes)."""
    flag = np.asarray(flag, float)
    u, inv = np.unique(np.asarray(unit, str), return_inverse=True)
    s = np.bincount(inv, weights=flag, minlength=len(u))
    n = np.bincount(inv, minlength=len(u)).astype(float)
    if len(u) < 5:
        return np.nan, np.nan
    pick = rng.integers(0, len(u), size=(n_boot, len(u)))
    v = s[pick].sum(1) / n[pick].sum(1)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", required=True)
    ap.add_argument("--arms", action="append", required=True)
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--rmsk-sql", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    tables, arms = pairs(a.table), pairs(a.arms)
    if set(tables) != set(arms):
        sys.exit("--table and --arms must name the same lines")
    os.makedirs(a.outdir, exist_ok=True)

    D, Q = {}, {}
    for line in tables:
        d = pd.read_csv(tables[line], sep="\t")
        if list(d.columns) != EXPECTED_COLUMNS:
            sys.exit(f"{tables[line]}: unexpected columns {list(d.columns)}")
        q = pd.read_csv(arms[line], sep="\t", dtype={"host": str, "host_strand": str, "strand_sense": str})
        miss = [c for c in ARM_COLS if c not in q.columns]
        if miss:
            sys.exit(f"{arms[line]}: missing columns {miss}")
        if len(q) != len(d):
            sys.exit(f"{line}: {len(d)} rows in the 02c input but {len(q)} in the step-04 table")
        for k in KEYS:
            x, y = d[k].to_numpy(), q[k].to_numpy()
            same = (x == y) | (pd.isna(x) & pd.isna(y)) if x.dtype.kind == "f" else (x == y)
            if not np.all(same):
                sys.exit(f"{line}: column {k} differs between the two tables at {int((~same).sum())} rows")
        D[line], Q[line] = d, q
        print(f"[{line}] {len(d):,} pairs checked row by row against {arms[line]}")
    chroms = sorted({c for d in D.values() for c in d["chrom"].astype(str).unique()})

    names = rmsk_columns(a.rmsk_sql)
    R, rtab = read_rmsk(a.rmsk, names, chroms)
    print(f"RepeatMasker: {len(rtab):,} elements on {len(R)} chromosomes of the pairs; "
          f"TE classes {', '.join(TE_CLASSES)}: {int(rtab['is_te'].sum()):,}")
    spans = gene_spans(a.gtf)
    print(f"GTF: {len(spans):,} gene records")

    rows, filt_counts, bg_rows = [], [], []
    for line, d in D.items():
        q = Q[line]
        intr = (d["klasa"] == INTRON_CLASS).to_numpy()
        act = intr & (d["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL).to_numpy()
        group = np.where(~intr, "intergenic", np.where(act, "intron_active", "intron_silent"))
        ok = (q["status"] == "ok").to_numpy() & q["pos_sense"].notna().to_numpy() & q["pos_antisense"].notna().to_numpy()
        n = len(d)
        te_j = np.full(n, -1)
        te_ov = np.zeros(n, int)
        any_rep = np.zeros(n, bool)
        s_in, a_in = np.zeros(n, bool), np.zeros(n, bool)
        for i in np.flatnonzero(ok):
            c = str(d["chrom"].iat[i])
            ps, pa = int(q["pos_sense"].iat[i]), int(q["pos_antisense"].iat[i])
            s, e = min(ps, pa) - PAD, max(ps, pa) + PAD + 1
            te_j[i], te_ov[i] = best_overlap(R, c, s, e, True)
            any_rep[i] = in_any(R, c, s, e, False)
            s_in[i] = in_any(R, c, ps - PAD, ps + PAD + 1, True)
            a_in[i] = in_any(R, c, pa - PAD, pa + PAD + 1, True)
        for i in range(n):
            c = str(d["chrom"].iat[i])
            j = te_j[i]
            rec = {"line": line, "row": i, "chrom": c, "mid": int(d["mid"].iat[i]), "group": group[i],
                   "host": q["host"].iat[i], "host_strand": q["host_strand"].iat[i],
                   "strand_sense": q["strand_sense"].iat[i], "coords_ok": bool(ok[i]),
                   "sens": d["sens"].iat[i], "anty": d["anty"].iat[i], "DI": d["DI"].iat[i],
                   "in_te": bool(j >= 0), "in_any_repeat": bool(any_rep[i]),
                   "sense_start_in_te": bool(s_in[i]), "antisense_start_in_te": bool(a_in[i]),
                   "te_class": "", "te_family": "", "te_name": "", "te_strand": "", "te_overlap_nt": 0,
                   "te_orientation": "", "DI_te": np.nan}
            if j >= 0:
                A = R[c]
                ts = str(A["strand"][j])
                rec.update(te_class=str(A["cls"][j]), te_family=str(A["repFamily"][j]), te_name=str(A["repName"][j]),
                           te_strand=ts, te_overlap_nt=int(te_ov[i]))
                ref = rec["host_strand"] if group[i] != "intergenic" else rec["strand_sense"]
                rec["te_orientation"] = "same" if ts == ref else "opposite"
                on_te_is_sense = ts == rec["strand_sense"]
                on, other = (rec["sens"], rec["anty"]) if on_te_is_sense else (rec["anty"], rec["sens"])
                rec["DI_te"] = float(np.log2((on + 1) / (other + 1)))
            rows.append(rec)

        # filtered 02c tables
        for tag, keep in (("noTE", ok & (te_j < 0)), ("noRepeat", ok & ~any_rep)):
            out = os.path.join(a.outdir, f"istart_directionality_{line}.{tag}.tsv")
            d.loc[keep].to_csv(out, sep="\t", index=False)
            for g in ("intron_active", "intergenic"):
                m = (group == g) if g != "intergenic" else ~intr
                filt_counts.append({"line": line, "filter": tag, "group": g, "n_all": int(m.sum()),
                                    "n_no_coords": int((m & ~ok).sum()), "n_kept": int((m & keep).sum()),
                                    "pct_removed": 100 * (1 - (m & keep).sum() / m.sum())})
            print(f"[{line}] {tag}: kept {int(keep.sum()):,} of {n:,}; written {out}")

        # background orientation of TEs in the bodies of the host genes of intronic pairs
        hosts = pd.unique(q.loc[intr, "host"].astype(str))
        hs = q.loc[intr, ["host", "host_strand"]].drop_duplicates("host").set_index("host")["host_strand"].to_dict()
        found = [h for h in hosts if h in spans]
        if len(found) < 0.9 * len(hosts):
            strip = {k.split(".")[0]: v for k, v in spans.items()}
            alt = [h for h in hosts if h.split(".")[0] in strip]
            print(f"[{line}] host ids matched to GTF gene_id: exact {len(found):,}, without version {len(alt):,} "
                  f"of {len(hosts):,}")
            if len(alt) > len(found):
                spans_use, key = strip, (lambda h: h.split(".")[0])
                found = alt
            else:
                spans_use, key = spans, (lambda h: h)
        else:
            spans_use, key = spans, (lambda h: h)
        print(f"[{line}] host genes of intronic pairs: {len(hosts):,}; with a GTF span: {len(found):,}")
        cnt = defaultdict(lambda: [0, 0])
        for h in found:
            c, s0, e0 = spans_use[key(h)]
            if c not in R:
                continue
            A = R[c]
            sel = np.flatnonzero((A["mid"] >= s0) & (A["mid"] < e0) & A["is_te"])
            for j in sel:
                k = (str(A["cls"][j]), str(A["repFamily"][j]))
                cnt[k][0 if str(A["strand"][j]) == hs[h] else 1] += 1
        for (cl, fam), (same, opp) in cnt.items():
            bg_rows.append({"line": line, "te_class": cl, "te_family": fam, "n_same": same, "n_opposite": opp})

    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(a.outdir, "23c_pair_te.tsv"), sep="\t", index=False)
    F = pd.DataFrame(filt_counts)
    F.to_csv(os.path.join(a.outdir, "23c_filter_counts.tsv"), sep="\t", index=False, float_format="%.4g")
    BG = pd.DataFrame(bg_rows)
    BG.to_csv(os.path.join(a.outdir, "23c_background_orientation.tsv"), sep="\t", index=False)

    # 1. overlap by group and class
    S = []
    for (line, g), x in P[P["coords_ok"]].groupby(["line", "group"]):
        rec = {"line": line, "group": g, "n_pairs": len(x), "pct_in_te": 100 * x["in_te"].mean(),
               "pct_in_any_repeat": 100 * x["in_any_repeat"].mean(),
               "pct_sense_start_in_te": 100 * x["sense_start_in_te"].mean(),
               "pct_antisense_start_in_te": 100 * x["antisense_start_in_te"].mean()}
        for cl in TE_CLASSES:
            rec[f"pct_{cl}"] = 100 * (x["te_class"] == cl).mean()
        S.append(rec)
    S = pd.DataFrame(S)
    S.to_csv(os.path.join(a.outdir, "23c_overlap.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\n== 1. pairs whose element overlaps a TE (% of pairs with coordinates), by group and class ==")
    print(S.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    # 2. orientation of TEs at intronic pairs vs background in the same genes
    O = []
    for line in P["line"].unique():
        x = P[(P["line"] == line) & P["in_te"] & (P["group"] != "intergenic")]
        bg = BG[BG["line"] == line]
        for cl in ("all",) + TE_CLASSES:
            xx = x if cl == "all" else x[x["te_class"] == cl]
            b = bg if cl == "all" else bg[bg["te_class"] == cl]
            if not len(xx):
                continue
            opp = (xx["te_orientation"] == "opposite").to_numpy()
            lo, hi = boot_frac_ci_clustered(opp, xx["host"].astype(str), a.n_boot, rng)
            nb = b["n_same"].sum() + b["n_opposite"].sum()
            O.append({"line": line, "te_class": cl, "n_pairs": len(xx), "pct_opposite_pairs": 100 * opp.mean(),
                      "lo95": 100 * lo, "hi95": 100 * hi, "n_te_in_host_bodies": int(nb),
                      "pct_opposite_background": 100 * b["n_opposite"].sum() / nb if nb else np.nan})
    O = pd.DataFrame(O)
    O.to_csv(os.path.join(a.outdir, "23c_orientation.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\n== 2. intronic pairs in TEs: % with the TE on the strand opposite to the host gene "
          "(95% CI, resampling host genes), and the same % for all TEs in the bodies of these host genes ==")
    print(O.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    # 3. direction relative to the TE
    T = []
    for (line, g), x in P[P["in_te"]].groupby(["line", "group"]):
        for cl in ("all",) + TE_CLASSES:
            xx = x if cl == "all" else x[x["te_class"] == cl]
            if len(xx) < 5:
                continue
            lo, hi = boot_median_ci(xx["DI_te"], a.n_boot, rng)
            T.append({"line": line, "group": g, "te_class": cl, "n_pairs": len(xx),
                      "median_DI_te": float(xx["DI_te"].median()), "lo95": lo, "hi95": hi,
                      "pct_te_start_stronger": 100 * (xx["DI_te"] > 0).mean()})
    T = pd.DataFrame(T)
    T.to_csv(os.path.join(a.outdir, "23c_te_direction.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\n== 3. asymmetry relative to the TE: DI_TE = log2((reads at the start on the TE strand + 1)/(other + 1)); "
          "median with 95% CI (bootstrap over pairs) ==")
    print(T.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # 4. intronic active pairs: measured DI by TE status and orientation (no depth matching; descriptive)
    M = []
    for line in P["line"].unique():
        x = P[(P["line"] == line) & (P["group"] == "intron_active") & P["coords_ok"]]
        cat = np.where(~x["in_te"], "no TE", np.where(x["te_orientation"] == "same", "TE same strand", "TE opposite strand"))
        for c in ("no TE", "TE same strand", "TE opposite strand"):
            xx = x[cat == c]
            lo, hi = boot_median_ci(xx["DI"], a.n_boot, rng)
            M.append({"line": line, "category": c, "n_pairs": len(xx), "median_DI": float(xx["DI"].median()),
                      "lo95": lo, "hi95": hi, "median_depth": float((xx["sens"] + xx["anty"]).median())})
    M = pd.DataFrame(M)
    M.to_csv(os.path.join(a.outdir, "23c_di_by_te.tsv"), sep="\t", index=False, float_format="%.4g")
    print("\n== 4. intronic pairs in active hosts: measured DI by TE status (descriptive; depth not matched) ==")
    print(M.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n== 5. filtered tables for 02c ==")
    print(F.to_string(index=False, float_format=lambda v: f"{v:.1f}"))


if __name__ == "__main__":
    main()
