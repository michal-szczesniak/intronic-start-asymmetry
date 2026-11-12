#!/usr/bin/env python3
"""Step 12b: is the step-12 result robust to mutation-rate heterogeneity beyond the trinucleotide?

Step 12 matched each variant to neutral variants with the same trinucleotide context.
Two known sources of mutation-rate variation are not captured by that:
  - CpG sites: their mutation rate depends on methylation, and unmethylated CpGs
    (CpG islands, often near predicted start sites) mutate less, which by itself
    raises the share of singletons;
  - regional base composition (local GC), which tracks many other rate differences.
This script loads the data once and repeats the step-12 analysis in four settings:
  base        exactly step 12 (same random seed; its tables must be byte-identical
              to the step-12 output, which checks this script)
  noCpG       variants at CpG sites removed everywhere (classes, neutral reference,
              splice-site control)
  gc5         expectation matched on trinucleotide context AND local GC (+/-50 nt,
              5 bins with edges at the quintiles of all scored variants)
  noCpG_gc5   both
Readout: a compact table across settings for the positive control, the H3 pairs
(creates / gains, sense vs antisense) and the "destroys" classes, whose step-12
excess was suspected to be a CpG-methylation artefact.

Usage: as step 12, plus --out-dir results/robust
"""
import argparse
import gzip
import os
import re
import sys
import zlib

import numpy as np
import pandas as pd

THRESHOLDS = (0.0439, 0.1432)
STRONG, MID, OTHER = 0.01, 0.003, 0.003
BLOCK = 1_000_000
CODE = np.full(256, -1, dtype=np.int16)
for _i, _b in enumerate("ACGT"):
    CODE[ord(_b)] = _i
    CODE[ord(_b.lower())] = _i
BASES = np.array(list("ACGT"))


def context_codes(up, pos0, alt_code, minus):
    """Trinucleotide + alternative base as an integer 0..255 (tri*4 + alt), read on
    the plus strand, or on the minus strand where `minus` is True. -1 where any base is N."""
    b1, b2, b3 = CODE[up[pos0 - 1]], CODE[up[pos0]], CODE[up[pos0 + 1]]
    bad = (b1 < 0) | (b2 < 0) | (b3 < 0) | (alt_code < 0)
    c1, c2, c3, ca = b1.astype(np.int32), b2.astype(np.int32), b3.astype(np.int32), alt_code.astype(np.int32)
    r1, r2, r3, ra = 3 - c3, 3 - c2, 3 - c1, 3 - ca     # reverse complement: A<->T (0<->3), C<->G (1<->2)
    t1, t2, t3, ta = (np.where(minus, r1, c1), np.where(minus, r2, c2),
                      np.where(minus, r3, c3), np.where(minus, ra, ca))
    code = ((t1 * 4 + t2) * 4 + t3) * 4 + ta
    return np.where(bad, -1, code)


def code_name(c):
    t, a = divmod(int(c), 4)
    t1, r = divmod(t, 16)
    t2, t3 = divmod(r, 4)
    return f"{BASES[t1]}{BASES[t2]}{BASES[t3]}>{BASES[a]}"


def canonical_splice_positions(gtf, chroms):
    """0-based positions of canonical splice sites (+1,+2 donor; -1,-2 acceptor) and of
    donor +3..+6, with the gene strand and gene id, from canonical protein-coding
    transcripts (all protein-coding transcripts if the GTF has no Ensembl_canonical tags)."""
    tx = {}
    any_tag = False
    re_t = re.compile(r'transcript_id "([^"]+)"')
    re_g = re.compile(r'gene_id "([^"]+)"')
    re_tt = re.compile(r'transcript_type "([^"]+)"')
    re_gt = re.compile(r'gene_type "([^"]+)"')
    with gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in chroms:
                continue
            a = f[8]
            gt = re_gt.search(a)
            tt = re_tt.search(a)
            if (gt and gt.group(1) != "protein_coding") or (tt and tt.group(1) != "protein_coding"):
                continue
            canon = 'tag "Ensembl_canonical"' in a
            any_tag |= canon
            key = re_t.search(a).group(1) if re_t.search(a) else re_g.search(a).group(1)
            d = tx.setdefault(key, {"chrom": f[0], "strand": f[6], "gene": re_g.search(a).group(1),
                                    "canon": canon, "exons": []})
            d["exons"].append((int(f[3]) - 1, int(f[4])))
    rows = []
    for d in tx.values():
        if any_tag and not d["canon"]:
            continue
        ex = sorted(d["exons"])
        for (s1, e1), (s2, e2) in zip(ex[:-1], ex[1:]):
            s, e = e1, s2                       # intron [s, e)
            if e - s < 20:
                continue
            if d["strand"] == "+":
                canon_pos = [s, s + 1, e - 2, e - 1]
                donor36 = [s + 2, s + 3, s + 4, s + 5]
            else:
                canon_pos = [e - 1, e - 2, s, s + 1]
                donor36 = [e - 3, e - 4, e - 5, e - 6]
            rows += [(d["chrom"], p, d["strand"], d["gene"], "splice_canonical") for p in canon_pos]
            rows += [(d["chrom"], p, d["strand"], d["gene"], "splice_donor_3to6") for p in donor36]
    df = pd.DataFrame(rows, columns=["chrom", "pos0", "strand", "unit", "cls"]).drop_duplicates(["chrom", "pos0", "cls"])
    return df, any_tag


def boot_matrix(units, n_boot, rng):
    """(n_boot, U) resampling counts over units."""
    U = len(units)
    idx = rng.integers(0, U, size=(n_boot, U))
    W = np.zeros((n_boot, U), dtype=np.int32)
    for b in range(n_boot):
        W[b] = np.bincount(idx[b], minlength=U)
    return W



SETTINGS = (("base", False, 1), ("noCpG", True, 1), ("gc5", False, 5), ("noCpG_gc5", True, 5))
GC_HALF = 50
G_, C_ = ord("G"), ord("C")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", nargs="+", required=True)
    ap.add_argument("--gnomad", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--an-frac", type=float, default=0.8)
    ap.add_argument("--neutral", type=float, default=0.001)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)

    parts = []
    for f in a.scores:
        d = pd.read_csv(f, sep="\t", keep_default_na=False, na_values=[],
                        dtype={"chrom": str, "ref": str, "alt": str, "set": str, "gene_id": str, "strand": str})
        for k in ("chrom", "ref", "alt", "set", "gene_id", "strand"):
            d[k] = d[k].astype("category")
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    for k in ("chrom", "ref", "alt", "set", "gene_id", "strand"):
        d[k] = d[k].astype(str).astype("category")
    n0 = len(d)
    d = d.drop_duplicates(["chrom", "pos", "alt"]).reset_index(drop=True)
    chroms = sorted(d["chrom"].astype(str).unique())
    max_an = int(d["AN"].max())
    d = d[d["AN"] >= a.an_frac * max_an].reset_index(drop=True)
    print(f"scores: {n0:,} rows, {len(d):,} after removing duplicates and sites with AN < {a.an_frac} x {max_an:,}")

    ctx = np.empty(len(d), dtype=np.int32)
    rep = np.zeros(len(d), dtype=bool)
    cpg = np.zeros(len(d), dtype=bool)
    gc = np.zeros(len(d), dtype=np.float32)
    seqs, gcs = {}, {}

    def site_props(raw, cs, p0):
        up_ = np.where(raw >= 97, raw - 32, raw)
        b1, b2, b3 = up_[p0 - 1], up_[p0], up_[p0 + 1]
        is_cpg = ((b2 == C_) & (b3 == G_)) | ((b2 == G_) & (b1 == C_))
        g = (cs[np.minimum(p0 + GC_HALF + 1, len(cs) - 1)] - cs[np.maximum(p0 - GC_HALF, 0)]) / (2 * GC_HALF + 1)
        return is_cpg, g.astype(np.float32)

    for c in chroms:
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        seqs[c] = raw
        isgc = (raw == G_) | (raw == C_) | (raw == G_ + 32) | (raw == C_ + 32)
        cs = np.concatenate([[0], np.cumsum(isgc, dtype=np.int32)])
        gcs[c] = cs
        m = (d["chrom"] == c).to_numpy()
        p0 = d.loc[m, "pos"].to_numpy() - 1
        alt_code = CODE[np.frombuffer("".join(d.loc[m, "alt"].astype(str)).encode(), dtype=np.uint8)]
        minus = (d.loc[m, "strand"].astype(str).to_numpy() == "-") & (d.loc[m, "set"].astype(str).to_numpy() != "intergenic")
        ctx[m] = context_codes(raw, p0, alt_code, minus)
        rep[m] = raw[p0] >= 97
        cpg[m], gc[m] = site_props(raw, cs, p0)
    d["ctx"], d["rep"], d["cpg"], d["gc"] = ctx, rep, cpg, gc
    d = d[d["ctx"] >= 0].reset_index(drop=True)
    plus_gene = (d["strand"].astype(str) == "+") | (d["set"].astype(str) == "intergenic")
    for side, (cp, cm) in {"ref": ("ref_plus", "ref_minus"), "alt": ("alt_plus", "alt_minus")}.items():
        d[f"{side}_s"] = np.where(plus_gene, d[cp], d[cm])
        d[f"{side}_a"] = np.where(plus_gene, d[cm], d[cp])
    d["ds"] = d["alt_s"] - d["ref_s"]
    d["da"] = d["alt_a"] - d["ref_a"]
    d["single"] = (d["AC"] == 1).astype(np.int8)
    d["neutral"] = (d["ds"].abs() < a.neutral) & (d["da"].abs() < a.neutral)
    d["unit"] = np.where(d["set"].astype(str) == "intergenic",
                         d["chrom"].astype(str) + ":" + (d["pos"] // BLOCK).astype(str), d["gene_id"].astype(str))
    edges = np.quantile(d["gc"], [0.2, 0.4, 0.6, 0.8])
    print(f"variants with a valid context: {len(d):,}; at CpG sites: {100 * d['cpg'].mean():.1f}%; "
          f"local GC quintile edges: " + ", ".join(f"{e:.3f}" for e in edges))

    sp, tagged = canonical_splice_positions(a.gtf, set(chroms))
    ctrl = []
    for c in chroms:
        v = pd.read_csv(os.path.join(a.gnomad, f"{c}.snv_pass.tsv.gz"), sep="\t", header=None,
                        names=["chrom", "pos", "ref", "alt", "AC", "AN"],
                        dtype={"chrom": str, "pos": np.int64, "ref": "category", "alt": "category",
                               "AC": np.int64, "AN": np.int64})
        v = v[v["pos"].isin(sp.loc[sp["chrom"] == c, "pos0"] + 1)]
        v = v[(v["AN"] >= a.an_frac * max_an) & (v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
        v["pos0"] = v["pos"] - 1
        v = v.merge(sp[sp["chrom"] == c], on=["chrom", "pos0"])
        if not len(v):
            continue
        alt_code = CODE[np.frombuffer("".join(v["alt"].astype(str)).encode(), dtype=np.uint8)]
        v["ctx"] = context_codes(seqs[c], v["pos0"].to_numpy(), alt_code, (v["strand"] == "-").to_numpy())
        v["single"] = (v["AC"] == 1).astype(np.int8)
        v["cpg"], v["gc"] = site_props(seqs[c], gcs[c], v["pos0"].to_numpy())
        ctrl.append(v[v["ctx"] >= 0])
    C_all = pd.concat(ctrl, ignore_index=True) if ctrl else None

    pairs = [(f"creates_sense_{T}", f"creates_anti_{T}") for T in THRESHOLDS] + \
            [(f"creates_sense_{T}_unique", f"creates_anti_{T}_unique") for T in THRESHOLDS] + \
            [(f"destroys_sense_{T}", f"destroys_anti_{T}") for T in THRESHOLDS] + \
            [("gain_sense_strong", "gain_anti_strong"), ("gain_sense_mid", "gain_anti_mid"),
             ("loss_sense_strong", "loss_anti_strong")]

    def classes(g):
        ds, da, rs, ra, as_, aa = g["ds"], g["da"], g["ref_s"], g["ref_a"], g["alt_s"], g["alt_a"]
        out = {"neutral": g["neutral"]}
        for T in THRESHOLDS:
            out[f"creates_sense_{T}"] = (rs < T) & (as_ >= T)
            out[f"creates_anti_{T}"] = (ra < T) & (aa >= T)
            out[f"destroys_sense_{T}"] = (rs >= T) & (as_ < T)
            out[f"destroys_anti_{T}"] = (ra >= T) & (aa < T)
        out["gain_sense_strong"] = (ds >= STRONG) & (da < OTHER)
        out["gain_anti_strong"] = (da >= STRONG) & (ds < OTHER)
        out["gain_sense_mid"] = (ds >= MID) & (ds < STRONG) & (da.abs() < OTHER)
        out["gain_anti_mid"] = (da >= MID) & (da < STRONG) & (ds.abs() < OTHER)
        out["loss_sense_strong"] = (ds <= -STRONG) & (da.abs() < OTHER)
        out["loss_anti_strong"] = (da <= -STRONG) & (ds.abs() < OTHER)
        for T in THRESHOLDS:
            out[f"creates_sense_{T}_unique"] = out[f"creates_sense_{T}"] & ~g["rep"]
            out[f"creates_anti_{T}_unique"] = out[f"creates_anti_{T}"] & ~g["rep"]
        return out

    compact = []
    for name, no_cpg, K in SETTINGS:
        rng = np.random.default_rng(a.seed)   # same stream per setting: "base" reproduces step 12 exactly
        D = d[~d["cpg"]].reset_index(drop=True) if no_cpg else d
        D = D.assign(key=D["ctx"].to_numpy() * K + (np.searchsorted(edges, D["gc"].to_numpy()) if K > 1 else 0))
        print(f"\n######## setting {name}: {len(D):,} variants", flush=True)

        def p_table(sub):
            n = np.bincount(sub["key"], minlength=256 * K).astype(float)
            s = np.bincount(sub["key"], weights=sub["single"], minlength=256 * K)
            overall = s.sum() / max(1.0, n.sum())
            return np.where(n >= 20, s / np.maximum(n, 1), overall)

        sets = [s for s in ("intron_protein_coding", "intron_lncRNA", "intergenic") if (D["set"] == s).any()]
        ptab = {s: p_table(D[(D["set"] == s) & D["neutral"]]) for s in sets}
        for s in sets:
            sub = D[(D["set"] == s) & D["neutral"]]
            half = sub["unit"].map(lambda u: zlib.crc32(str(u).encode()) % 2).to_numpy() == 0
            if sub.loc[half, "unit"].nunique() < 10 or sub.loc[~half, "unit"].nunique() < 10:
                continue
            pa = p_table(sub[half])
            b = sub[~half]
            print(f"  calibration {s:22s} excess {b['single'].mean() - pa[b['key']].mean():+.4f}")

        rows, drows = [], []

        def summarise(label, g, cls_masks, p):
            units, uidx = np.unique(g["unit"].to_numpy(), return_inverse=True)
            U = len(units)
            idx = rng.integers(0, U, size=(a.n_boot, U))
            W = np.zeros((a.n_boot, U), dtype=np.int32)
            for b in range(a.n_boot):
                W[b] = np.bincount(idx[b], minlength=U)
            expv = p[g["key"].to_numpy()]
            sing = g["single"].to_numpy().astype(float)
            agg = {}
            for k, mk in cls_masks.items():
                mk = np.asarray(mk, dtype=bool)
                n_u = np.bincount(uidx[mk], minlength=U).astype(float)
                s_u = np.bincount(uidx[mk], weights=sing[mk], minlength=U)
                e_u = np.bincount(uidx[mk], weights=expv[mk], minlength=U)
                agg[k] = (n_u, s_u, e_u)
                n = n_u.sum()
                if n == 0:
                    continue
                ex = (s_u.sum() - e_u.sum()) / n
                nb, sb, eb = W @ n_u, W @ s_u, W @ e_u
                ok = nb > 0
                exb = (sb[ok] - eb[ok]) / nb[ok]
                rows.append({"set": label, "class": k, "n": int(n), "units": int((n_u > 0).sum()),
                             "obs_singleton": s_u.sum() / n, "exp_singleton": e_u.sum() / n, "excess": ex,
                             "excess_lo95": float(np.percentile(exb, 2.5)) if len(exb) else np.nan,
                             "excess_hi95": float(np.percentile(exb, 97.5)) if len(exb) else np.nan})
            return agg, W

        def diff(label, agg, W, k1, k2):
            if k1 not in agg or k2 not in agg:
                return
            (n1, s1, e1), (n2, s2, e2) = agg[k1], agg[k2]
            if n1.sum() == 0 or n2.sum() == 0:
                return
            d0 = (s1.sum() - e1.sum()) / n1.sum() - (s2.sum() - e2.sum()) / n2.sum()
            nb1, nb2 = W @ n1, W @ n2
            ok = (nb1 > 0) & (nb2 > 0)
            db = ((W @ s1 - W @ e1)[ok] / nb1[ok]) - ((W @ s2 - W @ e2)[ok] / nb2[ok])
            p = 2 * min((db <= 0).mean(), (db >= 0).mean())
            drows.append({"set": label, "sense_class": k1, "anti_class": k2, "n_sense": int(n1.sum()),
                          "n_anti": int(n2.sum()), "excess_diff": d0, "diff_lo95": float(np.percentile(db, 2.5)),
                          "diff_hi95": float(np.percentile(db, 97.5)), "p_boot": max(p, 1 / a.n_boot)})

        if C_all is not None and "intron_protein_coding" in ptab:
            C = C_all[~C_all["cpg"]].reset_index(drop=True) if no_cpg else C_all
            C = C.assign(key=C["ctx"].to_numpy() * K + (np.searchsorted(edges, C["gc"].to_numpy()) if K > 1 else 0))
            summarise("positive_control", C, {k: (C["cls"] == k) for k in ("splice_canonical", "splice_donor_3to6")},
                      ptab["intron_protein_coding"])
        for s in sets:
            g = D[D["set"] == s].reset_index(drop=True)
            agg, W = summarise(s, g, classes(g), ptab[s])
            for k1, k2 in pairs:
                diff(s, agg, W, k1, k2)
        R, Dp = pd.DataFrame(rows), pd.DataFrame(drows)
        R.to_csv(os.path.join(a.out_dir, f"{name}_classes.tsv"), sep="\t", index=False)
        Dp.to_csv(os.path.join(a.out_dir, f"{name}_pairs.tsv"), sep="\t", index=False)

        def get(df, s, col, k, field):
            x = df[(df["set"] == s) & (df[col] == k)]
            return x.iloc[0][field] if len(x) else np.nan

        for s in ("positive_control", "intron_protein_coding", "intergenic"):
            for k in ("splice_canonical", "splice_donor_3to6", "creates_sense_0.0439", "creates_anti_0.0439",
                      "gain_sense_strong", "gain_anti_strong", "destroys_sense_0.0439", "destroys_anti_0.0439",
                      "loss_sense_strong", "loss_anti_strong"):
                if (s == "positive_control") != k.startswith("splice"):
                    continue
                if np.isnan(get(R, s, "class", k, "excess")):
                    continue
                compact.append({"setting": name, "set": s, "quantity": f"excess {k}", "n": get(R, s, "class", k, "n"),
                                "value": get(R, s, "class", k, "excess"),
                                "lo95": get(R, s, "class", k, "excess_lo95"), "hi95": get(R, s, "class", k, "excess_hi95"),
                                "p_boot": np.nan})
            if s == "positive_control":
                continue
            for k1 in ("creates_sense_0.0439", "creates_sense_0.1432", "gain_sense_strong", "destroys_sense_0.0439"):
                if np.isnan(get(Dp, s, "sense_class", k1, "excess_diff")):
                    continue
                compact.append({"setting": name, "set": s, "quantity": f"diff {k1} - anti",
                                "n": get(Dp, s, "sense_class", k1, "n_sense") + get(Dp, s, "sense_class", k1, "n_anti"),
                                "value": get(Dp, s, "sense_class", k1, "excess_diff"),
                                "lo95": get(Dp, s, "sense_class", k1, "diff_lo95"),
                                "hi95": get(Dp, s, "sense_class", k1, "diff_hi95"),
                                "p_boot": get(Dp, s, "sense_class", k1, "p_boot")})

    Cm = pd.DataFrame(compact)
    Cm.to_csv(os.path.join(a.out_dir, "compact.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print("\n== compact comparison across settings ==")
    for (s, q), g in Cm.groupby(["set", "quantity"], sort=False):
        line = "; ".join(f"{r.setting} {r.value:+.4f} [{r.lo95:+.4f}, {r.hi95:+.4f}]" + (f" p={r.p_boot:.3g}" if r.p_boot == r.p_boot else "")
                         for r in g.itertuples())
        print(f"  {s:22s} {q:38s} {line}")
    print(f"\nwritten: {a.out_dir}/<setting>_classes.tsv, <setting>_pairs.tsv, compact.tsv")


if __name__ == "__main__":
    main()
