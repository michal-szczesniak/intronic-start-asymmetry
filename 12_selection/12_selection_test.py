#!/usr/bin/env python3
"""Step 12: are variants that create a start on the gene's strand selected against?

Population-genetic test of H3. Harmful variants are kept rare, so among variants
of a class under purifying selection a larger share is seen only once
(singletons, AC = 1 in gnomAD). The share of singletons also depends on how often
a site mutates (recurrent, fast-mutating sites such as CpG are less often
singletons), and in genes the two strands mutate differently (transcription-
associated mutation). So each class is compared with an expectation built from
variants that change nothing, with the same mutation in the same neighbouring
bases, read on the gene's own strand:

  expected(class) = mean over its variants of p(context),
  p(context)      = share of singletons among NEUTRAL variants of the same set
                    with that context (trinucleotide + alternative base, 192 kinds,
                    oriented to the gene's strand; plus strand for intergenic),
  excess          = observed share of singletons - expected      (the MAPS idea).

Neutral = predicted change below --neutral on both strands (step 11d scores).

Classes, within each set (protein-coding introns, lncRNA introns, intergenic):
  creates_sense_T / creates_anti_T    reference below threshold T, alternative at or above,
                                      on the gene's strand / the opposite strand
  destroys_sense_T / destroys_anti_T  the reverse
  gain_sense_strong / gain_anti_strong     change >= 0.01 on one strand, < 0.003 on the other
  gain_sense_mid / gain_anti_mid           0.003-0.01 on one strand, |change| < 0.003 on the other
  loss_sense_strong / loss_anti_strong     change <= -0.01 on one strand, |change| < 0.003 on the other
  ..._unique                          the same, outside soft-masked repeats
(For intergenic sequence "sense" is the plus strand: an orientation-free control.)

The test of H3: excess(creates_sense) - excess(creates_anti) > 0 in protein-coding
introns, absent in the intergenic control. Every difference is bootstrapped over
genes (1 Mb blocks for intergenic), resampling the same units for both members of
a pair.

Controls printed first:
  - calibration: expectation from the neutral variants of half the genes, applied
    to the neutral variants of the other half -> excess should be ~0;
  - positive control: SNVs at canonical splice sites (first two and last two
    intronic bases of canonical protein-coding transcripts) and at donor +3..+6,
    taken from the gnomAD tables of step 11b, with the protein-coding intronic
    expectation -> a clear positive excess shows the test can see selection.

Only sites with AN >= --an-frac x the largest AN are used (well-covered sites).

Usage:
  python3 12_selection_test.py --scores chr8.scores.tsv.gz chr9.scores.tsv.gz \
      --gnomad /scratch/.../gnomad_v4.1 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --out-prefix results/selection_chr8_chr9 [--n-boot 1000] [--seed 1]
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", nargs="+", required=True)
    ap.add_argument("--gnomad", required=True, help="directory with <chrom>.snv_pass.tsv.gz (step 11b)")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--an-frac", type=float, default=0.8)
    ap.add_argument("--neutral", type=float, default=0.001)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

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

    # context, repeat status, sense/antisense
    ctx = np.empty(len(d), dtype=np.int32)
    rep = np.zeros(len(d), dtype=bool)
    seqs = {}
    for c in chroms:
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        seqs[c] = raw
        m = (d["chrom"] == c).to_numpy()
        p0 = d.loc[m, "pos"].to_numpy() - 1
        alt_code = CODE[np.frombuffer("".join(d.loc[m, "alt"].astype(str)).encode(), dtype=np.uint8)]
        minus = (d.loc[m, "strand"].astype(str).to_numpy() == "-") & (d.loc[m, "set"].astype(str).to_numpy() != "intergenic")
        ctx[m] = context_codes(raw, p0, alt_code, minus)
        rep[m] = raw[p0] >= 97
    d["ctx"] = ctx
    d["rep"] = rep
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
    print(f"variants with a valid context: {len(d):,}; neutral (|change| < {a.neutral} on both strands): "
          f"{100 * d['neutral'].mean():.1f}%; singletons overall {100 * d['single'].mean():.1f}%")

    def p_table(sub):
        n = np.bincount(sub["ctx"], minlength=256).astype(float)
        s = np.bincount(sub["ctx"], weights=sub["single"], minlength=256)
        overall = s.sum() / max(1.0, n.sum())
        return np.where(n >= 20, s / np.maximum(n, 1), overall)

    sets = [s for s in ("intron_protein_coding", "intron_lncRNA", "intergenic") if (d["set"] == s).any()]
    ptab = {s: p_table(d[(d["set"] == s) & d["neutral"]]) for s in sets}

    # calibration: expectation from half of the units, applied to the neutral variants of the other half
    print("\n== calibration: neutral variants, expectation from the other half of the units (should be ~0) ==")
    for s in sets:
        sub = d[(d["set"] == s) & d["neutral"]]
        half = sub["unit"].map(lambda u: zlib.crc32(str(u).encode()) % 2).to_numpy() == 0
        if sub.loc[half, "unit"].nunique() < 10 or sub.loc[~half, "unit"].nunique() < 10:
            print(f"  {s:22s} too few units to split ({sub['unit'].nunique()})")
            continue
        pa = p_table(sub[half])
        b = sub[~half]
        print(f"  {s:22s} n = {len(b):,}: observed {b['single'].mean():.4f}, expected {pa[b['ctx']].mean():.4f}, "
              f"excess {b['single'].mean() - pa[b['ctx']].mean():+.4f}")

    # classes
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

    pairs = [(f"creates_sense_{T}", f"creates_anti_{T}") for T in THRESHOLDS] + \
            [(f"creates_sense_{T}_unique", f"creates_anti_{T}_unique") for T in THRESHOLDS] + \
            [(f"destroys_sense_{T}", f"destroys_anti_{T}") for T in THRESHOLDS] + \
            [("gain_sense_strong", "gain_anti_strong"), ("gain_sense_mid", "gain_anti_mid"),
             ("loss_sense_strong", "loss_anti_strong")]

    rows, drows = [], []

    def summarise(label, g, cls_masks, p, unit_col):
        """Per-class excess with unit bootstrap; returns dict of (n, s, e) per unit for pairs."""
        units, uidx = np.unique(g[unit_col].to_numpy(), return_inverse=True)
        W = boot_matrix(units, a.n_boot, rng)
        expv = p[g["ctx"].to_numpy()]
        sing = g["single"].to_numpy().astype(float)
        agg = {}
        for k, mk in cls_masks.items():
            mk = np.asarray(mk, dtype=bool)
            U = len(units)
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

    # positive control: canonical splice sites (gnomAD tables, protein-coding intronic expectation)
    sp, tagged = canonical_splice_positions(a.gtf, set(chroms))
    print(f"\nsplice-site control positions: {len(sp):,} ({'Ensembl_canonical transcripts' if tagged else 'all protein-coding transcripts'})")
    ctrl = []
    for c in chroms:
        gpath = os.path.join(a.gnomad, f"{c}.snv_pass.tsv.gz")
        v = pd.read_csv(gpath, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "AC", "AN"],
                        dtype={"chrom": str, "pos": np.int64, "ref": "category", "alt": "category",
                               "AC": np.int64, "AN": np.int64})
        v = v[v["pos"].isin(sp.loc[sp["chrom"] == c, "pos0"] + 1)]
        v = v[(v["AN"] >= a.an_frac * max_an) & (v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
        v["pos0"] = v["pos"] - 1
        v = v.merge(sp[sp["chrom"] == c], on=["chrom", "pos0"])
        if not len(v):
            continue
        raw = seqs[c]
        alt_code = CODE[np.frombuffer("".join(v["alt"].astype(str)).encode(), dtype=np.uint8)]
        v["ctx"] = context_codes(raw, v["pos0"].to_numpy(), alt_code, (v["strand"] == "-").to_numpy())
        v["single"] = (v["AC"] == 1).astype(np.int8)
        ctrl.append(v[v["ctx"] >= 0])
    if ctrl and "intron_protein_coding" in ptab:
        C = pd.concat(ctrl, ignore_index=True)
        summarise("positive_control", C, {k: (C["cls"] == k) for k in ("splice_canonical", "splice_donor_3to6")},
                  ptab["intron_protein_coding"], "unit")

    for s in sets:
        g = d[d["set"] == s].reset_index(drop=True)
        agg, W = summarise(s, g, classes(g), ptab[s], "unit")
        for k1, k2 in pairs:
            diff(s, agg, W, k1, k2)

    R = pd.DataFrame(rows)
    D = pd.DataFrame(drows)
    R.to_csv(a.out_prefix + "_classes.tsv", sep="\t", index=False)
    D.to_csv(a.out_prefix + "_pairs.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    fmt = lambda x: f"{x:.4g}"
    print("\n== excess of singletons by class (observed - expected from matched neutral variants) ==")
    print(R.to_string(index=False, float_format=fmt))
    print("\n== H3: excess(sense class) - excess(antisense class); > 0 = stronger selection against the sense version ==")
    print(D.to_string(index=False, float_format=fmt))
    print(f"\nCIs and p: bootstrap over genes (1 Mb blocks for intergenic), {a.n_boot} resamples, seed {a.seed}; "
          f"sites with AN >= {a.an_frac} x {max_an:,}.")
    print(f"written: {a.out_prefix}_classes.tsv, {a.out_prefix}_pairs.tsv")


if __name__ == "__main__":
    main()
