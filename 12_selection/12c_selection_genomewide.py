#!/usr/bin/env python3
"""Step 12c: the selection test of steps 12/12b on many chromosomes, read one at a time.

Steps 12 and 12b load every scored variant into one table. That works for chr8+chr9
(23 M variants) but not for the whole genome (~200 M). This script computes exactly
the same quantities, in the same four settings as 12b (base, noCpG, gc5, noCpG_gc5),
with the same random seed, but it reads the score files one chromosome at a time
(in chunks) and keeps only a compact record per variant:
  AN, set, unit (gene, or 1 Mb block for intergenic), trinucleotide+alt context,
  local GC count (+/-50 nt), CpG / repeat / singleton flags, and a bit mask of the
  step-12 classes the variant belongs to (computed from the scores on the fly).
About 17 bytes per variant, so ~4 GB for the genome.

Checks built in:
  - On chr8+chr9 the four <setting>_classes.tsv, <setting>_pairs.tsv and compact.tsv
    must be byte-identical to the step-12b output (same rows, same order of summation,
    same random stream). Run that comparison before the genome-wide run.
  - The splice-site positive control and the half-split calibration, as before.

New output (genome-wide only makes sense with it):
  groups.tsv   the H3 pairs and the positive control per chromosome and per group of
               chromosomes, each with its own bootstrap over units (--n-boot-group):
                 heldout_test = chr8, chr9 (Puffin-D test chromosomes),
                 heldout_valid = chr10 (validation), training = all others.
               Puffin-D saw the reference sequence of the training chromosomes during
               training; if that mattered, the result should differ between groups.
  calibration.tsv, loaded.tsv (per-chromosome row counts).

Usage:
  python3 12c_selection_genomewide.py --scores-dir $SCRATCH/11d_scores --chroms chr8,chr9 \
      --gnomad $SCRATCH/gnomad_v4.1 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --out-dir results/gw_check
Files read: <scores-dir>/<chrom>.scores.tsv.gz and <gnomad>/<chrom>.snv_pass.tsv.gz for
each chromosome named in --chroms (no globbing).
"""
import argparse
import gzip
import os
import re
import resource
import sys
import time
import zlib

import numpy as np
import pandas as pd

THRESHOLDS = (0.0439, 0.1432)
STRONG, MID, OTHER = 0.01, 0.003, 0.003
BLOCK = 1_000_000
SETTINGS = (("base", False, 1), ("noCpG", True, 1), ("gc5", False, 5), ("noCpG_gc5", True, 5))
GC_HALF = 50
G_, C_ = ord("G"), ord("C")
SETS = ("intron_protein_coding", "intron_lncRNA", "intergenic")
CODE = np.full(256, -1, dtype=np.int16)
for _i, _b in enumerate("ACGT"):
    CODE[ord(_b)] = _i
    CODE[ord(_b.lower())] = _i

# class order = the order of step 12's classes() dict (it fixes the row order of the tables)
CLASSES = ["neutral"]
for _T in THRESHOLDS:
    CLASSES += [f"creates_sense_{_T}", f"creates_anti_{_T}", f"destroys_sense_{_T}", f"destroys_anti_{_T}"]
CLASSES += ["gain_sense_strong", "gain_anti_strong", "gain_sense_mid", "gain_anti_mid",
            "loss_sense_strong", "loss_anti_strong"]
for _T in THRESHOLDS:
    CLASSES += [f"creates_sense_{_T}_unique", f"creates_anti_{_T}_unique"]
BIT = {k: np.uint32(1 << i) for i, k in enumerate(CLASSES)}
PAIRS = [(f"creates_sense_{T}", f"creates_anti_{T}") for T in THRESHOLDS] + \
        [(f"creates_sense_{T}_unique", f"creates_anti_{T}_unique") for T in THRESHOLDS] + \
        [(f"destroys_sense_{T}", f"destroys_anti_{T}") for T in THRESHOLDS] + \
        [("gain_sense_strong", "gain_anti_strong"), ("gain_sense_mid", "gain_anti_mid"),
         ("loss_sense_strong", "loss_anti_strong")]
CTRL = ("splice_canonical", "splice_donor_3to6")
GROUPS = (("heldout_test", ("chr8", "chr9")), ("heldout_valid", ("chr10",)))
F_CPG, F_REP, F_SINGLE = np.uint8(1), np.uint8(2), np.uint8(4)
UPPER = np.where(np.arange(256) >= 97, np.arange(256) - 32, np.arange(256)).astype(np.uint8)   # 12b: raw >= 97 -> raw - 32


def context_codes(up, pos0, alt_code, minus):
    """Trinucleotide + alternative base as an integer 0..255 (tri*4 + alt), read on
    the plus strand, or on the minus strand where `minus` is True. -1 where any base is N."""
    b1, b2, b3 = CODE[up[pos0 - 1]], CODE[up[pos0]], CODE[up[pos0 + 1]]
    bad = (b1 < 0) | (b2 < 0) | (b3 < 0) | (alt_code < 0)
    c1, c2, c3, ca = b1.astype(np.int32), b2.astype(np.int32), b3.astype(np.int32), alt_code.astype(np.int32)
    r1, r2, r3, ra = 3 - c3, 3 - c2, 3 - c1, 3 - ca
    t1, t2, t3, ta = (np.where(minus, r1, c1), np.where(minus, r2, c2),
                      np.where(minus, r3, c3), np.where(minus, ra, ca))
    code = ((t1 * 4 + t2) * 4 + t3) * 4 + ta
    return np.where(bad, -1, code)


def site_props(raw, cs, p0):
    """CpG status of the site and the GC count in the 101 nt around it (0..101), as 12b."""
    b1, b2, b3 = UPPER[raw[p0 - 1]], UPPER[raw[p0]], UPPER[raw[p0 + 1]]
    is_cpg = ((b2 == C_) & (b3 == G_)) | ((b2 == G_) & (b1 == C_))
    cnt = cs[np.minimum(p0 + GC_HALF + 1, len(cs) - 1)] - cs[np.maximum(p0 - GC_HALF, 0)]
    return is_cpg, cnt.astype(np.uint8)


def canonical_splice_positions(gtf, chroms):
    """As in steps 12/12b (unchanged)."""
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
            s, e = e1, s2
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


def mem_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2


class Units:
    """Global unit ids keyed by (set code, name); remembers each unit's chromosome."""

    def __init__(self):
        self.ids, self.names, self.chrom = {}, [], []

    def get(self, setc, name, chrom):
        k = (setc, name)
        u = self.ids.get(k)
        if u is None:
            u = self.ids[k] = len(self.names)
            self.names.append(name)
            self.chrom.append(chrom)
        return u


def load_scores(path, chrom, raw, cs, units, chunk, neutral):
    """Compact per-variant records of one chromosome's score file (duplicates removed,
    keeping the first, as step 12; variants with an N in the context dropped)."""
    use = ["chrom", "pos", "alt", "AC", "AN", "set", "gene_id", "strand", "ref_plus", "ref_minus", "alt_plus", "alt_minus"]
    dt = {"chrom": str, "alt": str, "set": str, "gene_id": str, "strand": str,
          "pos": np.int64, "AC": np.int64, "AN": np.int64}
    parts = {k: [] for k in ("dkey", "AN", "setc", "uid", "ctx", "gc", "flags", "cls")}
    n_in = 0
    for ch in pd.read_csv(path, sep="\t", usecols=use, dtype=dt, keep_default_na=False, na_values=[], chunksize=chunk):
        n_in += len(ch)
        cc = ch["chrom"].to_numpy(dtype=object)
        if not (cc == chrom).all():
            sys.exit(f"STOP: {path} contains chromosome(s) other than {chrom}: {sorted(set(cc) - {chrom})}")
        pos = ch["pos"].to_numpy()
        p0 = pos - 1
        alt = ch["alt"].to_numpy(dtype=object)
        alt1 = np.array([x if len(x) == 1 else "N" for x in alt], dtype="S1") if any(len(x) != 1 for x in alt) \
            else np.frombuffer("".join(alt).encode(), dtype="S1")
        alt_code = CODE[alt1.view(np.uint8)]
        setn = ch["set"].to_numpy(dtype=object)
        setc = np.full(len(ch), len(SETS), dtype=np.int8)        # len(SETS) = any other set name
        for i, s in enumerate(SETS):
            setc[setn == s] = i
        is_ig = setc == SETS.index("intergenic")
        strand = ch["strand"].to_numpy(dtype=object)
        minus = (strand == "-") & ~is_ig
        ctx = context_codes(raw, p0, alt_code, minus)
        rep = raw[p0] >= 97
        cpg, gcc = site_props(raw, cs, p0)
        plus_gene = (strand == "+") | is_ig
        rp, rm = ch["ref_plus"].to_numpy(), ch["ref_minus"].to_numpy()
        ap_, am = ch["alt_plus"].to_numpy(), ch["alt_minus"].to_numpy()
        rs, ra = np.where(plus_gene, rp, rm), np.where(plus_gene, rm, rp)
        as_, aa = np.where(plus_gene, ap_, am), np.where(plus_gene, am, ap_)
        ds, da = as_ - rs, aa - ra
        m = {"neutral": (np.abs(ds) < neutral) & (np.abs(da) < neutral)}
        for T in THRESHOLDS:
            m[f"creates_sense_{T}"] = (rs < T) & (as_ >= T)
            m[f"creates_anti_{T}"] = (ra < T) & (aa >= T)
            m[f"destroys_sense_{T}"] = (rs >= T) & (as_ < T)
            m[f"destroys_anti_{T}"] = (ra >= T) & (aa < T)
        m["gain_sense_strong"] = (ds >= STRONG) & (da < OTHER)
        m["gain_anti_strong"] = (da >= STRONG) & (ds < OTHER)
        m["gain_sense_mid"] = (ds >= MID) & (ds < STRONG) & (np.abs(da) < OTHER)
        m["gain_anti_mid"] = (da >= MID) & (da < STRONG) & (np.abs(ds) < OTHER)
        m["loss_sense_strong"] = (ds <= -STRONG) & (np.abs(da) < OTHER)
        m["loss_anti_strong"] = (da <= -STRONG) & (np.abs(ds) < OTHER)
        for T in THRESHOLDS:
            m[f"creates_sense_{T}_unique"] = m[f"creates_sense_{T}"] & ~rep
            m[f"creates_anti_{T}_unique"] = m[f"creates_anti_{T}"] & ~rep
        cls = np.zeros(len(ch), dtype=np.uint32)
        for k in CLASSES:
            cls |= np.where(m[k], BIT[k], np.uint32(0)).astype(np.uint32)
        # units: gene id (gene sets) or chrom:1Mb-block (intergenic), as step 12
        gid = ch["gene_id"].to_numpy(dtype=object)
        gcode, guniq = pd.factorize(gid)
        G = len(guniq) + 1
        blk = pos // BLOCK
        comb = np.where(is_ig, -1 - blk, setc.astype(np.int64) * G + gcode)
        cu, inv = np.unique(comb, return_inverse=True)
        lut = np.empty(len(cu), dtype=np.int32)
        for j, v in enumerate(cu):
            if v < 0:
                lut[j] = units.get(SETS.index("intergenic"), f"{chrom}:{-1 - v}", chrom)
            else:
                sc, gc_ = divmod(int(v), G)
                lut[j] = units.get(sc, str(guniq[gc_]), chrom)
        uid = lut[inv.ravel()]
        flags = (np.where(cpg, F_CPG, 0) | np.where(rep, F_REP, 0) | np.where(ch["AC"].to_numpy() == 1, F_SINGLE, 0)).astype(np.uint8)
        parts["dkey"].append(pos * 5 + (alt_code.astype(np.int64) + 1))
        parts["AN"].append(ch["AN"].to_numpy().astype(np.int32))
        parts["setc"].append(setc)
        parts["uid"].append(uid)
        parts["ctx"].append(ctx.astype(np.int16))
        parts["gc"].append(gcc)
        parts["flags"].append(flags)
        parts["cls"].append(cls)
    if n_in == 0:
        sys.exit(f"STOP: {path} has no rows")
    R = {k: np.concatenate(v) for k, v in parts.items()}
    _, first = np.unique(R.pop("dkey"), return_index=True)
    first.sort()
    max_an = int(R["AN"].max())
    keep = first[R["ctx"][first] >= 0]
    return {k: v[keep] for k, v in R.items()}, n_in, max_an, len(first)


def load_control(path, chrom, sp_c, raw, cs, chunk):
    """gnomAD SNVs at the splice-control positions of one chromosome (AN filter applied later)."""
    if not len(sp_c):
        return None
    want = np.unique(sp_c["pos0"].to_numpy() + 1)
    got = []
    for ch in pd.read_csv(path, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "AC", "AN"],
                          dtype={"chrom": str, "pos": np.int64, "ref": str, "alt": str, "AC": np.int64, "AN": np.int64},
                          keep_default_na=False, na_values=[], chunksize=chunk):
        got.append(ch[np.isin(ch["pos"].to_numpy(), want)])
    v = pd.concat(got, ignore_index=True)
    v = v[(v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
    v = v.assign(pos0=v["pos"] - 1)
    v = v.merge(sp_c, on=["chrom", "pos0"])
    if not len(v):
        return None
    p0 = v["pos0"].to_numpy()
    alt_code = CODE[np.frombuffer("".join(v["alt"].astype(str)).encode(), dtype=np.uint8)]
    ctx = context_codes(raw, p0, alt_code, (v["strand"] == "-").to_numpy())
    cpg, gcc = site_props(raw, cs, p0)
    v = v.assign(ctx=ctx, single=(v["AC"] == 1).astype(np.int8), cpg=cpg, gcc=gcc)
    return v[v["ctx"] >= 0][["chrom", "pos", "AN", "unit", "cls", "ctx", "single", "cpg", "gcc"]]


def boot_W(U, n_boot, rng):
    idx = rng.integers(0, U, size=(n_boot, U))
    W = np.zeros((n_boot, U), dtype=np.int32)
    for b in range(n_boot):
        W[b] = np.bincount(idx[b], minlength=U)
    return W


def excess_stats(n_u, s_u, e_u, W):
    n = n_u.sum()
    ex = (s_u.sum() - e_u.sum()) / n
    nb, sb, eb = W @ n_u, W @ s_u, W @ e_u
    ok = nb > 0
    exb = (sb[ok] - eb[ok]) / nb[ok]
    return {"n": int(n), "units": int((n_u > 0).sum()), "obs_singleton": s_u.sum() / n,
            "exp_singleton": e_u.sum() / n, "excess": ex,
            "excess_lo95": float(np.percentile(exb, 2.5)) if len(exb) else np.nan,
            "excess_hi95": float(np.percentile(exb, 97.5)) if len(exb) else np.nan}


def diff_stats(A1, A2, W, n_boot):
    (n1, s1, e1), (n2, s2, e2) = A1, A2
    d0 = (s1.sum() - e1.sum()) / n1.sum() - (s2.sum() - e2.sum()) / n2.sum()
    nb1, nb2 = W @ n1, W @ n2
    ok = (nb1 > 0) & (nb2 > 0)
    db = ((W @ s1 - W @ e1)[ok] / nb1[ok]) - ((W @ s2 - W @ e2)[ok] / nb2[ok])
    if not len(db):
        return {"n_sense": int(n1.sum()), "n_anti": int(n2.sum()), "excess_diff": d0,
                "diff_lo95": np.nan, "diff_hi95": np.nan, "p_boot": np.nan}
    p = 2 * min((db <= 0).mean(), (db >= 0).mean())
    return {"n_sense": int(n1.sum()), "n_anti": int(n2.sum()), "excess_diff": d0,
            "diff_lo95": float(np.percentile(db, 2.5)), "diff_hi95": float(np.percentile(db, 97.5)),
            "p_boot": max(p, 1 / n_boot)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores-dir", required=True)
    ap.add_argument("--chroms", required=True, help="comma-separated, processed in this order")
    ap.add_argument("--gnomad", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--an-frac", type=float, default=0.8)
    ap.add_argument("--neutral", type=float, default=0.001)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-boot-group", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=4_000_000, help="rows read at a time")
    a = ap.parse_args()
    chroms = [c for c in a.chroms.split(",") if c]
    if len(set(chroms)) != len(chroms):
        sys.exit("STOP: a chromosome is listed twice")
    for c in chroms:
        for f in (os.path.join(a.scores_dir, f"{c}.scores.tsv.gz"), os.path.join(a.gnomad, f"{c}.snv_pass.tsv.gz")):
            if not os.path.isfile(f):
                sys.exit(f"STOP: missing {f}")
    os.makedirs(a.out_dir, exist_ok=True)
    t0 = time.time()

    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)
    sp, tagged = canonical_splice_positions(a.gtf, set(chroms))
    print(f"splice-site control positions: {len(sp):,} ({'Ensembl_canonical transcripts' if tagged else 'all protein-coding transcripts'})", flush=True)

    units = Units()
    store, ctrl, loaded = [], {}, []
    max_an = 0
    for c in chroms:
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        isgc = (raw == G_) | (raw == C_) | (raw == G_ + 32) | (raw == C_ + 32)
        cs = np.concatenate([[0], np.cumsum(isgc, dtype=np.int32)])
        del isgc
        R, n_in, mx, n_dedup = load_scores(os.path.join(a.scores_dir, f"{c}.scores.tsv.gz"), c, raw, cs, units,
                                           a.chunk, a.neutral)
        max_an = max(max_an, mx)
        store.append(R)
        ctrl[c] = load_control(os.path.join(a.gnomad, f"{c}.snv_pass.tsv.gz"), c, sp[sp["chrom"] == c].reset_index(drop=True),
                               raw, cs, max(a.chunk, 1_000_000))
        loaded.append({"chrom": c, "rows": n_in, "after_dedup": n_dedup, "valid_context": len(R["AN"]),
                       "max_AN": mx, "splice_ctrl_rows": 0 if ctrl[c] is None else len(ctrl[c])})
        print(f"  {c}: {n_in:,} rows, {n_dedup:,} after duplicates, {len(R['AN']):,} with a valid context; "
              f"max AN {mx:,}; splice-control variants {loaded[-1]['splice_ctrl_rows']:,}  "
              f"[{(time.time() - t0) / 60:.1f} min, peak {mem_gb():.1f} GB]", flush=True)
        del raw, cs
    pd.DataFrame(loaded).to_csv(os.path.join(a.out_dir, "loaded.tsv"), sep="\t", index=False)

    # AN filter (largest AN over all chromosomes), then one table of compact columns
    thr = a.an_frac * max_an
    cols = {}
    for k in store[0]:
        cols[k] = np.concatenate([R[k][R["AN"] >= thr] for R in store])
        for R in store:
            if k != "AN":
                R[k] = None
    del store
    N = len(cols["AN"])
    setc, uid, ctx, cls = cols["setc"], cols["uid"], cols["ctx"].astype(np.int32), cols["cls"]
    cpg = (cols["flags"] & F_CPG) > 0
    single = ((cols["flags"] & F_SINGLE) > 0).astype(np.int8)
    unit_names = np.array(units.names, dtype=object)
    unit_chrom = np.array(units.chrom, dtype=object)
    # local GC: edges at the quintiles of all variants, as 12b (np.quantile on the same float32 values)
    gc_vals = (np.arange(2 * GC_HALF + 2) / (2 * GC_HALF + 1)).astype(np.float32)
    hist = np.bincount(cols["gc"], minlength=len(gc_vals))
    edges = np.quantile(np.repeat(gc_vals, hist), [0.2, 0.4, 0.6, 0.8])
    gc_bin_lut = np.searchsorted(edges, gc_vals).astype(np.int32)
    gcbin = gc_bin_lut[cols["gc"]]
    print(f"\nkept {N:,} variants with AN >= {a.an_frac} x {max_an:,}; at CpG sites {100 * cpg.mean():.1f}%; "
          f"neutral {100 * ((cls & BIT['neutral']) > 0).mean():.1f}%; singletons {100 * single.mean():.1f}%; "
          f"local GC quintile edges: " + ", ".join(f"{e:.3f}" for e in edges) +
          f"  [{(time.time() - t0) / 60:.1f} min, peak {mem_gb():.1f} GB]", flush=True)

    Cc = [ctrl[c] for c in sorted(chroms) if ctrl[c] is not None]       # 12b: control in sorted chromosome order
    C_all = pd.concat(Cc, ignore_index=True) if Cc else None
    if C_all is not None:
        C_all = C_all[C_all["AN"] >= thr].reset_index(drop=True)
        C_all["gcbin"] = gc_bin_lut[C_all["gcc"].to_numpy()]

    def unit_index(rows_uid):
        present = np.unique(rows_uid)
        names = unit_names[present]
        order = np.argsort(names, kind="stable")
        rank = np.full(len(unit_names), -1, dtype=np.int32)
        rank[present[order]] = np.arange(len(present), dtype=np.int32)
        return rank[rows_uid], names[order], unit_chrom[present[order]]

    group_defs = [(c, (c,)) for c in chroms] + [g for g in GROUPS if any(x in chroms for x in g[1])] + \
                 [("training", tuple(c for c in chroms if c not in {x for g in GROUPS for x in g[1]}))]
    group_defs = [g for g in group_defs if g[1]]

    compact, calib, grows = [], [], []
    for si, (name, no_cpg, K) in enumerate(SETTINGS):
        rng = np.random.default_rng(a.seed)       # same stream per setting, as 12b
        rng_g = np.random.default_rng([a.seed, 1000 + si])
        sel = ~cpg if no_cpg else np.ones(N, dtype=bool)
        key_all = ctx * np.int32(K) + (gcbin if K > 1 else np.int32(0))
        print(f"\n######## setting {name}: {int(sel.sum()):,} variants  [{(time.time() - t0) / 60:.1f} min]", flush=True)

        def p_table(keys, sing):
            n = np.bincount(keys, minlength=256 * K).astype(float)
            s = np.bincount(keys, weights=sing, minlength=256 * K)
            overall = s.sum() / max(1.0, n.sum())
            return np.where(n >= 20, s / np.maximum(n, 1), overall)

        sets_here = [s for s in SETS if (sel & (setc == SETS.index(s))).any()]
        rowsets, ptab = {}, {}
        for s in sets_here:
            r = np.flatnonzero(sel & (setc == SETS.index(s)))
            rowsets[s] = r
            nm = (cls[r] & BIT["neutral"]) > 0
            ptab[s] = p_table(key_all[r][nm], single[r][nm])
            # calibration: expectation from half of the units, applied to the neutral variants of the other half
            ur = uid[r][nm]
            pres = np.unique(ur)
            hb = np.zeros(len(unit_names), dtype=bool)
            hb[pres] = [zlib.crc32(str(u).encode()) % 2 == 0 for u in unit_names[pres]]
            half = hb[ur]
            if len(np.unique(ur[half])) < 10 or len(np.unique(ur[~half])) < 10:
                continue
            kk, ss = key_all[r][nm], single[r][nm]
            pa = p_table(kk[half], ss[half])
            ob, ex = ss[~half].mean(), pa[kk[~half]].mean()
            calib.append({"setting": name, "set": s, "n": int((~half).sum()), "observed": ob, "expected": ex,
                          "excess": ob - ex})
            print(f"  calibration {s:22s} excess {ob - ex:+.4f}")

        rows, drows = [], []
        aggs = {}

        def summarise(label, g_uid_rank, n_units, keys, sing, masks, p):
            W = boot_W(n_units, a.n_boot, rng)
            expv = p[keys]
            singf = sing.astype(float)
            agg = {}
            for k, mk in masks:
                n_u = np.bincount(g_uid_rank[mk], minlength=n_units).astype(float)
                s_u = np.bincount(g_uid_rank[mk], weights=singf[mk], minlength=n_units)
                e_u = np.bincount(g_uid_rank[mk], weights=expv[mk], minlength=n_units)
                agg[k] = (n_u, s_u, e_u)
                if n_u.sum() == 0:
                    continue
                rows.append({"set": label, "class": k, **excess_stats(n_u, s_u, e_u, W)})
            return agg, W

        if C_all is not None and "intron_protein_coding" in ptab:
            C = C_all[~C_all["cpg"].to_numpy()].reset_index(drop=True) if no_cpg else C_all
            ckey = C["ctx"].to_numpy().astype(np.int64) * K + (C["gcbin"].to_numpy() if K > 1 else 0)
            cu, cidx = np.unique(C["unit"].to_numpy(dtype=object), return_inverse=True)
            cidx = cidx.ravel()
            agg, W = summarise("positive_control", cidx, len(cu), ckey, C["single"].to_numpy(),
                               [(k, (C["cls"] == k).to_numpy()) for k in CTRL], ptab["intron_protein_coding"])
            cchrom = C.groupby(cidx)["chrom"].first().reindex(range(len(cu))).to_numpy(dtype=object)
            aggs["positive_control"] = (agg, cchrom)
        for s in sets_here:
            r = rowsets[s]
            rank, unames, uchrom = unit_index(uid[r])
            cr = cls[r]
            agg, W = summarise(s, rank, len(unames), key_all[r], single[r],
                               [(k, (cr & BIT[k]) > 0) for k in CLASSES], ptab[s])
            for k1, k2 in PAIRS:
                if k1 not in agg or k2 not in agg or agg[k1][0].sum() == 0 or agg[k2][0].sum() == 0:
                    continue
                drows.append({"set": s, "sense_class": k1, "anti_class": k2, **diff_stats(agg[k1], agg[k2], W, a.n_boot)})
            aggs[s] = (agg, uchrom)
        R = pd.DataFrame(rows)
        Dp = pd.DataFrame(drows)
        R.to_csv(os.path.join(a.out_dir, f"{name}_classes.tsv"), sep="\t", index=False)
        Dp.to_csv(os.path.join(a.out_dir, f"{name}_pairs.tsv"), sep="\t", index=False)

        # per chromosome and per group of chromosomes, each with its own bootstrap over its units
        for gname, gch in group_defs:
            for s, (agg, uchrom) in aggs.items():
                cols_g = np.flatnonzero(np.isin(uchrom, gch))
                if not len(cols_g):
                    continue
                Wg = boot_W(len(cols_g), a.n_boot_group, rng_g)
                sub = {k: tuple(x[cols_g] for x in v) for k, v in agg.items()}
                base = {"setting": name, "group": gname, "chroms": ",".join(gch), "set": s}
                if s == "positive_control":
                    for k in CTRL:
                        if k in sub and sub[k][0].sum() > 0:
                            st = excess_stats(*sub[k], Wg)
                            grows.append({**base, "quantity": f"excess {k}", "n": st["n"], "units": st["units"],
                                          "value": st["excess"], "lo95": st["excess_lo95"], "hi95": st["excess_hi95"],
                                          "p_boot": np.nan})
                    continue
                for k1, k2 in PAIRS:
                    if k1 not in sub or k2 not in sub or sub[k1][0].sum() == 0 or sub[k2][0].sum() == 0:
                        continue
                    st = diff_stats(sub[k1], sub[k2], Wg, a.n_boot_group)
                    grows.append({**base, "quantity": f"diff {k1} - anti", "n": st["n_sense"] + st["n_anti"],
                                  "units": int(((sub[k1][0] + sub[k2][0]) > 0).sum()), "value": st["excess_diff"],
                                  "lo95": st["diff_lo95"], "hi95": st["diff_hi95"], "p_boot": st["p_boot"]})

        def get(df, s, col, k, field):
            x = df[(df["set"] == s) & (df[col] == k)] if len(df) else df
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
    pd.DataFrame(calib).to_csv(os.path.join(a.out_dir, "calibration.tsv"), sep="\t", index=False)
    Gm = pd.DataFrame(grows)
    Gm.to_csv(os.path.join(a.out_dir, "groups.tsv"), sep="\t", index=False)

    print("\n== compact comparison across settings ==")
    for (s, q), g in Cm.groupby(["set", "quantity"], sort=False):
        line = "; ".join(f"{r.setting} {r.value:+.4f} [{r.lo95:+.4f}, {r.hi95:+.4f}]" + (f" p={r.p_boot:.3g}" if r.p_boot == r.p_boot else "")
                         for r in g.itertuples())
        print(f"  {s:22s} {q:38s} {line}")

    if len(Gm):
        print("\n== by chromosome / group (own bootstrap per row; settings base and noCpG_gc5) ==")
        show = [("positive_control", "excess splice_canonical"),
                ("intron_protein_coding", "diff creates_sense_0.0439 - anti"),
                ("intron_protein_coding", "diff gain_sense_strong - anti"),
                ("intergenic", "diff creates_sense_0.0439 - anti")]
        for s, q in show:
            print(f"  {s} | {q}")
            for gname, _ in group_defs:
                parts = []
                for st in ("base", "noCpG_gc5"):
                    x = Gm[(Gm["setting"] == st) & (Gm["group"] == gname) & (Gm["set"] == s) & (Gm["quantity"] == q)]
                    if len(x):
                        r = x.iloc[0]
                        parts.append(f"{st} {r.value:+.4f} [{r.lo95:+.4f}, {r.hi95:+.4f}] n={int(r.n):,}")
                if parts:
                    print(f"    {gname:14s} " + "; ".join(parts))
    print(f"\nwritten to {a.out_dir}: <setting>_classes.tsv, <setting>_pairs.tsv, compact.tsv, groups.tsv, "
          f"calibration.tsv, loaded.tsv  [{(time.time() - t0) / 60:.1f} min, peak {mem_gb():.1f} GB]")


if __name__ == "__main__":
    main()
