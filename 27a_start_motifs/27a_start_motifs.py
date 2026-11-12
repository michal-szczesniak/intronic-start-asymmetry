#!/usr/bin/env python3
"""Step 27a: sequence around the two starts of each pair - logos and pre-specified core-promoter elements.

Question. The sense start of intronic pairs is weaker than its antisense partner, more so than in intergenic
pairs of the same read depth, and ProCapNet predicts 30-42% of this from the element's sequence. Do known
core-promoter elements differ between the sense and the antisense start of intronic pairs, beyond what is seen
between the two (arbitrarily labelled) starts of intergenic pairs?

Pairs. The depth-matched pairs of step 23b (23b_all_matched_pairs.tsv): every intronic pair in an active host
with its intergenic partner of nearly the same depth (K562 1,610, HCT116 2,489). Start positions and strands
come from the iSTART arm tables (arms_<line>.tsv: pos_sense, pos_antisense, strand_sense, strand_antisense;
joined on chrom + mid; positions are 0-based genome indices, as used in steps 07, 16b and 20a). For intergenic
pairs the "sense" arm is the plus-strand arm (arbitrary), so their two arms are called plus and minus here.

Each start is read on its own strand, 5'->3' in the direction of transcription, position 0 = the arm position
(first transcribed base), from -F to +F (--flank). Groups: IS / IA = sense / antisense start of intronic
pairs; GP / GM = plus / minus start of the matched intergenic pairs.

Pre-specified features (motif start positions relative to position 0; exact matches, IUPAC):
  TATA   TATAWAWR, start in [-36, -24]           Inr   BBCABW with its A in [-3, +3]
  DPE    RGWYV, start in [+24, +32]              YY1f  AANATGGC anywhere in [-60, +52] (transcribed strand)
  YY1r   GCCATNTT anywhere in [-60, +52]         GC    G+C fraction in [-50, +50]     CpG  CG count in [-50, +50]
  sensitivity: TATA1 = TATAWAWR with at most one mismatch, same positions.
Composition-expected value of each motif feature: its mean over --n-shuf dinucleotide shuffles (step 10's
Altschul-Erickson shuffle) of the window [-60, +60] of the same start; adj = observed - expected.

Per pair: d = F(sense start) - F(antisense start) for intronic pairs, c = F(plus) - F(minus) for the matched
intergenic partner (expected 0). Reported: mean F per group; D_intron = mean d; D_inter = mean c; Diff =
D_intron - D_inter, each with a 95% CI from --n-boot resamples of host genes (each intronic pair travels with
its partner) and a two-sided bootstrap p. PRIMARY: Diff of the composition-adjusted motif features (TATA, Inr,
DPE, YY1f, YY1r), Holm-corrected over these five within each cell line. Positive control: Spearman correlation
between c and the partner's DI (plus/minus asymmetry) in intergenic pairs; real core-promoter features should
correlate positively. Secondary: the same correlation between d and DI in intronic pairs.

Also written: base counts per position and group (for logos), and the frequency of a pyrimidine-purine (YR)
step at positions -1/0 shifted by -10..+10 nt (a check of how precisely the arm positions mark the start).

Usage:
  python3 27a_start_motifs.py --arms K562=arms_K562.tsv --arms HCT116=arms_HCT116.tsv \\
      --matched 23b_all_matched_pairs.tsv --genome hg38.fa --step10 10_shuffle_test.py --out-prefix results/27a \\
      [--flank 40] [--n-shuf 20] [--n-boot 2000] [--seed 1]
  or in two steps (only the first needs the genome):
  python3 27a_start_motifs.py --arms ... --matched ... --genome hg38.fa --extract-windows 27a_windows.fa.gz
  python3 27a_start_motifs.py --arms ... --matched ... --windows 27a_windows.fa.gz --step10 ... --out-prefix ...
"""
import argparse
import importlib.util
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ARM_COLS = ["chrom", "mid", "klasa", "host", "pos_sense", "pos_antisense", "strand_sense", "strand_antisense",
            "status"]
MATCH_COLS = ["line", "chrom", "mid", "host", "DI", "ctrl_chrom", "ctrl_mid", "ctrl_DI"]
INTRON = "intron genu"
W = 60                       # half-width of the scanned window
COMP = str.maketrans("ACGTN", "TGCAN")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "W": "[AT]", "R": "[AG]", "Y": "[CT]", "B": "[CGT]",
         "V": "[ACG]", "N": "[ACGT]"}
MOTIFS = {  # name: (IUPAC consensus, first allowed start, last allowed start), positions relative to 0
    "TATA": ("TATAWAWR", -36, -24),
    "Inr": ("BBCABW", -6, 0),
    "DPE": ("RGWYV", 24, 32),
    "YY1f": ("AANATGGC", -60, 52),
    "YY1r": ("GCCATNTT", -60, 52),
}
PRIMARY = ["TATA", "Inr", "DPE", "YY1f", "YY1r"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def rx(cons):
    return re.compile("(?=(" + "".join(IUPAC[c] for c in cons) + "))")


RX = {k: rx(v[0]) for k, v in MOTIFS.items()}
SETS = {k: [set(IUPAC[c].strip("[]")) for c in v[0]] for k, v in MOTIFS.items()}


def has_motif(seq, name, mism=0):
    """seq: oriented window of length 2W+1 (index W = position 0). True if the motif starts in its range."""
    cons, lo, hi = MOTIFS[name]
    a, b = W + lo, W + hi
    if mism == 0:
        return any(a <= m.start() <= b for m in RX[name].finditer(seq))
    sets, k = SETS[name], len(SETS[name])
    for s0 in range(max(a, 0), min(b, len(seq) - k) + 1):
        if sum(seq[s0 + j] not in sets[j] for j in range(k)) <= mism:
            return True
    return False


def gc_cpg(seq):
    s = seq[W - 50:W + 51]
    acgt = sum(s.count(x) for x in "ACGT")
    return (s.count("G") + s.count("C")) / max(acgt, 1), float(s.count("CG"))


def window(fa, chrom, p, strand):
    """Oriented window of 2W+1 nt with the arm position at index W; None near a chromosome end.
    fa: a pyfaidx Fasta of the genome, or a dict {"chrom:lo-hi": plus-strand sequence} (--windows)."""
    lo, hi = p - W, p + W + 1
    if isinstance(fa, dict):
        s = fa.get(f"{chrom}:{lo}-{hi}")
        if s is None:
            return None
        s = s.upper()
    else:
        if chrom not in fa:
            sys.exit(f"STOP: {chrom} not in the genome")
        if lo < 0 or hi > len(fa[chrom]):
            return None
        s = str(fa[chrom][lo:hi]).upper()
    return s if strand == "+" else s.translate(COMP)[::-1]


def read_windows(path):
    """FASTA written by --extract-windows: >chrom:lo-hi (0-based, half-open), plus-strand sequence."""
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    out, name, buf = {}, None, []
    with op(path, "rt") as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = ln[1:], []
            elif ln:
                buf.append(ln)
    if name:
        out[name] = "".join(buf)
    return out


def starts_of(arms_specs, matched):
    """Yield (line, group, chrom, pos, strand) for the four starts of every matched pair (checks as in main)."""
    m = pd.read_csv(matched, sep="\t", dtype={"chrom": str, "ctrl_chrom": str, "host": str})
    for spec in arms_specs:
        line, path = spec.split("=", 1)
        arms = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], dtype={"chrom": str, "host": str})
        A = arms.set_index(["chrom", "mid"])
        mm = m[m["line"] == line]
        for keys, grp in ((list(zip(mm["chrom"], mm["mid"])), ("IS", "IA")),
                          (list(zip(mm["ctrl_chrom"], mm["ctrl_mid"])), ("GP", "GM"))):
            sub = A.loc[keys].reset_index()
            for c, ps, ss, pa, sa in zip(sub["chrom"], sub["pos_sense"], sub["strand_sense"], sub["pos_antisense"],
                                         sub["strand_antisense"]):
                yield line, grp[0], str(c), int(ps), str(ss)
                yield line, grp[1], str(c), int(pa), str(sa)


def features(seq, rng, s10, n_shuf):
    out = {k: float(has_motif(seq, k)) for k in MOTIFS}
    out["TATA1"] = float(has_motif(seq, "TATA", 1))
    out["GC"], out["CpG"] = gc_cpg(seq)
    exp = dict.fromkeys(list(MOTIFS) + ["TATA1"], 0.0)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    for _ in range(n_shuf):
        sh = s10.dinuc_shuffle(arr, rng).tobytes().decode()
        for k in MOTIFS:
            exp[k] += has_motif(sh, k)
        exp["TATA1"] += has_motif(sh, "TATA", 1)
    for k in exp:
        out[k + "_exp"] = exp[k] / n_shuf
    return out


def boot(units, fn, n_boot, rng):
    u, inv = np.unique(units, return_inverse=True)
    idx = [np.flatnonzero(inv == i) for i in range(len(u))]
    vals = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, len(u), len(u))
        vals[b] = fn(np.concatenate([idx[i] for i in pick]))
    return vals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", action="append", required=True, help="LINE=arms_<line>.tsv")
    ap.add_argument("--matched", required=True)
    ap.add_argument("--genome", default="", help="indexed genome FASTA (pyfaidx)")
    ap.add_argument("--windows", default="", help="instead of --genome: FASTA written by --extract-windows")
    ap.add_argument("--extract-windows", default="",
                    help="only write the +/-60-nt plus-strand windows of all starts to this FASTA(.gz) and stop")
    ap.add_argument("--step10", default="")
    ap.add_argument("--out-prefix", default="")
    ap.add_argument("--flank", type=int, default=40)
    ap.add_argument("--n-shuf", type=int, default=20)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.flank > W:
        sys.exit(f"STOP: --flank must be <= {W}")
    if a.extract_windows:
        import gzip
        need = {}
        for line, grp, c, p, strand in starts_of(a.arms, a.matched):
            need.setdefault(c, set()).add((p - W, p + W + 1))
        n_need = sum(len(v) for v in need.values())
        seen, n_edge = 0, 0
        op = gzip.open if a.extract_windows.endswith(".gz") else open
        with op(a.extract_windows, "wt") as fo:
            def emit(c, get, length):
                nonlocal seen, n_edge
                for lo, hi in sorted(need.get(c, ())):
                    if lo < 0 or hi > length:
                        n_edge += 1
                        continue
                    fo.write(f">{c}:{lo}-{hi}\n{get(lo, hi)}\n")
                    seen += 1
            if a.genome.endswith(".gz"):
                # plain-gzip FASTA (e.g. UCSC hg38.fa.gz): one pass, one chromosome in memory at a time
                found, name, buf = set(), None, bytearray()
                with gzip.open(a.genome, "rb") as f:
                    for ln in f:
                        if ln.startswith(b">"):
                            if name in need:
                                seq = buf.decode()
                                emit(name, lambda lo, hi: seq[lo:hi], len(seq))
                                found.add(name)
                            name, buf = ln[1:].split()[0].decode(), bytearray()
                        elif name in need:
                            buf.extend(ln.rstrip())
                    if name in need:
                        seq = buf.decode()
                        emit(name, lambda lo, hi: seq[lo:hi], len(seq))
                        found.add(name)
                lost = sorted(set(need) - found)
                if lost:
                    sys.exit(f"STOP: chromosomes not in the genome: {lost}")
            else:
                from pyfaidx import Fasta
                fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)
                for c in need:
                    if c not in fa:
                        sys.exit(f"STOP: {c} not in the genome")
                    emit(c, lambda lo, hi, c=c: str(fa[c][lo:hi]), len(fa[c]))
        print(f"written: {a.extract_windows}: {seen:,} of {n_need:,} windows of {2 * W + 1} nt "
              f"({n_edge} skipped at chromosome ends)")
        return
    if not (a.step10 and a.out_prefix and (a.genome or a.windows)):
        sys.exit("STOP: --step10, --out-prefix and --genome or --windows are required")
    s10 = load_module("s10", a.step10)
    if a.windows:
        fa = read_windows(a.windows)
        print(f"windows: {len(fa):,} read from {a.windows}")
    else:
        from pyfaidx import Fasta
        fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)
    m = pd.read_csv(a.matched, sep="\t", dtype={"chrom": str, "ctrl_chrom": str, "host": str})
    miss = [c for c in MATCH_COLS if c not in m.columns]
    if miss:
        sys.exit(f"STOP: {a.matched} lacks columns {miss}")
    rng = np.random.default_rng(a.seed)
    starts, pfm, yr, tests = [], [], [], []
    for spec in a.arms:
        line, path = spec.split("=", 1)
        arms = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], dtype={"chrom": str, "host": str})
        miss = [c for c in ARM_COLS if c not in arms.columns]
        if miss:
            sys.exit(f"STOP: {path} lacks columns {miss}")
        if arms.duplicated(["chrom", "mid"]).any():
            sys.exit(f"STOP: {path}: chrom + mid is not unique")
        mm = m[m["line"] == line].reset_index(drop=True)
        if not len(mm):
            sys.exit(f"STOP: no matched pairs for {line}")
        A = arms.set_index(["chrom", "mid"])
        ki = list(zip(mm["chrom"], mm["mid"]))
        kc = list(zip(mm["ctrl_chrom"], mm["ctrl_mid"]))
        for k in ki + kc:
            if k not in A.index:
                sys.exit(f"STOP: {line}: pair {k} of the matched table is not in {path}")
        ai, ac = A.loc[ki].reset_index(), A.loc[kc].reset_index()
        if not ((ai["klasa"] == INTRON).all() and (ac["klasa"] != INTRON).all()):
            sys.exit(f"STOP: {line}: classes of matched pairs are not intronic / intergenic")
        if not (ac["strand_sense"] == "+").all():
            sys.exit(f"STOP: {line}: an intergenic 'sense' arm is not on the plus strand")
        if not ((ai["strand_sense"] != ai["strand_antisense"]).all() and (ac["strand_sense"] != ac["strand_antisense"]).all()):
            sys.exit(f"STOP: {line}: the two arms of a pair are on the same strand")
        rows = []
        for i in range(len(mm)):
            rec = {"line": line, "pair": i, "host": str(mm["host"].iat[i]), "DI": float(mm["DI"].iat[i]),
                   "ctrl_DI": float(mm["ctrl_DI"].iat[i])}
            ok = True
            for grp, src, pcol, scol in (("IS", ai, "pos_sense", "strand_sense"), ("IA", ai, "pos_antisense", "strand_antisense"),
                                          ("GP", ac, "pos_sense", "strand_sense"), ("GM", ac, "pos_antisense", "strand_antisense")):
                seq = window(fa, str(src["chrom"].iat[i]), int(src[pcol].iat[i]), str(src[scol].iat[i]))
                if seq is None or seq.count("N") > 10:
                    ok = False
                    break
                rec[grp + "_seq"] = seq
            if ok:
                rows.append(rec)
        d = pd.DataFrame(rows)
        print(f"\n######## {line}: {len(mm):,} matched pairs; {len(d):,} with all four windows readable "
              f"(+/-{W} nt, <= 10 N)")
        # features
        for g in ("IS", "IA", "GP", "GM"):
            F = pd.DataFrame([features(s, rng, s10, a.n_shuf) for s in d[g + "_seq"]])
            for col in F.columns:
                d[f"{g}_{col}"] = F[col].to_numpy()
            # base counts for logos and the YR check
            arr = np.array([list(s) for s in d[g + "_seq"]])
            for pos in range(-a.flank, a.flank + 1):
                col = arr[:, W + pos]
                pfm.append({"line": line, "group": g, "pos": pos, **{b: int((col == b).sum()) for b in "ACGT"},
                            "n": len(col)})
            for off in range(-10, 11):
                y = np.isin(arr[:, W + off - 1], ["C", "T"]) & np.isin(arr[:, W + off], ["A", "G"])
                yr.append({"line": line, "group": g, "offset": off, "frac_YR": float(y.mean()), "n": len(y)})
        for g in ("IS", "IA", "GP", "GM"):
            starts.append(d[[c for c in d.columns if c.startswith(g + "_") and not c.endswith("_seq")] + ["pair", "host"]]
                          .rename(columns=lambda c: c[len(g) + 1:] if c.startswith(g + "_") else c)
                          .assign(line=line, group=g))
        units = d["host"].astype(str).to_numpy()
        feats = PRIMARY + ["TATA1", "GC", "CpG"]
        res = []
        for f in feats:
            for meas in (("obs", "adj") if f not in ("GC", "CpG") else ("obs",)):
                def val(g):
                    v = d[f"{g}_{f}"].to_numpy(float)
                    return v - d[f"{g}_{f}_exp"].to_numpy(float) if meas == "adj" else v
                dI = val("IS") - val("IA")
                dC = val("GP") - val("GM")
                est = {"D_intron": dI.mean(), "D_inter": dC.mean(), "Diff": dI.mean() - dC.mean()}
                bD = boot(units, lambda k: dI[k].mean(), a.n_boot, rng)
                bF = boot(units, lambda k: dI[k].mean() - dC[k].mean(), a.n_boot, rng)
                bC = boot(units, lambda k: dC[k].mean(), a.n_boot, rng)
                r = {"line": line, "feature": f, "measure": meas, "n_pairs": len(d), "n_host_genes": len(np.unique(units)),
                     **{f"mean_{g}": float(val(g).mean()) for g in ("IS", "IA", "GP", "GM")}, **est}
                for nm, b in (("D_intron", bD), ("D_inter", bC), ("Diff", bF)):
                    r[nm + "_lo95"], r[nm + "_hi95"] = np.percentile(b, [2.5, 97.5])
                    r[nm + "_p"] = min(1.0, 2 * min((b <= 0).mean(), (b >= 0).mean()))
                rho_c = spearmanr(dC, d["ctrl_DI"]).statistic if np.ptp(dC) > 0 else np.nan
                rho_i = spearmanr(dI, d["DI"]).statistic if np.ptp(dI) > 0 else np.nan
                r["rho_inter_c_vs_DI"], r["rho_intron_d_vs_DI"] = rho_c, rho_i
                res.append(r)
        R = pd.DataFrame(res)
        prim = (R["feature"].isin(PRIMARY)) & (R["measure"] == "adj")
        p = R.loc[prim, "Diff_p"].to_numpy()
        order = np.argsort(p)
        holm = np.empty_like(p)
        run = 0.0
        for rank, j in enumerate(order):
            run = max(run, min(1.0, (len(p) - rank) * p[j]))
            holm[j] = run
        R["Diff_p_holm"] = np.nan
        R.loc[prim, "Diff_p_holm"] = holm
        R["primary"] = prim
        tests.append(R)
        show = R[["feature", "measure", "mean_IS", "mean_IA", "mean_GP", "mean_GM", "D_intron", "D_intron_lo95",
                  "D_intron_hi95", "Diff", "Diff_lo95", "Diff_hi95", "Diff_p", "Diff_p_holm", "rho_inter_c_vs_DI",
                  "rho_intron_d_vs_DI"]]
        print(show.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
        Y = pd.DataFrame(yr)
        Y = Y[Y["line"] == line].pivot(index="offset", columns="group", values="frac_YR")
        print("YR step frequency by shift of the arm position (peak at 0 = arms mark the first transcribed base):")
        print(Y.to_string(float_format=lambda x: f"{x:.3f}"))
    pd.DataFrame(pfm).to_csv(a.out_prefix + "_pfm.tsv", sep="\t", index=False)
    pd.DataFrame(yr).to_csv(a.out_prefix + "_yr_offset.tsv", sep="\t", index=False, float_format="%.5g")
    pd.concat(starts).to_csv(a.out_prefix + "_starts.tsv.gz", sep="\t", index=False, float_format="%.5g")
    pd.concat(tests).to_csv(a.out_prefix + "_tests.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples of host genes (seed {a.seed}); p two-sided bootstrap; "
          f"expected values from {a.n_shuf} dinucleotide shuffles of +/-{W} nt per start.")
    print(f"written: {a.out_prefix}_pfm.tsv, _yr_offset.tsv, _starts.tsv.gz, _tests.tsv")


if __name__ == "__main__":
    main()
