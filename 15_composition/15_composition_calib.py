#!/usr/bin/env python3
"""Step 15: is the part of the sense weakness that Puffin-D does not explain carried by
strand composition that the model under-reads?

The sense strand of genes is T- and G-rich (transcription-coupled mutation in the germ
line), whether or not the gene is active in a given cell line. Step 10 showed that the
sense deficit of potential start sites is carried by composition, and step 14 found the
observed sense weakness as strong in silent as in active hosts. Puffin-D's outputs are
compressed, so it may capture only part of what composition does to the balance of two
divergent arms. This step measures that directly where composition varies at random:

  1. Intergenic pairs (orientation arbitrary: "sense" = plus arm). For each pair, the
     composition of a window around it is read on the sense strand as odd (orientation-
     sensitive) features: skews of T vs A and G vs C, and of the six non-palindromic
     dinucleotide pairs (TT/AA, TG/CA, TC/GA, GT/AC, GG/CC, CT/AG) - the same
     composition that the step-10 shuffle preserved. Observed DI is regressed on these
     (and on Puffin-D's pred_DI), without intercept: taking every control pair in both
     orientations, as 07b, gives exactly this fit, and orientation-free properties such
     as GC level drop out by symmetry.
  2. The fitted relation is applied to intronic pairs, read on the host's sense strand:
     expected DI shift = beta . mean(features); residual = observed - expected.
Models: pred (= 07b), mono, mono+di, pred+mono, pred+mono+di. Windows: from the outer arm
minus W to the outer arm plus W, W = --windows.

Reading: if adding composition moves the residual to ~0, the unexplained part is mostly
composition the model under-reads; if the residual stays, look elsewhere. Out-of-fold
R^2 in the controls (no-intercept, around 0) shows whether composition predicts arm
balance at all where there is no gene. A caveat to keep in mind: transcription in the
germ line itself skews composition, so in intergenic pairs a habitually stronger arm may
have made its own strand T/G-rich (the reverse direction); the sign of the fitted
effect tells which way the controls lean.

CIs: bootstrap over host genes for intronic pairs and over pairs for controls, the fit
redone in every resample. Subsets: all chromosomes, and the held-out chromosomes (fit
and application both within them, as 07b).

Usage:
  python3 15_composition_calib.py --pairs-dir ../07_puffin_pairs/results --lines K562,HCT116 \
      --genome hg38.fa --out-prefix results/15_comp [--windows 100,300,1000] [--n-boot 1000]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

INTRON = "intron genu"
ACTIVE = 5.0
OBS, PRED = "DI", "pred_DI"
LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i
A_, C_, G_, T_ = 0, 1, 2, 3
DI_PAIRS = (("TT", "AA"), ("TG", "CA"), ("TC", "GA"), ("GT", "AC"), ("GG", "CC"), ("CT", "AG"))
MONO = ["skew_TA", "skew_GC"]
DIN = [f"skew_{x}_{y}" for x, y in DI_PAIRS]
MODELS = {"pred": [PRED], "mono": MONO, "mono+di": MONO + DIN, "pred+mono": [PRED] + MONO,
          "pred+mono+di": [PRED] + MONO + DIN}


def dcode(s):
    return "ACGT".index(s[0]) * 4 + "ACGT".index(s[1])


def skew(a, b):
    t = a + b
    return np.where(t > 0, (a - b) / np.maximum(t, 1), 0.0)


def window_features(seq_codes, lo, hi):
    """Plus-strand odd composition features of seq[lo:hi] (codes 0..3, -1 = N)."""
    x = seq_codes[lo:hi]
    v = x[x >= 0]
    m = np.bincount(v, minlength=4).astype(float)
    ok = (x[:-1] >= 0) & (x[1:] >= 0)
    dn = np.bincount((x[:-1][ok].astype(np.int16) * 4 + x[1:][ok]), minlength=16).astype(float)
    f = [skew(m[T_], m[A_]), skew(m[G_], m[C_])]
    f += [skew(dn[dcode(p)], dn[dcode(q)]) for p, q in DI_PAIRS]
    return np.array(f, dtype=float), len(v) / max(1, len(x))


def add_features(d, fa, windows):
    """Features for every pair and window, read on the pair's sense strand."""
    sign = np.where(d["strand_sense"].astype(str) == "-", -1.0, 1.0)
    lo_arm = np.minimum(d["pos_plus"], d["pos_minus"]).to_numpy(np.int64)
    hi_arm = np.maximum(d["pos_plus"], d["pos_minus"]).to_numpy(np.int64)
    out = {W: np.full((len(d), len(MONO) + len(DIN)), np.nan) for W in windows}
    acgt = {W: np.full(len(d), np.nan) for W in windows}
    for c, idx in d.groupby("chrom").indices.items():
        codes = LUT[np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)]
        L = len(codes)
        for i in idx:
            for W in windows:
                lo, hi = max(0, int(lo_arm[i]) - W), min(L, int(hi_arm[i]) + W)
                f, frac = window_features(codes, lo, hi)
                out[W][i] = sign[i] * f
                acgt[W][i] = frac
    return out, acgt


def ols(XtX, Xty):
    try:
        return np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(XtX, Xty, rcond=None)[0]


def boot_counts(n, n_boot, rng):
    idx = rng.integers(0, n, size=(n_boot, n))
    W = np.zeros((n_boot, n), dtype=np.float64)
    for b in range(n_boot):
        W[b] = np.bincount(idx[b], minlength=n)
    return W


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) else (np.nan, np.nan)


def cv_r2(X, y, folds, rng):
    """Out-of-fold R^2 of a no-intercept fit, relative to 0 (the symmetric mean)."""
    perm = rng.permutation(len(y))
    parts = np.array_split(perm, folds)
    ss = 0.0
    for k in range(folds):
        te = parts[k]
        tr = np.concatenate([parts[j] for j in range(folds) if j != k])
        b = ols(X[tr].T @ X[tr], X[tr].T @ y[tr])
        ss += float(np.sum((y[te] - X[te] @ b) ** 2))
    return 1 - ss / float(np.sum(y * y))


def analyse(ctrl, cases, units, cols, n_boot, rng, folds):
    """Fit on controls (no intercept), apply to each case group; bootstrap both."""
    Xc, yc = ctrl[cols].to_numpy(float), ctrl[OBS].to_numpy(float)
    beta = ols(Xc.T @ Xc, Xc.T @ yc)
    r2 = cv_r2(Xc, yc, folds, rng)
    Wc = boot_counts(len(yc), n_boot, rng)
    P = Xc[:, :, None] * Xc[:, None, :]
    XtX_b = (Wc @ P.reshape(len(yc), -1)).reshape(n_boot, len(cols), len(cols))
    Xty_b = Wc @ (Xc * yc[:, None])
    beta_b = np.array([ols(XtX_b[b], Xty_b[b]) for b in range(n_boot)])
    uu, ui = np.unique(units, return_inverse=True)
    U = len(uu)
    Wg = boot_counts(U, n_boot, rng)
    res = []
    for gname, g in cases.items():
        if len(g) < 20:
            continue
        gi = ui[g.index.to_numpy()]
        X, y = g[cols].to_numpy(float), g[OBS].to_numpy(float)
        n_u = np.bincount(gi, minlength=U).astype(float)
        sx = np.stack([np.bincount(gi, weights=X[:, j], minlength=U) for j in range(len(cols))], axis=1)
        sy = np.bincount(gi, weights=y, minlength=U)
        nb = Wg @ n_u
        with np.errstate(invalid="ignore", divide="ignore"):
            mx_b = (Wg @ sx) / nb[:, None]
            my_b = (Wg @ sy) / nb
        exp_b = np.sum(beta_b * mx_b, axis=1)
        mx, my = X.mean(0), y.mean()
        e = float(beta @ mx)
        res.append({"group": gname, "n_case": len(g), "n_ctrl": len(yc), "cv_r2_ctrl": r2, "obs_mean_DI": my,
                    "expected_DI": e, "expected_lo95": ci(exp_b)[0], "expected_hi95": ci(exp_b)[1],
                    "residual_DI": my - e, "residual_lo95": ci(my_b - exp_b)[0], "residual_hi95": ci(my_b - exp_b)[1],
                    "share": e / my if my != 0 else np.nan, "share_lo95": ci(exp_b / my_b)[0],
                    "share_hi95": ci(exp_b / my_b)[1]})
    coefs = [{"feature": c, "beta": float(beta[j]), "lo95": ci(beta_b[:, j])[0], "hi95": ci(beta_b[:, j])[1]}
             for j, c in enumerate(cols)]
    return res, coefs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-dir", required=True, help="directory with puffin_rc_pairs_<line>.tsv (step 07c)")
    ap.add_argument("--lines", default="K562,HCT116")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--windows", default="100,300,1000")
    ap.add_argument("--heldout", default="chr8,chr9")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    lines = [x for x in a.lines.split(",") if x]
    windows = [int(x) for x in a.windows.split(",") if x]
    heldout = set(a.heldout.split(","))
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    rng = np.random.default_rng(a.seed)
    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)
    pd.set_option("display.width", 250)
    fmt = lambda v: f"{v:.3g}"

    rows, crows = [], []
    for L in lines:
        f = os.path.join(a.pairs_dir, f"puffin_rc_pairs_{L}.tsv")
        if not os.path.isfile(f):
            sys.exit(f"STOP: missing {f}")
        d = pd.read_csv(f, sep="\t", keep_default_na=False, na_values=["", "NA", "nan", "NaN"],
                        dtype={"host": str, "klasa": str, "strand_sense": str})
        fin = [c for c in [OBS, "gc_local", PRED, "pred_DI_point", "pred_DI_wide", "pred_DI_sum"] if c in d.columns]
        d = d[(d["pred_status"] == "ok") & np.isfinite(d[fin].to_numpy(float)).all(axis=1)].reset_index(drop=True)
        feats, acgt = add_features(d, fa, windows)
        intr = (d["klasa"] == INTRON).to_numpy()
        hs = pd.to_numeric(d["host_signal"], errors="coerce").to_numpy(float)
        units = np.where(intr, d["host"].fillna("").astype(str), "")
        units = np.where(units == "", "pair" + d.index.astype(str), units)
        print(f"\n######## {L}: {intr.sum():,} intronic, {(~intr).sum():,} intergenic pairs")
        for W in windows:
            F = feats[W]
            print(f"  window +/-{W} nt: mean sense-strand skew, intronic vs intergenic: " +
                  "; ".join(f"{nm} {np.nanmean(F[intr, j]):+.4f} vs {np.nanmean(F[~intr, j]):+.4f}"
                            for j, nm in enumerate(MONO)) + f"; min ACGT fraction {np.nanmin(acgt[W]):.2f}")
        for subset in ("all", "heldout"):
            sub = np.ones(len(d), bool) if subset == "all" else d["chrom"].isin(heldout).to_numpy()
            for W in windows:
                D = d[[OBS, PRED, "chrom"]].copy()
                for j, nm in enumerate(MONO + DIN):
                    D[nm] = feats[W][:, j]
                ok = sub & np.isfinite(D[MONO + DIN].to_numpy()).all(axis=1)
                ctrl = D[ok & ~intr]
                cases = {"intron_all": D[ok & intr]}
                if subset == "all":
                    cases["intron_active"] = D[ok & intr & (hs >= ACTIVE)]
                    cases["intron_silent"] = D[ok & intr & (hs < ACTIVE)]
                if len(ctrl) < 50:
                    continue
                for mname, cols in MODELS.items():
                    if mname == "pred" and W != windows[0]:
                        continue                     # does not depend on the window
                    r, cf = analyse(ctrl, cases, units, cols, a.n_boot, rng, a.folds)
                    wl = "-" if mname == "pred" else str(W)
                    rows += [{"line": L, "subset": subset, "window": wl, "model": mname, **x} for x in r]
                    crows += [{"line": L, "subset": subset, "window": wl, "model": mname, **x} for x in cf]
    R = pd.DataFrame(rows)
    Cf = pd.DataFrame(crows)
    R.to_csv(a.out_prefix + "_models.tsv", sep="\t", index=False)
    Cf.to_csv(a.out_prefix + "_coefs.tsv", sep="\t", index=False)

    show = ["window", "model", "n_case", "cv_r2_ctrl", "obs_mean_DI", "expected_DI", "residual_DI", "residual_lo95",
            "residual_hi95", "share", "share_lo95", "share_hi95"]
    for L in lines:
        for subset, group in (("all", "intron_all"), ("all", "intron_silent"), ("all", "intron_active"),
                              ("heldout", "intron_all")):
            x = R[(R.line == L) & (R.subset == subset) & (R.group == group)]
            if len(x):
                print(f"\n== {L} | {subset} chromosomes | {group}: observed vs expected sense-minus-antisense DI ==")
                print(x[show].to_string(index=False, float_format=fmt))
        x = Cf[(Cf.line == L) & (Cf.subset == "all") & (Cf.model.isin(["pred+mono", "mono"]))]
        if len(x):
            print(f"\n== {L} | coefficients fitted on intergenic pairs (all chromosomes) ==")
            print(x[["window", "model", "feature", "beta", "lo95", "hi95"]].to_string(index=False, float_format=fmt))
    print("\nDI = sense - antisense (< 0: sense weaker). expected = fit on intergenic pairs applied to intronic ones; "
          "share = expected / observed; cv_r2_ctrl = out-of-fold R^2 in intergenic pairs.")
    print(f"CIs: bootstrap ({a.n_boot}) over host genes (intronic) and pairs (intergenic), fit redone per resample.")
    print(f"written: {a.out_prefix}_models.tsv, {a.out_prefix}_coefs.tsv")


if __name__ == "__main__":
    main()
