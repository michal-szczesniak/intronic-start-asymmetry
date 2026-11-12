#!/usr/bin/env python3
"""Step 07: does Puffin-D reproduce the iSTART strand asymmetry from reference sequence alone?

This is the go/no-go test of PURGE. The observed asymmetry (the arm on the host
strand is weaker than the arm against it) could come from chromatin or from the
DNA sequence itself. Puffin-D predicts PRO-cap initiation from sequence only, with
no cell-type information, so a deficit in its predictions can only come from
sequence.

For each divergent pair from iSTART, one 100 kb window centred on the pair's
midpoint is scored once; both arms lie inside it (arms are at most ~800 nt apart).
Predicted initiation is read at each arm on its own strand and summarised as

    pred_DI = (pred_sense - pred_antisense) / ln(2)

which equals log2((count_sense + 1) / (count_antisense + 1)) if the model emits
ln(count + 1), as the authors state for Puffin. It is therefore directly
comparable with the observed DI from iSTART, which is log2((sens + 1)/(anty + 1)).
The raw difference is also written out, so the summary does not depend on that
assumption.

Two assumptions are CHECKED, not taken on trust, and reported before any result:

  1. Coordinate alignment. For a sample of arms, the position of the predicted
     maximum within +/-500 nt of the arm is reported. A median offset near 0 means
     output position i corresponds to genomic position window_start + i.
  2. Track identity. Tracks 4 and 5 are said to be PRO-cap on the forward and
     reverse strand of the input sequence. Since the input is always the genomic
     plus strand, track 4 should dominate at plus-strand arms and track 5 at
     minus-strand arms. The fraction of arms where this holds is reported; it
     should be well above 50% if the tracks are what we think they are.

Grouping follows iSTART (intronic = klasa 'intron genu'; active host =
host_signal >= 5). Intergenic pairs are the built-in negative control: their
"sense" arm is the plus-strand arm by convention, so their predicted asymmetry
should centre on zero.

Usage:
  python3 07_puffin_pairs.py --arms results/arms_K562.qc.tsv --line K562 \
      --genome /path/hg38.fa --puffin-dir /path/tools/puffin \
      --out-pred results/puffin_pairs_K562.tsv --out-summary results/puffin_summary_K562.tsv \
      [--signal-window 50] [--limit 300] [--no-cuda] [--n-boot 10000] [--seed 1]
"""
import argparse
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
WINDOW = 100000          # Puffin-D input length
TRACK_PRO_CAP_PLUS = 4   # per the authors' track order
TRACK_PRO_CAP_MINUS = 5
CALIB_RADIUS = 500       # +/- nt searched when calibrating coordinate alignment


def boot_median_ci(x, n_boot, rng):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan, np.nan
    meds = np.median(x[rng.integers(0, len(x), size=(n_boot, len(x)))], axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", required=True, help="output of 04_arm_qc.py (or 03_arm_coordinates.py)")
    ap.add_argument("--line", required=True)
    ap.add_argument("--genome", required=True, help="hg38 fasta (UCSC naming, soft-masked is fine)")
    ap.add_argument("--puffin-dir", required=True, help="clone of jzhoulab/puffin holding resources/puffin_D.pth")
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--signal-window", type=int, default=50,
                    help="half-width in nt for reading predicted signal at an arm; 50 matches iSTART")
    ap.add_argument("--limit", type=int, default=0, help="score only the first N pairs (0 = all)")
    ap.add_argument("--no-cuda", action="store_true")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    arms_path = os.path.abspath(a.arms)
    genome_path = os.path.abspath(a.genome)
    out_pred = os.path.abspath(a.out_pred)
    out_summary = os.path.abspath(a.out_summary)
    puffin_dir = os.path.abspath(a.puffin_dir)

    d = pd.read_csv(arms_path, sep="\t")
    for col in ("chrom", "mid", "pos_sense", "pos_antisense", "strand_sense", "klasa", "status"):
        if col not in d.columns:
            sys.exit(f"{arms_path}: missing column {col}")
    d = d[d["status"] == "ok"].reset_index(drop=True)
    if a.limit:
        d = d.head(a.limit).copy()
    print(f"[{a.line}] {len(d)} pairs to score from {arms_path}")

    from pyfaidx import Fasta
    fa = Fasta(genome_path, as_raw=True, sequence_always_upper=False)

    # the model resolves resources/ relative to the working directory
    os.chdir(puffin_dir)
    sys.path.insert(0, puffin_dir)
    from puffin_D import Puffin_D
    t0 = time.time()
    model = Puffin_D(use_cuda=not a.no_cuda)
    print(f"[{a.line}] model loaded in {time.time() - t0:.1f} s (cuda={not a.no_cuda})")

    hw = a.signal_window
    rows = []
    t_start = time.time()
    for i, r in enumerate(d.itertuples(index=False)):
        chrom = r.chrom
        if chrom not in fa:
            rows.append({"pred_status": "no_such_contig"})
            continue
        clen = len(fa[chrom])
        start = int(r.mid) - WINDOW // 2
        start = max(0, min(start, clen - WINDOW))
        if clen < WINDOW:
            rows.append({"pred_status": "contig_shorter_than_window"})
            continue
        seq = str(fa[chrom][start:start + WINDOW])
        # hg38 from UCSC is soft-masked: repeats are lowercase and would otherwise
        # be encoded as unknown bases. Introns are repeat-rich, so this matters.
        seq = seq.upper()
        n_frac = seq.count("N") / len(seq)

        i_sense = int(r.pos_sense) - start
        i_anti = int(r.pos_antisense) - start
        if not (hw <= i_sense < WINDOW - hw and hw <= i_anti < WINDOW - hw):
            rows.append({"pred_status": "arm_outside_window"})
            continue

        arr = np.asarray(model.predict(seq), dtype=np.float32)
        if arr.ndim == 3:
            arr = arr[0]
        plus, minus = arr[TRACK_PRO_CAP_PLUS], arr[TRACK_PRO_CAP_MINUS]

        t_sense = plus if r.strand_sense == "+" else minus
        t_anti = minus if r.strand_sense == "+" else plus
        w_s = t_sense[i_sense - hw:i_sense + hw + 1]
        w_a = t_anti[i_anti - hw:i_anti + hw + 1]
        pred_sense, pred_anti = float(w_s.max()), float(w_a.max())
        # three further readouts, so the conclusion does not rest on one choice of metric
        sum_sense, sum_anti = float(w_s.sum()), float(w_a.sum())
        pt_sense, pt_anti = float(t_sense[i_sense]), float(t_anti[i_anti])
        WW = 250
        wide_sense = float(t_sense[max(0, i_sense - WW):i_sense + WW + 1].max())
        wide_anti = float(t_anti[max(0, i_anti - WW):i_anti + WW + 1].max())

        # calibration: where is the predicted maximum relative to the arm, and does
        # the expected track dominate at each arm?
        lo_s, lo_a = max(0, i_sense - CALIB_RADIUS), max(0, i_anti - CALIB_RADIUS)
        off_sense = int(np.argmax(t_sense[lo_s:i_sense + CALIB_RADIUS + 1])) + lo_s - i_sense
        off_anti = int(np.argmax(t_anti[lo_a:i_anti + CALIB_RADIUS + 1])) + lo_a - i_anti
        i_plus = int(r.pos_plus) - start if not pd.isna(r.pos_plus) else None
        i_minus = int(r.pos_minus) - start if not pd.isna(r.pos_minus) else None
        plus_wins = minus_wins = np.nan
        if i_plus is not None and hw <= i_plus < WINDOW - hw:
            plus_wins = float(plus[i_plus - hw:i_plus + hw + 1].max()
                              > minus[i_plus - hw:i_plus + hw + 1].max())
        if i_minus is not None and hw <= i_minus < WINDOW - hw:
            minus_wins = float(minus[i_minus - hw:i_minus + hw + 1].max()
                               > plus[i_minus - hw:i_minus + hw + 1].max())

        local = seq[max(0, i_sense - 500):i_sense + 500]
        gc = (local.count("G") + local.count("C")) / max(1, len(local))

        rows.append({
            "pred_status": "ok", "window_start": start, "n_frac": n_frac, "gc_local": gc,
            "pred_sense": pred_sense, "pred_antisense": pred_anti,
            "pred_diff": pred_sense - pred_anti,
            "pred_DI": (pred_sense - pred_anti) / math.log(2),
            "pred_sense_sum": sum_sense, "pred_antisense_sum": sum_anti,
            "pred_DI_sum": math.log2((sum_sense + 1) / (sum_anti + 1)),
            "pred_sense_point": pt_sense, "pred_antisense_point": pt_anti,
            "pred_DI_point": (pt_sense - pt_anti) / math.log(2),
            "pred_sense_wide": wide_sense, "pred_antisense_wide": wide_anti,
            "pred_DI_wide": (wide_sense - wide_anti) / math.log(2),
            "pred_window_sum_plus": float(plus.sum()), "pred_window_sum_minus": float(minus.sum()),
            "calib_offset_sense": off_sense, "calib_offset_antisense": off_anti,
            "track_plus_wins_at_plus_arm": plus_wins, "track_minus_wins_at_minus_arm": minus_wins,
        })
        if (i + 1) % 500 == 0:
            el = time.time() - t_start
            print(f"  {i + 1}/{len(d)} pairs, {el:.0f} s elapsed, {el / (i + 1) * 1000:.0f} ms per pair", flush=True)

    pred = pd.concat([d.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    pred.to_csv(out_pred, sep="\t", index=False)
    ok = pred[pred["pred_status"] == "ok"].copy()
    print(f"\n[{a.line}] scored {len(ok)} / {len(pred)} pairs in {time.time() - t_start:.0f} s")
    if pred["pred_status"].nunique() > 1:
        print(pred["pred_status"].value_counts().to_string())

    print(f"\n== check 1: coordinate alignment (median offset of the predicted maximum) ==")
    for c in ("calib_offset_sense", "calib_offset_antisense"):
        v = ok[c]
        print(f"  {c}: median {v.median():.0f} nt, IQR {v.quantile(.25):.0f} to {v.quantile(.75):.0f}, "
              f"|offset| < 100 nt for {100 * (v.abs() < 100).mean():.1f}% of arms")
    print(f"\n== check 2: track identity (does the expected strand dominate at each arm?) ==")
    print(f"  track 4 (plus) stronger at plus-strand arms:  {100 * ok['track_plus_wins_at_plus_arm'].mean():.1f}%")
    print(f"  track 5 (minus) stronger at minus-strand arms: {100 * ok['track_minus_wins_at_minus_arm'].mean():.1f}%")
    print(f"\n== check 3: do predictions track the measured signal at all? ==")
    for obs, prd, lab in (("sens", "pred_sense", "sense arm"), ("anty", "pred_antisense", "antisense arm")):
        if obs in ok.columns:
            rho, p = spearmanr(ok[obs], ok[prd])
            print(f"  Spearman(observed {lab}, predicted {lab}): rho = {rho:.3f}, p = {p:.1e}")
    if "DI" in ok.columns:
        rho, p = spearmanr(ok["DI"], ok["pred_DI"])
        print(f"  Spearman(observed DI, predicted DI): rho = {rho:.3f}, p = {p:.1e}")

    groups = {
        "intergenic": ok[ok["klasa"] != INTRON_CLASS],
        "intron_all": ok[ok["klasa"] == INTRON_CLASS],
        "intron_silent": ok[(ok["klasa"] == INTRON_CLASS) & (ok["host_signal"] < ACTIVE_MIN_HOST_SIGNAL)],
        "intron_active": ok[(ok["klasa"] == INTRON_CLASS) & (ok["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL)],
    }
    srows = []
    ref = groups["intergenic"]["pred_DI"]
    for name, g in groups.items():
        if not len(g):
            continue
        lo, hi = boot_median_ci(g["pred_DI"], a.n_boot, rng)
        med = float(g["pred_DI"].median())
        p = (mannwhitneyu(g["pred_DI"], ref, alternative="two-sided").pvalue
             if name != "intergenic" and len(g) >= 5 and len(ref) >= 5 else np.nan)
        srows.append({
            "line": a.line, "group": name, "n": len(g),
            "median_pred_DI": med, "pred_DI_lo95": lo, "pred_DI_hi95": hi,
            "pred_sense_weaker_pct": 100 * (1 - 2 ** med),
            "pred_sense_weaker_pct_lo95": 100 * (1 - 2 ** hi),
            "pred_sense_weaker_pct_hi95": 100 * (1 - 2 ** lo),
            "pct_pred_DI_negative": 100 * float((g["pred_DI"] < 0).mean()),
            "median_pred_DI_sum": float(g["pred_DI_sum"].median()),
            "median_pred_DI_point": float(g["pred_DI_point"].median()),
            "median_pred_DI_wide": float(g["pred_DI_wide"].median()),
            "median_obs_DI": float(g["DI"].median()) if "DI" in g.columns else np.nan,
            "mannwhitney_p_vs_intergenic": p,
            "median_gc_local": float(g["gc_local"].median()),
            "test": "two-sided Mann-Whitney U on pred_DI vs intergenic; CI percentile bootstrap",
        })
    summary = pd.DataFrame(srows)
    summary.to_csv(out_summary, sep="\t", index=False)
    print("\n== predicted strand asymmetry by group ==")
    print(summary.drop(columns=["test"]).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nCIs: percentile bootstrap, {a.n_boot} resamples, seed {a.seed}.")
    print(f"written: {out_pred}; {out_summary}")


if __name__ == "__main__":
    main()
