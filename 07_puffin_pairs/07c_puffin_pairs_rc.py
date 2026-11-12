#!/usr/bin/env python3
"""Step 07c: step 07 again, with strand-symmetric (reverse-complement averaged) predictions.

Step 07b found that Puffin-D reproduces part of the intronic strand asymmetry, but
only for pairs in minus-strand genes; for plus-strand genes the predicted effect
was about zero, while the measured effect is present on both strands. Genes have
no preferred orientation in the genome, so a real sequence effect should look the
same on both strands. The likely cause is the model itself: it reads the plus
strand of the window and emits both strands, and nothing forces it to treat a
sequence and its reverse complement identically.

The standard remedy is used here: each window is scored twice, as given and as its
reverse complement, the second output is mapped back to genomic coordinates
(reversed, strands swapped), and the two are averaged. The averaged prediction is
strand-symmetric by construction, so any remaining difference between plus- and
minus-strand genes cannot come from the model's orientation.

Mapping of the reverse-complement pass back to the genome, for window length L:
    genomic plus strand,  position i  <-  RC track 5 (reverse strand of input) at L-1-i
    genomic minus strand, position i  <-  RC track 4 (forward strand of input) at L-1-i

Checks reported before results:
  1. Alignment of the RC pass: the lag (in nt) that best aligns the forward and
     the mapped-back RC track around each sense arm. Median 0 = mapping is right.
  2. RC consistency of the model: Spearman between DI from the forward pass and DI
     from the RC pass. Near 1 = the model is already nearly symmetric.
  3. The diagnostic that motivated this step: median predicted DI by group and by
     sense-arm strand, for the forward pass, the RC pass and their average.
  4. As in step 07: offset of the predicted maximum and track identity, now on the
     averaged tracks.

The per-pair table uses the same column names as step 07 for the averaged
prediction (pred_DI, pred_DI_point, pred_DI_wide, pred_DI_sum), so 07b_calibrate.py
runs on it unchanged; pred_DI_fwd and pred_DI_rc are kept for diagnostics.

Usage:
  python3 07c_puffin_pairs_rc.py --arms ../02_istart_baseline/results/arms_K562.qc.tsv --line K562 \
      --genome /path/hg38.fa --puffin-dir /path/tools/puffin \
      --out-pred results/puffin_rc_pairs_K562.tsv [--signal-window 50] [--limit 300] [--no-cuda]
"""
import argparse
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

INTRON_CLASS = "intron genu"
ACTIVE_MIN_HOST_SIGNAL = 5
WINDOW = 100000
TRACK_PRO_CAP_PLUS = 4
TRACK_PRO_CAP_MINUS = 5
CALIB_RADIUS = 500
WIDE = 250
MAX_LAG = 5
LAG_HALF = 200
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def tracks(arr):
    arr = np.asarray(arr, dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[0]
    return arr[TRACK_PRO_CAP_PLUS], arr[TRACK_PRO_CAP_MINUS]


def readouts(plus, minus, strand_sense, i_s, i_a, hw):
    t_s = plus if strand_sense == "+" else minus
    t_a = minus if strand_sense == "+" else plus
    w_s, w_a = t_s[i_s - hw:i_s + hw + 1], t_a[i_a - hw:i_a + hw + 1]
    mx_s, mx_a = float(w_s.max()), float(w_a.max())
    su_s, su_a = float(w_s.sum()), float(w_a.sum())
    pt_s, pt_a = float(t_s[i_s]), float(t_a[i_a])
    wd_s = float(t_s[max(0, i_s - WIDE):i_s + WIDE + 1].max())
    wd_a = float(t_a[max(0, i_a - WIDE):i_a + WIDE + 1].max())
    return {
        "pred_sense": mx_s, "pred_antisense": mx_a,
        "pred_DI": (mx_s - mx_a) / math.log(2),
        "pred_DI_sum": math.log2((su_s + 1) / (su_a + 1)) if su_s > -1 and su_a > -1 else np.nan,
        "pred_DI_point": (pt_s - pt_a) / math.log(2),
        "pred_DI_wide": (wd_s - wd_a) / math.log(2),
    }, t_s, t_a


def best_lag(x, y, i):
    """Lag (nt) of y relative to x that maximises their correlation around position i."""
    lo, hi = i - LAG_HALF, i + LAG_HALF
    if lo - MAX_LAG < 0 or hi + MAX_LAG > len(x):
        return np.nan
    xs = x[lo:hi]
    if xs.std() == 0:
        return np.nan
    best, best_r = np.nan, -np.inf
    for lag in range(-MAX_LAG, MAX_LAG + 1):
        ys = y[lo + lag:hi + lag]
        if ys.std() == 0:
            continue
        r = float(np.corrcoef(xs, ys)[0, 1])
        if r > best_r:
            best, best_r = lag, r
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", required=True)
    ap.add_argument("--line", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--puffin-dir", required=True)
    ap.add_argument("--out-pred", required=True)
    ap.add_argument("--signal-window", type=int, default=50)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-cuda", action="store_true")
    a = ap.parse_args()

    arms_path = os.path.abspath(a.arms)
    genome_path = os.path.abspath(a.genome)
    out_pred = os.path.abspath(a.out_pred)
    puffin_dir = os.path.abspath(a.puffin_dir)

    d = pd.read_csv(arms_path, sep="\t")
    for col in ("chrom", "mid", "pos_sense", "pos_antisense", "strand_sense", "klasa", "status"):
        if col not in d.columns:
            sys.exit(f"{arms_path}: missing column {col}")
    d = d[d["status"] == "ok"].reset_index(drop=True)
    if a.limit:
        d = d.head(a.limit).copy()
    print(f"[{a.line}] {len(d)} pairs to score (forward + reverse complement) from {arms_path}")

    from pyfaidx import Fasta
    fa = Fasta(genome_path, as_raw=True, sequence_always_upper=False)
    os.chdir(puffin_dir)  # the model resolves resources/ relative to the working directory
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
        if clen < WINDOW:
            rows.append({"pred_status": "contig_shorter_than_window"})
            continue
        start = max(0, min(int(r.mid) - WINDOW // 2, clen - WINDOW))
        seq = str(fa[chrom][start:start + WINDOW]).upper()  # hg38 is soft-masked
        i_s, i_a = int(r.pos_sense) - start, int(r.pos_antisense) - start
        if not (hw <= i_s < WINDOW - hw and hw <= i_a < WINDOW - hw):
            rows.append({"pred_status": "arm_outside_window"})
            continue

        plus_f, minus_f = tracks(model.predict(seq))
        rc_plus_in, rc_minus_in = tracks(model.predict(revcomp(seq)))
        plus_r = rc_minus_in[::-1].copy()   # RC track 5 -> genomic plus
        minus_r = rc_plus_in[::-1].copy()   # RC track 4 -> genomic minus
        plus_a, minus_a = (plus_f + plus_r) / 2, (minus_f + minus_r) / 2

        ro_a, t_s, t_a = readouts(plus_a, minus_a, r.strand_sense, i_s, i_a, hw)
        ro_f, tf_s, _ = readouts(plus_f, minus_f, r.strand_sense, i_s, i_a, hw)
        ro_r, tr_s, _ = readouts(plus_r, minus_r, r.strand_sense, i_s, i_a, hw)

        lo_s, lo_a = max(0, i_s - CALIB_RADIUS), max(0, i_a - CALIB_RADIUS)
        off_s = int(np.argmax(t_s[lo_s:i_s + CALIB_RADIUS + 1])) + lo_s - i_s
        off_a = int(np.argmax(t_a[lo_a:i_a + CALIB_RADIUS + 1])) + lo_a - i_a
        i_plus = int(r.pos_plus) - start if not pd.isna(r.pos_plus) else None
        i_minus = int(r.pos_minus) - start if not pd.isna(r.pos_minus) else None
        plus_wins = minus_wins = np.nan
        if i_plus is not None and hw <= i_plus < WINDOW - hw:
            plus_wins = float(plus_a[i_plus - hw:i_plus + hw + 1].max()
                              > minus_a[i_plus - hw:i_plus + hw + 1].max())
        if i_minus is not None and hw <= i_minus < WINDOW - hw:
            minus_wins = float(minus_a[i_minus - hw:i_minus + hw + 1].max()
                               > plus_a[i_minus - hw:i_minus + hw + 1].max())

        local = seq[max(0, i_s - 500):i_s + 500]
        rows.append({
            "pred_status": "ok", "window_start": start, "n_frac": seq.count("N") / len(seq),
            "gc_local": (local.count("G") + local.count("C")) / max(1, len(local)),
            **ro_a,
            "pred_DI_fwd": ro_f["pred_DI"], "pred_DI_rc": ro_r["pred_DI"],
            "rc_lag_sense": best_lag(tf_s, tr_s, i_s),
            "calib_offset_sense": off_s, "calib_offset_antisense": off_a,
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

    print("\n== check 1: alignment of the reverse-complement pass (lag of RC vs forward track at sense arms) ==")
    lag = ok["rc_lag_sense"].dropna()
    print(f"  n = {len(lag)}, median lag {lag.median():.0f} nt, lag == 0 for {100 * (lag == 0).mean():.1f}%, "
          f"|lag| <= 1 for {100 * (lag.abs() <= 1).mean():.1f}%")
    print("\n== check 2: is the model already strand-symmetric? ==")
    rho = spearmanr(ok["pred_DI_fwd"], ok["pred_DI_rc"])[0]
    print(f"  Spearman(DI forward pass, DI reverse-complement pass) = {rho:.3f}")
    if "DI" in ok.columns:
        for c in ("pred_DI_fwd", "pred_DI_rc", "pred_DI"):
            print(f"  Spearman(observed DI, {c:11s}) = {spearmanr(ok['DI'], ok[c])[0]:.3f}")

    print("\n== check 3: median predicted DI by group and sense-arm strand (x 1000) ==")
    groups = {
        "intergenic": ok["klasa"] != INTRON_CLASS,
        "intron_active": (ok["klasa"] == INTRON_CLASS) & (ok["host_signal"] >= ACTIVE_MIN_HOST_SIGNAL),
        "intron_silent": (ok["klasa"] == INTRON_CLASS) & (ok["host_signal"] < ACTIVE_MIN_HOST_SIGNAL),
    }
    tab = []
    for gname, gm in groups.items():
        for st in ("+", "-"):
            g = ok[gm & (ok["strand_sense"] == st)]
            if not len(g):
                continue
            tab.append({"group": gname, "sense_strand": st, "n": len(g),
                        "fwd": 1000 * g["pred_DI_fwd"].median(), "rc": 1000 * g["pred_DI_rc"].median(),
                        "avg": 1000 * g["pred_DI"].median(),
                        "obs": g["DI"].median() if "DI" in g.columns else np.nan})
    print(pd.DataFrame(tab).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("  (fwd/rc/avg are predicted DI x 1000; obs is the measured DI, not scaled.)")

    print("\n== check 4: offset of the predicted maximum and track identity, averaged tracks ==")
    for c in ("calib_offset_sense", "calib_offset_antisense"):
        v = ok[c]
        print(f"  {c}: median {v.median():.0f} nt, IQR {v.quantile(.25):.0f} to {v.quantile(.75):.0f}")
    print(f"  plus track stronger at plus-strand arms:  {100 * ok['track_plus_wins_at_plus_arm'].mean():.1f}%")
    print(f"  minus track stronger at minus-strand arms: {100 * ok['track_minus_wins_at_minus_arm'].mean():.1f}%")
    print(f"\nwritten: {out_pred}")


if __name__ == "__main__":
    main()
