#!/usr/bin/env python3
"""Step 20a: where in the sequence of an iSTART pair does ProCapNet find the weaker sense start?

ProCapNet, reading only the 2114 nt around a pair and never trained on its chromosome, predicts
about half of the measured weakness of sense starts in intronic pairs (16b + 07b). This step
asks whether that prediction comes from base composition or from the arrangement of the
sequence, and whether it sits in the element itself or in its surroundings. Each pair's input
window is re-scored after a shuffle that keeps the count of every base and of every
dinucleotide (Altschul-Erickson, step 10's dinuc_shuffle), in one of three regions:
  window   the whole 2114-nt input (composition kept, arrangement lost everywhere);
  element  from 100 nt (--element-pad) before the first arm to 100 nt after the second, i.e. both
           core promoters and the space between them; the flanks stay as they are;
  flanks   everything except the element (left and right flank shuffled separately).
Runs of N are left in place (only stretches of A/C/G/T are shuffled). Each region is shuffled
with --seeds different seeds; the predicted counts are averaged over the seeds, so the readout
is the expected prediction for a sequence of that composition. "none" = the original window,
scored in the same run; it must reproduce 16b's held-out prediction (check printed).

Every pair is scored by the ProCapNet model that never saw its chromosome (as 16b), as given and
as its reverse complement (16b's merged_counts), and read out with 16b's readouts (pred_DI etc.).
Output: <outdir>/procapnet_shuf_<region>_pairs_<line>.tsv for region in none, window, element,
flanks: the rows and columns of 16b's held-out file with the predictions replaced, so
07b_calibrate.py runs on them unchanged.

Usage (GPU):
  python3 20a_pairs_shuffle_procapnet.py --pairs-dir ../16_procapnet/results --lines K562,HCT116 \\
      --genome hg38.fa --models-dir .../ENCFF976FHE --cv-dir .../ENCFF033LJA \\
      --step16b ../16_procapnet/16b_score_pairs_procapnet.py --step10 ../10_shuffle/10_shuffle_test.py \\
      --outdir results [--seeds 5] [--element-pad 100] [--batch 256]
"""
import argparse
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

REGIONS = ("none", "window", "element", "flanks")
DROP = ("pred_DI_fwd", "pred_DI_rc", "pc_offset_sense", "pc_offset_antisense", "pc_plus_wins_at_plus_arm",
        "pc_minus_wins_at_minus_arm")
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def shuffle_runs(seq, lo, hi, rng, s10):
    """Dinucleotide shuffle of seq[lo:hi] (uint8 ACGTN), run by run of A/C/G/T, in place."""
    seg = seq[lo:hi]
    isb = np.isin(seg, BASES)
    edges = np.flatnonzero(np.diff(np.r_[0, isb.astype(np.int8), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        if b - a >= 3:
            seg[a:b] = s10.dinuc_shuffle(seg[a:b], rng)


def dinuc_counts(x):
    m = np.isin(x[:-1], BASES) & np.isin(x[1:], BASES)
    return np.bincount((np.searchsorted(BASES, x[:-1][m]) * 4 + np.searchsorted(BASES, x[1:][m])), minlength=16)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-dir", required=True, help="directory with 16b's procapnet_heldout_pairs_<line>.tsv")
    ap.add_argument("--lines", default="K562,HCT116")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--cv-dir", required=True)
    ap.add_argument("--step16b", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--element-pad", type=int, default=100)
    ap.add_argument("--signal-window", type=int, default=50)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--fake-model", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args()
    s16 = load_module("s16b", os.path.abspath(a.step16b))
    s10 = load_module("s10", os.path.abspath(a.step10))
    IN_W, OUT_W, TRIM = s16.IN_W, s16.OUT_W, s16.TRIM
    hw, E = a.signal_window, a.element_pad
    os.makedirs(a.outdir, exist_ok=True)
    folds = s16.load_folds(a.models_dir, a.cv_dir)
    print("folds (test chromosomes): " + "; ".join(f"{f['k']}: {','.join(sorted(f['test']))}" for f in folds))
    if a.fake_model:
        get_pred = lambda f: s16.fake_predictor(f["path"])
        print("device: cpu (fake model)")
    else:
        import torch
        dev = ("cuda" if torch.cuda.is_available() else "cpu") if a.device == "auto" else a.device
        get_pred = lambda f: s16.torch_predictor(f["path"], dev, 0)
        print(f"torch {torch.__version__}, device {dev}" + (f" ({torch.cuda.get_device_name(0)})" if dev == "cuda" else ""))
    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=True)
    fold_of = {c: f for f in folds for c in f["test"]}

    for line in [x for x in a.lines.split(",") if x]:
        t0 = time.time()
        f_in = os.path.join(a.pairs_dir, f"procapnet_heldout_pairs_{line}.tsv")
        if not os.path.isfile(f_in):
            sys.exit(f"STOP: missing {f_in}")
        d = pd.read_csv(f_in, sep="\t", keep_default_na=False, na_values=["", "NA", "nan", "NaN"],
                        dtype={"host": str, "klasa": str, "strand_sense": str, "chrom": str})
        n = len(d)
        ok = np.flatnonzero((d["pred_status"] == "ok").to_numpy())
        ps, pa = d["pos_sense"].to_numpy(np.int64), d["pos_antisense"].to_numpy(np.int64)
        center = (ps + pa) // 2
        in_start = center - IN_W // 2
        i_s, i_a = ps - (in_start + TRIM), pa - (in_start + TRIM)
        chroms = d["chrom"].astype(str).to_numpy()
        ss = d["strand_sense"].astype(str).to_numpy()
        # windows (uint8 A/C/G/T/N) of the scored pairs
        W = np.full((n, IN_W), ord("N"), dtype=np.uint8)
        for c in np.unique(chroms[ok]):
            if c not in fold_of:
                sys.exit(f"STOP: {c} has status ok in 16b but is in no test fold")
            seqc = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
            for j in ok[chroms[ok] == c]:
                lo, hi = in_start[j], in_start[j] + IN_W
                x0, x1 = max(lo, 0), min(hi, len(seqc))
                W[j, x0 - lo:x1 - lo] = seqc[x0:x1]
        W[~np.isin(W, BASES)] = ord("N")
        el_lo = np.clip(np.minimum(i_s, i_a) + TRIM - E, 0, IN_W)
        el_hi = np.clip(np.maximum(i_s, i_a) + TRIM + E + 1, 0, IN_W)
        print(f"\n######## {line}: {n:,} pairs, {len(ok):,} scored by 16b; element (arms +/- {E} nt) median "
              f"{int(np.median((el_hi - el_lo)[ok]))} nt; {a.seeds} seeds per region", flush=True)

        C = {r: np.zeros((n, 2, OUT_W), np.float32) for r in REGIONS}
        changed = {r: [] for r in REGIONS[1:]}
        for f in folds:
            sel = ok[np.isin(chroms[ok], list(f["test"]))]
            if not len(sel):
                continue
            t1 = time.time()
            pred = get_pred(f)
            for r in REGIONS:
                for s in range(1 if r == "none" else a.seeds):
                    X = W[sel].copy()
                    if r != "none":
                        for i, j in enumerate(sel):
                            rng = np.random.default_rng([s, REGIONS.index(r), j])
                            if r == "window":
                                shuffle_runs(X[i], 0, IN_W, rng, s10)
                            elif r == "element":
                                shuffle_runs(X[i], el_lo[j], el_hi[j], rng, s10)
                            else:
                                shuffle_runs(X[i], 0, el_lo[j], rng, s10)
                                shuffle_runs(X[i], el_hi[j], IN_W, rng, s10)
                            if s == 0 and (not np.array_equal(dinuc_counts(X[i]), dinuc_counts(W[j])) or
                                           not np.array_equal(np.sort(X[i]), np.sort(W[j]))):
                                sys.exit(f"STOP: composition changed by the shuffle (pair {j}, region {r})")
                        changed[r].append((X != W[sel]).mean(axis=1))
                    oh = (X[:, None, :] == BASES[None, :, None]).astype(np.float32)
                    C[r][sel] += s16.merged_counts(pred, oh, a.batch)["avg"] / (1 if r == "none" else a.seeds)
            print(f"  fold {f['k']}: {len(sel):,} pairs  [{time.time() - t1:.0f} s]", flush=True)

        for r in REGIONS:
            ro = s16.readouts(C[r][ok], ss[ok], i_s[ok], i_a[ok], hw)
            D = d.drop(columns=[x for x in DROP if x in d.columns]).copy()
            for k, v in ro.items():
                col = np.full(n, np.nan)
                col[ok] = v
                D[k] = col
            D["pc_total_counts"] = np.where(D["pred_status"] == "ok", C[r].sum(axis=(1, 2)), np.nan)
            D["shuffle_region"] = r
            D["shuffle_seeds"] = 0 if r == "none" else a.seeds
            out = os.path.join(a.outdir, f"procapnet_shuf_{r}_pairs_{line}.tsv")
            D.to_csv(out, sep="\t", index=False)
            intr = (D["klasa"] == "intron genu").to_numpy()[ok]
            msg = f"  {r:8s} median pred_DI intronic {np.median(ro['pred_DI'][intr]):+.3f}, intergenic " \
                  f"{np.median(ro['pred_DI'][~intr]):+.3f}"
            if r != "none":
                msg += f"; positions changed {100 * np.mean(np.concatenate(changed[r])):.0f}% of the window"
            print(msg + f"  -> {out}")
            if r == "none":
                dev_ = np.abs(ro["pred_DI"] - d["pred_DI"].to_numpy(float)[ok])
                print(f"  check: 'none' vs 16b held-out pred_DI, max |difference| {np.nanmax(dev_):.2e} "
                      f"({'OK' if np.nanmax(dev_) < 1e-3 else 'FAILED'})")
        print(f"  [{line}] done in {(time.time() - t0) / 60:.1f} min")
    print("\nnext: 07b_calibrate.py on each procapnet_shuf_<region>_pairs_<line>.tsv")


if __name__ == "__main__":
    main()
