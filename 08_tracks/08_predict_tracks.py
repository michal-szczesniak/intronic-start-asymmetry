#!/usr/bin/env python3
"""Step 08: strand-symmetric Puffin-D PRO-cap predictions along whole chromosomes.

First step of the genome-wide H1 test. Pairs from iSTART were too few on the
chromosomes Puffin was not trained on (chr8, chr9) to tell a sequence effect from
memorisation of the training signal. Here the model scores every position of the
chosen chromosomes, so the next step can compare predicted initiation on the
sense and antisense strand along whole introns, with millions of positions.

Tiling. Windows of 100 kb are placed every 50 kb; from each window only the
central 50 kb is kept (the first window also keeps its left flank and the last
its right flank), so every position, except within 25 kb of a chromosome end, is
predicted with at least 25 kb of sequence context on both sides, and exactly once.

Strand symmetry. As in step 07c, each window is scored as given and as its reverse
complement; the second output is mapped back (reversed, strands swapped) and the
two are averaged.

Windows with more than --max-n of N bases (gaps, centromeres) are not scored and
are stored as NaN.

Output, per chromosome, in --outdir (large; belongs on scratch):
  {chrom}.procap_plus.f16.npy, {chrom}.procap_minus.f16.npy
      float16, one value per genomic position (0-based), NaN where not scored
  {chrom}.windows.tsv
      per window: start, kept range, N fraction, and the Pearson correlation
      between forward and reverse-complement predictions (model symmetry QC)

Optional cross-check (--check-pairs, the per-pair table of step 07c): the signal
at each iSTART arm on these chromosomes is read from the new tracks exactly as
in 07c and compared with the 07c values. Windows differ (tiled here, centred on
the pair in 07c), so values agree closely but not exactly.

After writing, every array is re-opened and re-read, so a storage fault shows up
here and not in the next step.

Usage:
  python3 08_predict_tracks.py --genome /path/hg38.fa --puffin-dir /path/tools/puffin \
      --chroms chr8,chr9 --outdir /scratch/.../08_tracks \
      [--check-pairs ../07_puffin_pairs/results/puffin_rc_pairs_K562.tsv] [--no-cuda]
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

WINDOW = 100000
STRIDE = 50000
FLANK = (WINDOW - STRIDE) // 2
TRACK_PRO_CAP_PLUS = 4
TRACK_PRO_CAP_MINUS = 5
HW = 50  # half-width for reading signal at an arm, as in steps 07 and 07c
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def tracks(arr):
    arr = np.asarray(arr, dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[0]
    return arr[TRACK_PRO_CAP_PLUS], arr[TRACK_PRO_CAP_MINUS]


def window_plan(clen):
    """List of (start, keep_lo, keep_hi); kept ranges tile [0, clen) exactly once."""
    if clen < WINDOW:
        return []
    starts = list(range(0, clen - WINDOW + 1, STRIDE))
    if starts[-1] != clen - WINDOW:
        starts.append(clen - WINDOW)
    plan, prev_hi = [], 0
    for k, s in enumerate(starts):
        hi = clen if k == len(starts) - 1 else s + FLANK + STRIDE
        plan.append((s, prev_hi, hi))
        prev_hi = hi
    return plan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--puffin-dir", required=True)
    ap.add_argument("--chroms", required=True, help="comma-separated, e.g. chr8,chr9")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-n", type=float, default=0.5, help="skip windows with a larger fraction of N")
    ap.add_argument("--check-pairs", default="", help="per-pair table of step 07c, for a cross-check")
    ap.add_argument("--limit-windows", type=int, default=0, help="score only the first N windows (testing)")
    ap.add_argument("--no-cuda", action="store_true")
    a = ap.parse_args()

    genome_path = os.path.abspath(a.genome)
    puffin_dir = os.path.abspath(a.puffin_dir)
    outdir = os.path.abspath(a.outdir)
    check_path = os.path.abspath(a.check_pairs) if a.check_pairs else ""
    os.makedirs(outdir, exist_ok=True)
    chroms = a.chroms.split(",")

    from pyfaidx import Fasta
    fa = Fasta(genome_path, as_raw=True, sequence_always_upper=False)
    for c in chroms:
        if c not in fa:
            sys.exit(f"{c} not in {genome_path}")
    os.chdir(puffin_dir)  # the model resolves resources/ relative to the working directory
    sys.path.insert(0, puffin_dir)
    from puffin_D import Puffin_D
    model = Puffin_D(use_cuda=not a.no_cuda)

    for chrom in chroms:
        clen = len(fa[chrom])
        plan = window_plan(clen)
        if a.limit_windows:
            plan = plan[:a.limit_windows]
        f_plus = os.path.join(outdir, f"{chrom}.procap_plus.f16.npy")
        f_minus = os.path.join(outdir, f"{chrom}.procap_minus.f16.npy")
        m_plus = np.lib.format.open_memmap(f_plus, mode="w+", dtype=np.float16, shape=(clen,))
        m_minus = np.lib.format.open_memmap(f_minus, mode="w+", dtype=np.float16, shape=(clen,))
        m_plus[:] = np.nan
        m_minus[:] = np.nan
        print(f"[{chrom}] length {clen:,} nt, {len(plan)} windows", flush=True)
        wrows, t0 = [], time.time()
        for k, (s, lo, hi) in enumerate(plan):
            seq = str(fa[chrom][s:s + WINDOW]).upper()  # hg38 is soft-masked
            n_frac = seq.count("N") / WINDOW
            if n_frac > a.max_n:
                wrows.append({"start": s, "keep_lo": lo, "keep_hi": hi, "n_frac": n_frac,
                              "scored": 0, "r_fwd_rc": np.nan})
                continue
            pf, mf = tracks(model.predict(seq))
            rp, rm = tracks(model.predict(revcomp(seq)))
            pr, mr = rm[::-1], rp[::-1]          # RC track 5 -> genomic plus, RC track 4 -> genomic minus
            plus, minus = (pf + pr) / 2, (mf + mr) / 2
            m_plus[lo:hi] = plus[lo - s:hi - s].astype(np.float16)
            m_minus[lo:hi] = minus[lo - s:hi - s].astype(np.float16)
            x = np.concatenate([pf, mf])
            y = np.concatenate([pr, mr])
            r = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan
            wrows.append({"start": s, "keep_lo": lo, "keep_hi": hi, "n_frac": n_frac,
                          "scored": 1, "r_fwd_rc": r})
            if (k + 1) % 500 == 0:
                el = time.time() - t0
                print(f"  {k + 1}/{len(plan)} windows, {el:.0f} s", flush=True)
        m_plus.flush()
        m_minus.flush()
        del m_plus, m_minus
        w = pd.DataFrame(wrows)
        w.to_csv(os.path.join(outdir, f"{chrom}.windows.tsv"), sep="\t", index=False)

        # read back from storage before reporting anything
        p = np.load(f_plus, mmap_mode="r")
        m = np.load(f_minus, mmap_mode="r")
        covered = int(np.isfinite(p).sum())
        kept = int(sum(hi - lo for _, lo, hi in plan))
        sc = w[w["scored"] == 1]
        print(f"[{chrom}] {len(sc)}/{len(w)} windows scored in {time.time() - t0:.0f} s; "
              f"positions with a value: {covered:,} of {clen:,} ({100 * covered / clen:.1f}%); "
              f"kept ranges cover {kept:,} nt")
        if len(sc):
            print(f"  forward vs reverse-complement agreement per window: median r = "
                  f"{sc['r_fwd_rc'].median():.3f} (IQR {sc['r_fwd_rc'].quantile(.25):.3f}-"
                  f"{sc['r_fwd_rc'].quantile(.75):.3f})")
            fin = np.isfinite(p)
            pp = np.asarray(p[fin], dtype=np.float64)
            mm = np.asarray(m[np.isfinite(m)], dtype=np.float64)
            print(f"  plus strand:  mean {pp.mean():.4f}, 99.9th pct {np.percentile(pp, 99.9):.3f}, max {pp.max():.3f}")
            print(f"  minus strand: mean {mm.mean():.4f}, 99.9th pct {np.percentile(mm, 99.9):.3f}, max {mm.max():.3f}")

        if check_path:
            d = pd.read_csv(check_path, sep="\t")
            d = d[(d["chrom"] == chrom) & (d["pred_status"] == "ok")]
            got_s, got_a, ref_s, ref_a = [], [], [], []
            for r in d.itertuples(index=False):
                ts = p if r.strand_sense == "+" else m
                ta = m if r.strand_sense == "+" else p
                i_s, i_a = int(r.pos_sense), int(r.pos_antisense)
                ws = np.asarray(ts[i_s - HW:i_s + HW + 1], dtype=np.float64)
                wa = np.asarray(ta[i_a - HW:i_a + HW + 1], dtype=np.float64)
                if not (np.isfinite(ws).all() and np.isfinite(wa).all()):
                    continue
                got_s.append(ws.max()); got_a.append(wa.max())
                ref_s.append(r.pred_sense); ref_a.append(r.pred_antisense)
            if got_s:
                gs, ga, rs, ra = map(np.asarray, (got_s, got_a, ref_s, ref_a))
                di_new, di_ref = gs - ga, rs - ra
                print(f"  cross-check vs step 07c on {len(gs)} pairs: "
                      f"Spearman sense arm {spearmanr(gs, rs)[0]:.3f}, antisense arm {spearmanr(ga, ra)[0]:.3f}, "
                      f"DI {spearmanr(di_new, di_ref)[0]:.3f}; "
                      f"median |difference| at arms {np.median(np.abs(np.r_[gs - rs, ga - ra])):.4f}")
            else:
                print("  cross-check: no pairs on this chromosome")
        del p, m
    print(f"\nwritten to {outdir}")


if __name__ == "__main__":
    main()
