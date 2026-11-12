#!/usr/bin/env python3
"""Step 16c: genome-wide ProCapNet tracks for one chromosome, real and shuffled, from the model
that never saw that chromosome.

Steps 13c/13d (Puffin-D) found that weak-site results differ between the chromosomes the
model was trained on and the held-out ones (chr8, 9, 10 are the top three of 22). ProCapNet
was trained 7 times with chromosome-split test sets, so every chromosome is a test
chromosome of exactly one model. This step redoes the step-10 input for one chromosome
with that model:

  1. eligible positions as in step 10 (its own functions: protein-coding / lncRNA gene
     bodies minus exons +/-100 nt, TSSs +/-1 kb, other genes; intergenic >10 kb from genes),
     with "usable" = not N (Puffin-D's track coverage no longer applies);
  2. the same composition-preserving shuffle (step 10's shuffle_chrom, random stream
     default_rng([seed, 0]) as in 13b), with the base and dinucleotide counts checked;
  3. predictions for the real and the shuffled chromosome: windows of 2114 nt whose
     central 1000 nt tile the chromosome without overlap, each scored as given and as its
     reverse complement and averaged (16b's merged_counts). Real: every window that is not
     all N (the site thresholds are quantiles of all real predictions on chr8+chr9, as for
     Puffin-D); shuffled: only windows with eligible positions (the rest is unchanged);
  4. tracks of predicted counts per base, plus and minus strand, saved as float16:
     <outdir>/<chrom>.pc_{real,shuf1}_{plus,minus}.f16.npy, and <chrom>.pc_done.tsv last.
Site counting (thresholds from chr8+chr9, as Puffin-D) is step 16d, on the CPU.

Usage (GPU):
  python3 16c_procapnet_tracks.py --chrom chr8 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --models-dir .../ENCFF976FHE --cv-dir .../ENCFF033LJA --outdir $SCRATCH/16c_tracks \
      --step10 ../10_shuffle/10_shuffle_test.py --step16b 16b_score_pairs_procapnet.py [--batch 256]
"""
import argparse
import hashlib
import importlib.util
import os
import sys
import time

import numpy as np


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--cv-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--step10", required=True, help="path to 10_shuffle_test.py (masks, shuffle)")
    ap.add_argument("--step16b", required=True, help="path to 16b_score_pairs_procapnet.py (model, RC averaging)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--chunk", type=int, default=4096, help="windows prepared at a time")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--fake-model", default="", help=argparse.SUPPRESS)
    a = ap.parse_args()
    c = a.chrom
    os.makedirs(a.outdir, exist_ok=True)
    s10 = load_module("s10", os.path.abspath(a.step10))
    s16 = load_module("s16b", os.path.abspath(a.step16b))
    IN_W, OUT_W, TRIM, LUT = s16.IN_W, s16.OUT_W, s16.TRIM, s16.LUT
    t0 = time.time()

    folds = s16.load_folds(a.models_dir, a.cv_dir)
    f = [x for x in folds if c in x["test"]]
    if len(f) != 1:
        sys.exit(f"STOP: {c} is in the test set of {len(f)} folds")
    f = f[0]
    if a.fake_model:
        predict = s16.fake_predictor(f["path"], motif=a.fake_model)
        dev = "cpu (fake model)"
    else:
        import torch
        dev = ("cuda" if torch.cuda.is_available() else "cpu") if a.device == "auto" else a.device
        predict = s16.torch_predictor(f["path"], dev, 0)
        if dev == "cuda":
            dev += f" ({torch.cuda.get_device_name(0)})"
    print(f"[{c}] held-out model: fold {f['k']} (test {','.join(sorted(f['test']))}), {os.path.basename(f['path'])}; "
          f"device {dev}", flush=True)

    from pyfaidx import Fasta
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
    is_rep = raw >= 97
    up = np.where(is_rep, raw - 32, raw).astype(np.uint8)
    del raw
    L = len(up)
    ok = up != ord("N")
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), {c})
    glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
    elig = ig.copy()
    for _, s, e, sel in glist:
        elig[s:e] |= sel
    print(f"[{c}] {L:,} nt; {len(glist):,} genes with eligible positions; eligible {int(elig.sum()):,} nt "
          f"(intergenic {int(ig.sum()):,})  [{time.time() - t0:.0f} s]", flush=True)

    rng = np.random.default_rng([a.seed, 0])
    shuf, shuffled_mask = s10.shuffle_chrom(up, is_rep, glist, ig, rng)
    m0, d0 = s10.composition(up)
    m1, d1 = s10.composition(shuf)
    same = bool(np.array_equal(m0, m1) and np.array_equal(d0, d1))
    print(f"[{c}] shuffled {int(shuffled_mask.sum()):,} nt; base and dinucleotide counts identical: {same}; "
          f"positions changed {100 * float((shuf != up)[shuffled_mask].mean()):.1f}%; "
          f"md5 {hashlib.md5(shuf.tobytes()).hexdigest()}  [{time.time() - t0:.0f} s]", flush=True)
    if not same:
        sys.exit("STOP: composition check failed")

    n_win = (L + OUT_W - 1) // OUT_W
    starts = np.arange(0, L, OUT_W)
    win_sets = {"real": np.flatnonzero(np.add.reduceat(ok.astype(np.int64), starts) > 0),
                "shuf1": np.flatnonzero(np.add.reduceat(elig.astype(np.int64), starts) > 0)}
    print(f"[{c}] {n_win:,} windows; scored: real {len(win_sets['real']):,} (not all N), "
          f"shuffled {len(win_sets['shuf1']):,} (with eligible positions)", flush=True)
    pad_l, pad_r = TRIM, TRIM + OUT_W
    for tag, seq in (("real", up), ("shuf1", shuf)):
        wins = win_sets[tag]
        t1 = time.time()
        codes = np.concatenate([np.full(pad_l, -1, np.int8), LUT[seq], np.full(pad_r, -1, np.int8)])
        P = np.full(L, np.nan, dtype=np.float32)
        M = np.full(L, np.nan, dtype=np.float32)
        ar = np.arange(IN_W)
        for i in range(0, len(wins), a.chunk):
            w = wins[i:i + a.chunk]
            cw = codes[(w * OUT_W)[:, None] + ar[None, :]]           # input start = w*OUT_W - TRIM (+ pad_l)
            X = (cw[:, None, :] == np.arange(4, dtype=np.int8)[None, :, None]).astype(np.float32)
            out = s16.merged_counts(predict, X, a.batch)["avg"]
            for j, wj in enumerate(w):
                lo, hi = wj * OUT_W, min(L, wj * OUT_W + OUT_W)
                P[lo:hi] = out[j, 0, :hi - lo]
                M[lo:hi] = out[j, 1, :hi - lo]
            if (i // a.chunk) % 10 == 0:
                print(f"  [{c} {tag}] {min(i + a.chunk, len(wins)):,}/{len(wins):,} windows  [{time.time() - t1:.0f} s]",
                      flush=True)
        for strand, T in (("plus", P), ("minus", M)):
            np.save(os.path.join(a.outdir, f"{c}.pc_{tag}_{strand}.f16.npy"), T.astype(np.float16))
        fin = np.isfinite(P)
        print(f"[{c} {tag}] predicted {int(fin.sum()):,} nt in {time.time() - t1:.0f} s; plus-strand counts "
              f"quantiles 50/99/99.9%: " + ", ".join(f"{np.quantile(P[fin], q):.3g}" for q in (0.5, 0.99, 0.999)), flush=True)
        del codes, P, M
    with open(os.path.join(a.outdir, f"{c}.pc_done.tsv"), "w") as fh:
        fh.write("chrom\tlength\tfold\twindows_scored\teligible_nt\tseconds\n"
                 f"{c}\t{L}\t{f['k']}\t{len(win_sets['real'])}\t{int(elig.sum())}\t{time.time() - t0:.0f}\n")
    print(f"[{c}] done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
