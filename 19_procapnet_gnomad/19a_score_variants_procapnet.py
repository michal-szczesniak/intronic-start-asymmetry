#!/usr/bin/env python3
"""Step 19a: score the gnomAD SNVs of step 11d with ProCapNet, using the model that never saw
the chromosome, so that the variant classes of the selection test can be defined by a model
that reads only the local sequence (Puffin-D's intronic predictions also carry the trace of the
host gene's transcription, step 17).

Variant set: exactly the rows of the step-11d score file of the chromosome (same eligible
positions, gene labels, strands and random lncRNA / intergenic subsets as the Puffin-D run);
only the predictions are replaced.

Prediction, as in step 16c: tiles whose central 1000 nt are predicted from 2114 nt of sequence
(N and the space beyond the chromosome ends = no base), each window as given and as its reverse
complement, averaged (16b's merged_counts), in float32. Readout, as in step 11d: the maximum of
the predicted counts within +/-50 nt of the variant, on the plus and on the minus strand.
Reference and alternative are always read from the same window:
  tile_kind 0: the 16c tile that contains the variant, when the readout fits in its central
               1000 nt (variant 50..949 nt into the tile); the reference then equals the 16c track;
  tile_kind 1: otherwise the tile shifted by 500 nt, in which the variant lies 450..549 nt in.
Every alternative allele is scored alone (one substitution per sequence); there is no screen.

Output <out>: gzipped table, one row per input row, in input order:
  chrom pos ref alt AC AN set gene_id strand ref_plus ref_minus alt_plus alt_minus fold tile_kind
(the column names of step 11d, so 12c / 19b read it unchanged). Restartable: each chunk of
tiles is written to <out>.parts/ and listed in <out>.done; the final file is assembled at the end.

Check (with --tracks-dir): for tile_kind 0 the reference readout must equal the readout of the
16c track (stored as float16) at the same position; quantiles of the relative difference printed.

Usage (GPU):
  python3 19a_score_variants_procapnet.py --chrom chr21 --puffin-scores $SCRATCH/11d_scores/chr21.scores.tsv.gz \\
      --genome hg38.fa --models-dir .../ENCFF976FHE --cv-dir .../ENCFF033LJA \\
      --step16b ../16_procapnet/16b_score_pairs_procapnet.py --out $SCRATCH/19a_scores/chr21.scores.tsv.gz \\
      [--tracks-dir $SCRATCH/16c_tracks] [--batch 256] [--chunk-tiles 2000] [--limit-chunks 0]
"""
import argparse
import gzip
import importlib.util
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

HW = 50
SHIFT = 500
OFFS = np.arange(-HW, HW + 1)
COLS_IN = ["chrom", "pos", "ref", "alt", "AC", "AN", "set", "gene_id", "strand"]
COLS_OUT = COLS_IN + ["ref_plus", "ref_minus", "alt_plus", "alt_minus", "fold", "tile_kind"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def readout(C, jj, rows=None):
    """Max of C[row, strand, jj-HW .. jj+HW] for plus and minus; rows maps each output to a row of C."""
    ix = jj[:, None] + OFFS[None, :]
    r = (np.arange(len(jj)) if rows is None else rows)[:, None]
    return C[r, 0, ix].max(axis=1), C[r, 1, ix].max(axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--puffin-scores", required=True, help="step-11d score file of this chromosome")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--cv-dir", required=True)
    ap.add_argument("--step16b", required=True, help="path to 16b_score_pairs_procapnet.py (model, RC averaging)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tracks-dir", default="", help="16c tracks, for the reference check (optional)")
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--chunk-tiles", type=int, default=2000, help="1000-nt tiles per restart unit")
    ap.add_argument("--block", type=int, default=2048, help="sequences prepared at a time")
    ap.add_argument("--limit-chunks", type=int, default=0, help="pilot: stop after this many chunks")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--seed", type=int, default=1, help="only for the sample of the track check")
    ap.add_argument("--fake-model", default="", help=argparse.SUPPRESS)
    a = ap.parse_args()
    c = a.chrom
    out = os.path.abspath(a.out)
    parts = out + ".parts"
    os.makedirs(parts, exist_ok=True)
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
    Lc = len(raw)
    PADL, PADR = TRIM + SHIFT, TRIM + OUT_W + SHIFT
    padded = np.concatenate([np.full(PADL, -1, np.int8), LUT[raw], np.full(PADR, -1, np.int8)])
    up = np.where(raw >= 97, raw - 32, raw).astype(np.uint8)
    del raw

    v = pd.read_csv(a.puffin_scores, sep="\t", usecols=COLS_IN, keep_default_na=False, na_values=[],
                    dtype={"chrom": "category", "pos": np.int64, "ref": "category", "alt": "category", "AC": np.int64,
                           "AN": np.int64, "set": "category", "gene_id": "category", "strand": "category"})
    v = v[COLS_IN]
    n = len(v)
    if n == 0:
        sys.exit(f"STOP: {a.puffin_scores} has no rows")
    if not (v["chrom"] == c).all():
        sys.exit(f"STOP: {a.puffin_scores} has rows of other chromosomes")
    pos = v["pos"].to_numpy()
    if pos.min() < 1 or pos.max() > Lc:
        sys.exit("STOP: positions outside the chromosome")
    if (np.diff(pos) < 0).any():
        sys.exit("STOP: input rows are not sorted by position (the output keeps the input order only then)")
    if any(len(x) != 1 for x in v["ref"].cat.categories) or any(len(x) != 1 for x in v["alt"].cat.categories):
        sys.exit("STOP: not all rows are single-nucleotide variants")
    p0 = pos - 1
    ref_b = np.frombuffer("".join(v["ref"].astype(str)).encode(), dtype=np.uint8)
    alt_b = np.frombuffer("".join(v["alt"].astype(str)).encode(), dtype=np.uint8)
    n_bad_ref = int((up[p0] != ref_b).sum())
    if n_bad_ref:
        sys.exit(f"STOP: reference allele differs from the genome at {n_bad_ref:,} rows")
    alt_code = LUT[alt_b]
    if (alt_code < 0).any() or (alt_b == ref_b).any():
        sys.exit("STOP: alternative allele not A/C/G/T or equal to the reference")

    k = p0 // OUT_W
    j = p0 - k * OUT_W
    kind0 = (j >= HW) & (j < OUT_W - HW)
    ostart = np.where(kind0, k * OUT_W, np.where(j < HW, k * OUT_W - SHIFT, k * OUT_W + SHIFT))
    jj = p0 - ostart
    if (jj < HW).any() or (jj >= OUT_W - HW).any():
        sys.exit("STOP: internal error, readout outside the tile")
    chunk = k // a.chunk_tiles
    chunk_ids = np.unique(chunk)
    print(f"[{c}] {Lc:,} nt; {n:,} variants (" + ", ".join(f"{s} {m:,}" for s, m in v["set"].value_counts().items()) +
          f"); tile_kind 1 (shifted tile): {100 * (~kind0).mean():.1f}%; {len(chunk_ids):,} chunks of "
          f"{a.chunk_tiles} tiles  [{time.time() - t0:.0f} s]", flush=True)

    tracks = None
    if a.tracks_dir:
        tp = [os.path.join(a.tracks_dir, f"{c}.pc_real_{s}.f16.npy") for s in ("plus", "minus")]
        if all(os.path.isfile(x) for x in tp):
            tracks = [np.load(x, mmap_mode="r") for x in tp]
            if len(tracks[0]) != Lc:
                sys.exit(f"STOP: 16c track length {len(tracks[0]):,} != chromosome length {Lc:,}")
        else:
            print(f"[{c}] no 16c tracks in {a.tracks_dir}; reference check skipped", flush=True)
    rng = np.random.default_rng(a.seed)
    chk = {"plus": [], "minus": []}

    def predict_codes(cw):
        X = (cw[:, None, :] == np.arange(4, dtype=np.int8)[None, :, None]).astype(np.float32)
        return s16.merged_counts(predict, X, a.batch)["avg"]

    done_path = out + ".done"
    done = set()
    if os.path.exists(done_path):
        done = {int(x) for x in open(done_path).read().split()}
    todo = [q for q in chunk_ids if q not in done]
    if a.limit_chunks:
        todo = todo[:a.limit_chunks]
    print(f"[{c}] {len(done):,} chunks already done; this run: {len(todo):,}", flush=True)
    order = np.argsort(chunk, kind="stable")
    bounds = np.searchsorted(chunk[order], chunk_ids)
    bounds = dict(zip(chunk_ids, zip(bounds, np.r_[bounds[1:], n])))
    ar = np.arange(IN_W)
    t1 = time.time()
    n_run = 0
    for ci, q in enumerate(todo):
        lo, hi = bounds[q]
        idx = np.sort(order[lo:hi])                          # input order within the chunk
        uo, inv = np.unique(ostart[idx], return_inverse=True)
        inv = inv.ravel()
        C_ref = np.concatenate([predict_codes(padded[(uo[i:i + a.block] - TRIM + PADL)[:, None] + ar[None, :]])
                                for i in range(0, len(uo), a.block)])
        rp, rm = readout(C_ref, jj[idx], rows=inv)
        ap_, am = np.empty(len(idx), np.float32), np.empty(len(idx), np.float32)
        for i in range(0, len(idx), a.block):
            sl = idx[i:i + a.block]
            cw = padded[(ostart[sl] - TRIM + PADL)[:, None] + ar[None, :]].copy()
            cw[np.arange(len(sl)), TRIM + jj[sl]] = alt_code[sl]
            ap_[i:i + len(sl)], am[i:i + len(sl)] = readout(predict_codes(cw), jj[sl])
        res = v.iloc[idx].copy()
        res["ref_plus"], res["ref_minus"], res["alt_plus"], res["alt_minus"] = rp, rm, ap_, am
        res["fold"] = f["k"]
        res["tile_kind"] = np.where(kind0[idx], 0, 1)
        vals = res[["ref_plus", "ref_minus", "alt_plus", "alt_minus"]].to_numpy()
        if not np.isfinite(vals).all():
            sys.exit(f"STOP: non-finite prediction in chunk {q}")
        tmp = os.path.join(parts, f"chunk_{q:06d}.tsv.gz.tmp")
        with gzip.open(tmp, "wt", compresslevel=3) as fh:
            res.to_csv(fh, sep="\t", header=False, index=False, float_format="%.6g")
        os.replace(tmp, os.path.join(parts, f"chunk_{q:06d}.tsv.gz"))
        with open(done_path, "a") as fh:
            fh.write(f"{q}\n")
        if tracks is not None:
            m0 = np.flatnonzero(kind0[idx])
            if len(m0):
                m0 = rng.choice(m0, size=min(len(m0), 5000), replace=False)
                ix = p0[idx][m0][:, None] + OFFS[None, :]
                for s, T, mine in (("plus", tracks[0], rp), ("minus", tracks[1], rm)):
                    tr = np.asarray(T[ix.ravel()], dtype=np.float32).reshape(ix.shape).max(axis=1)
                    chk[s].append(np.abs(mine[m0] - tr) / (np.abs(tr) + 1e-4))
        n_run += len(idx)
        if (ci + 1) % 10 == 0 or ci == len(todo) - 1:
            el = time.time() - t1
            print(f"  chunk {ci + 1}/{len(todo)}: {n_run:,} variants; {el / 60:.1f} min, "
                  f"{1000 * el / max(1, n_run):.2f} ms per variant", flush=True)

    if tracks is not None and chk["plus"]:
        print(f"\n[{c}] check, reference readout vs 16c track (tile_kind 0, relative difference):")
        for s in ("plus", "minus"):
            x = np.concatenate(chk[s])
            print(f"  {s:5s} n={len(x):,}  median {np.median(x):.2e}  99% {np.quantile(x, 0.99):.2e}  "
                  f"99.9% {np.quantile(x, 0.999):.2e}  max {x.max():.2e}  within 2e-3: {100 * (x <= 2e-3).mean():.2f}%")
        worst = max(np.quantile(np.concatenate(chk[s]), 0.99) for s in ("plus", "minus"))
        print(f"  CHECK {'OK' if worst <= 5e-3 else 'FAILED: the reference does not reproduce the 16c track'}")

    done = {int(x) for x in open(done_path).read().split()}
    missing = [q for q in chunk_ids if q not in done]
    if missing:
        print(f"\n[{c}] {len(missing):,} chunks still to do; resubmit to continue  [{(time.time() - t0) / 60:.1f} min]")
        return
    tmp = out + ".tmp"
    with open(tmp, "wb") as fo:
        fo.write(gzip.compress(("\t".join(COLS_OUT) + "\n").encode()))
        for q in chunk_ids:
            with open(os.path.join(parts, f"chunk_{q:06d}.tsv.gz"), "rb") as fi:
                shutil.copyfileobj(fi, fo)
    os.replace(tmp, out)
    d = pd.read_csv(out, sep="\t", usecols=["pos", "ref_plus", "ref_minus", "alt_plus", "alt_minus", "tile_kind"])
    if len(d) != n or not np.array_equal(d["pos"].to_numpy(), pos):
        sys.exit(f"STOP: assembled file has {len(d):,} rows or a different order (expected {n:,})")
    shutil.rmtree(parts)
    os.remove(done_path)
    dp = (d["alt_plus"] - d["ref_plus"]).abs().to_numpy()
    dm = (d["alt_minus"] - d["ref_minus"]).abs().to_numpy()
    dd = np.maximum(dp, dm)
    print(f"\n[{c}] written {out}: {len(d):,} rows, same order as the input")
    print(f"  reference readout quantiles 50/90/99/99.9%: plus " +
          ", ".join(f"{np.quantile(d['ref_plus'], q):.3g}" for q in (0.5, 0.9, 0.99, 0.999)) + "; minus " +
          ", ".join(f"{np.quantile(d['ref_minus'], q):.3g}" for q in (0.5, 0.9, 0.99, 0.999)))
    print("  largest change on either strand, quantiles 50/90/99/99.9%: " +
          ", ".join(f"{np.quantile(dd, q):.3g}" for q in (0.5, 0.9, 0.99, 0.999)))
    print(f"  done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
