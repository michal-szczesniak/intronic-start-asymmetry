#!/usr/bin/env python3
"""Step 11c: can many gnomAD variants be scored in one Puffin-D pass, and how fast?

Step 3 of the gnomAD test needs the effect of each of ~20 million intronic SNVs on
predicted initiation, on both strands. One forward pass per variant (x2 for the
reverse complement) would take days, so the plan is to put many variants into one
100 kb window at once, far enough apart that they do not affect each other's
readout, and to run windows in batches on the GPU. This script measures whether
that is safe and what it costs, on real gnomAD variants of one chromosome:

  1. Equivalence: the batched predictor used here gives the same output as the
     model's own predict() (forward and reverse-complement averaged, as in step 08).
  2. Local effect of single variants: for each variant alone, the change of the
     predicted maximum within +/-HW nt on each strand (alt - ref), and how often a
     variant creates a site (ref below a threshold, alt at or above it).
  3. Cross-talk with distance: from the same single-variant runs, the change of the
     readout centred d nt away from the variant. This sets the minimum spacing.
  4. Multiplexing: one window carrying all variants spaced >= D nt (D = 500, 1000,
     2000); each variant's change is compared with its single-variant value.
  5. Precision: bfloat16 against float32 on the single-variant effects.
  6. Speed: milliseconds per 100 kb pass for several batch sizes, float32 and bf16.

Thresholds for "site" are those of steps 09b/10 (top 1e-3 and 3e-4 of the real
predictions on chr8+chr9): 0.0439 and 0.1432 (override with --thresholds).

Usage:
  python3 11c_mux_test.py --genome /path/hg38.fa --puffin-dir /path/tools/puffin \
      --gnomad /scratch/.../gnomad_v4.1/chr8.snv_pass.tsv.gz --chrom chr8 \
      --out-prefix results/11c_mux_chr8 [--n-windows 20] [--seed 1] [--no-cuda]
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

WINDOW = 100000
CENTRAL = (25000, 75000)
HW = 50
TRACK_PRO_CAP_PLUS = 4
TRACK_PRO_CAP_MINUS = 5
DISTS = (100, 200, 300, 500, 1000, 2000, 5000)
SPACINGS = (500, 1000, 2000)
MAX_PER_WINDOW = 100
COMP = str.maketrans("ACGTN", "TGCAN")
LUT = np.full((256, 4), 0.25, dtype=np.float32)   # unknown bases -> 0.25 each
for i, b in enumerate("ACGT"):
    LUT[ord(b)] = 0
    LUT[ord(b), i] = 1
    LUT[ord(b.lower())] = LUT[ord(b)]


def revcomp(s):
    return s.translate(COMP)[::-1]


class Predictor:
    """Forward + reverse-complement averaged PRO-cap tracks for a list of sequences.
    Uses batched torch forward when available; falls back to model.predict()."""

    def __init__(self, model, use_cuda):
        self.model, self.use_cuda = model, use_cuda
        try:
            import torch
            self.torch = torch if hasattr(model, "forward") and isinstance(model, torch.nn.Module) else None
        except ImportError:
            self.torch = None

    def _raw(self, seqs, batch, bf16):
        if self.torch is None:
            return [np.asarray(self.model.predict(s), dtype=np.float32).reshape(-1, len(s))[[TRACK_PRO_CAP_PLUS,
                    TRACK_PRO_CAP_MINUS]] for s in seqs]
        torch, out = self.torch, []
        for i in range(0, len(seqs), batch):
            chunk = seqs[i:i + batch]
            x = np.stack([LUT[np.frombuffer(s.encode(), dtype=np.uint8)].T for s in chunk])
            xt = torch.from_numpy(x)
            if self.use_cuda:
                xt = xt.cuda()
            with torch.no_grad():
                if bf16:
                    with torch.autocast("cuda", dtype=torch.bfloat16):
                        y = self.model(xt)
                else:
                    y = self.model(xt)
            y = y[:, [TRACK_PRO_CAP_PLUS, TRACK_PRO_CAP_MINUS]].float().cpu().numpy()
            out.extend(list(y))
        return out

    def rc_avg(self, seqs, batch=16, bf16=False):
        raw = self._raw(list(seqs) + [revcomp(s) for s in seqs], batch, bf16)
        n = len(seqs)
        res = []
        for k in range(n):
            f, r = raw[k], raw[n + k]
            res.append(((f[0] + r[1][::-1]) / 2, (f[1] + r[0][::-1]) / 2))
        return res


def readout(track, x):
    return float(track[x - HW:x + HW + 1].max())


def greedy_spaced(pos, d, cap):
    sel, last = [], -10 ** 9
    for i, p in enumerate(pos):
        if p - last >= d:
            sel.append(i)
            last = p
            if len(sel) >= cap:
                break
    return sel


def mutate(seq, offs, alts):
    b = bytearray(seq.encode())
    for o, a in zip(offs, alts):
        b[o] = ord(a)
    return b.decode()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--puffin-dir", required=True)
    ap.add_argument("--gnomad", required=True, help="table from step 11b: chrom pos ref alt AC AN")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-windows", type=int, default=20)
    ap.add_argument("--thresholds", default="0.0439,0.1432")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-cuda", action="store_true")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    thr = [float(x) for x in a.thresholds.split(",")]
    genome, gnomad, prefix, puffin_dir = (os.path.abspath(x) for x in (a.genome, a.gnomad, a.out_prefix, a.puffin_dir))

    from pyfaidx import Fasta
    fa = Fasta(genome, as_raw=True, sequence_always_upper=False)
    clen = len(fa[a.chrom])
    cwd = os.getcwd()
    os.chdir(puffin_dir)  # the model resolves resources/ relative to the working directory
    sys.path.insert(0, puffin_dir)
    from puffin_D import Puffin_D
    model = Puffin_D(use_cuda=not a.no_cuda)
    os.chdir(cwd)
    P = Predictor(model, not a.no_cuda)
    print(f"batched torch path: {P.torch is not None}", flush=True)

    t0 = time.time()
    v = pd.read_csv(gnomad, sep="\t", header=None, usecols=[1, 2, 3], names=["pos", "ref", "alt"],
                    dtype={"pos": np.int64, "ref": "category", "alt": "category"})
    v = v[(v["ref"].str.len() == 1) & (v["alt"].str.len() == 1)]
    vpos = v["pos"].to_numpy() - 1
    print(f"{len(v):,} SNVs on {a.chrom} loaded in {time.time() - t0:.0f} s", flush=True)

    # windows with few N and some variants in the central part
    wins = []
    while len(wins) < a.n_windows:
        s = int(rng.integers(0, clen - WINDOW))
        seq = str(fa[a.chrom][s:s + WINDOW]).upper()
        if seq.count("N") > 0.01 * WINDOW:
            continue
        lo, hi = np.searchsorted(vpos, [s + CENTRAL[0] + HW, s + CENTRAL[1] - HW])
        if hi - lo < 50:
            continue
        wins.append((s, seq, lo, hi))

    # 1. equivalence with the model's own predict()
    print("\n== 1. equivalence: batched predictor vs model.predict(), forward+RC averaged ==")
    for s, seq, _, _ in wins[:2]:
        mine = P.rc_avg([seq])[0]
        own_f = np.asarray(model.predict(seq), dtype=np.float32).reshape(-1, WINDOW)
        own_r = np.asarray(model.predict(revcomp(seq)), dtype=np.float32).reshape(-1, WINDOW)
        own = ((own_f[4] + own_r[5][::-1]) / 2, (own_f[5] + own_r[4][::-1]) / 2)
        d = max(float(np.abs(mine[0] - own[0]).max()), float(np.abs(mine[1] - own[1]).max()))
        print(f"  window {a.chrom}:{s}: max |difference| = {d:.2e}")

    rows, mux_rows, dist_rows, bf_rows = [], [], [], []
    n_ref_ok = n_ref = 0
    t1 = time.time()
    for wi, (s, seq, lo, hi) in enumerate(wins):
        sub = v.iloc[lo:hi]
        offs_all = sub["pos"].to_numpy() - 1 - s
        refb = np.frombuffer(seq.encode(), dtype=np.uint8)[offs_all]
        match = refb == np.frombuffer("".join(sub["ref"]).upper().encode(), dtype=np.uint8)
        n_ref += len(match)
        n_ref_ok += int(match.sum())
        offs_all, alts_all = offs_all[match], sub["alt"].to_numpy()[match]
        order = rng.permutation(len(offs_all))
        offs_all, alts_all = offs_all[order], alts_all[order]
        srt = np.argsort(offs_all, kind="stable")
        offs_all, alts_all = offs_all[srt], alts_all[srt]
        s500 = greedy_spaced(offs_all, SPACINGS[0], MAX_PER_WINDOW)
        offs, alts = offs_all[s500], alts_all[s500]

        ref_p, ref_m = P.rc_avg([seq])[0]
        singles = P.rc_avg([mutate(seq, [o], [x]) for o, x in zip(offs, alts)])
        for k, (o, x) in enumerate(zip(offs, alts)):
            ap_, am_ = singles[k]
            r = {"window": wi, "offset": int(o)}
            for name, rt, at in (("plus", ref_p, ap_), ("minus", ref_m, am_)):
                r[f"ref_{name}"] = readout(rt, o)
                r[f"alt_{name}"] = readout(at, o)
                r[f"d_{name}"] = r[f"alt_{name}"] - r[f"ref_{name}"]
            rows.append(r)
            for d in DISTS:
                for sign in (-1, 1):
                    c = o + sign * d
                    if HW <= c < WINDOW - HW:
                        for name, rt, at in (("plus", ref_p, ap_), ("minus", ref_m, am_)):
                            dist_rows.append({"d": d, "abs_change": abs(readout(at, c) - readout(rt, c))})
        single = pd.DataFrame(rows[-len(offs):])

        # 4. multiplexing at several spacings
        idx = list(range(len(offs)))
        for D in SPACINGS:
            sel = greedy_spaced(offs, D, MAX_PER_WINDOW) if D > SPACINGS[0] else idx
            mp, mm = P.rc_avg([mutate(seq, offs[sel], alts[sel])])[0]
            for k in sel:
                o = offs[k]
                for name, rt, mt in (("plus", ref_p, mp), ("minus", ref_m, mm)):
                    mux_rows.append({"D": D, "window": wi, "strand": name,
                                     "d_single": single.iloc[k][f"d_{name}"],
                                     "d_mux": readout(mt, o) - readout(rt, o),
                                     "ref": readout(rt, o), "alt_single": single.iloc[k][f"alt_{name}"],
                                     "alt_mux": readout(mt, o)})

        # 5. bf16 vs fp32 on the first windows
        if wi < 3 and P.torch is not None and not a.no_cuda:
            rb = P.rc_avg([seq], bf16=True)[0]
            sb = P.rc_avg([mutate(seq, [o], [x]) for o, x in zip(offs, alts)], bf16=True)
            for k, o in enumerate(offs):
                for j, name in enumerate(("plus", "minus")):
                    bf_rows.append({"d_fp32": single.iloc[k][f"d_{name}"],
                                    "d_bf16": readout(sb[k][j], o) - readout(rb[j], o)})
        print(f"  window {wi + 1}/{len(wins)} {a.chrom}:{s}: {len(offs)} variants, {time.time() - t1:.0f} s", flush=True)

    R = pd.DataFrame(rows)
    M = pd.DataFrame(mux_rows)
    Dd = pd.DataFrame(dist_rows)
    R.to_csv(prefix + "_singles.tsv", sep="\t", index=False)
    M.to_csv(prefix + "_mux.tsv", sep="\t", index=False)

    print(f"\nreference allele matches the genome: {n_ref_ok}/{n_ref} ({100 * n_ref_ok / max(1, n_ref):.2f}%)")
    print(f"\n== 2. local effect of single variants (n = {len(R)} variants, both strands) ==")
    dd = np.abs(np.r_[R["d_plus"], R["d_minus"]])
    print("  |change| quantiles 50/90/99/99.9%: " + ", ".join(f"{np.quantile(dd, q):.4f}" for q in (.5, .9, .99, .999)))
    for t in (0.001, 0.005, 0.02):
        print(f"  |change| > {t}: {100 * (dd > t).mean():.2f}% of variant-strand pairs")
    for T in thr:
        cr = np.r_[(R["ref_plus"] < T) & (R["alt_plus"] >= T), (R["ref_minus"] < T) & (R["alt_minus"] >= T)]
        de = np.r_[(R["ref_plus"] >= T) & (R["alt_plus"] < T), (R["ref_minus"] >= T) & (R["alt_minus"] < T)]
        print(f"  threshold {T}: creates a site {100 * cr.mean():.3f}%, destroys one {100 * de.mean():.3f}% "
              f"of variant-strand pairs")

    print("\n== 3. cross-talk: change of the readout centred d nt away from a single variant ==")
    loc_med = float(np.median(dd[dd > 0])) if (dd > 0).any() else np.nan
    for d, g in Dd.groupby("d"):
        x = g["abs_change"].to_numpy()
        print(f"  d = {d:5d} nt: median {np.median(x):.2e}, 99th pct {np.quantile(x, .99):.2e}, "
              f"max {x.max():.2e}, > 0.001 in {100 * (x > 0.001).mean():.2f}%")
    print(f"  (for scale: median local |change| among nonzero = {loc_med:.2e}; site thresholds {thr})")

    print("\n== 4. multiplexed windows vs single-variant runs ==")
    for D, g in M.groupby("D"):
        err = np.abs(g["d_mux"] - g["d_single"])
        r = pearsonr(g["d_single"], g["d_mux"])[0] if g["d_single"].std() > 0 else np.nan
        dis = {T: int((((g["ref"] < T) & (g["alt_single"] >= T)) != ((g["ref"] < T) & (g["alt_mux"] >= T))).sum())
               for T in thr}
        print(f"  spacing >= {D:4d} nt: {len(g) // 2} variants ({len(g) // 2 / len(wins):.0f} per window); "
              f"Pearson r {r:.5f}; |error| 99th pct {np.quantile(err, .99):.2e}, max {err.max():.2e}; "
              f"site-creation calls that differ: " + ", ".join(f"{T}: {n}" for T, n in dis.items()))

    if bf_rows:
        B = pd.DataFrame(bf_rows)
        err = np.abs(B["d_bf16"] - B["d_fp32"])
        print(f"\n== 5. bfloat16 vs float32 (single-variant changes, n = {len(B)}) ==")
        print(f"  Pearson r {pearsonr(B['d_fp32'], B['d_bf16'])[0]:.5f}; |error| median {np.median(err):.2e}, "
              f"99th pct {np.quantile(err, .99):.2e}, max {err.max():.2e}")

    print("\n== 6. speed (ms per 100 kb forward pass) ==")
    seq = wins[0][1]
    for bf16 in ((False, True) if P.torch is not None and not a.no_cuda else (False,)):
        for bs in (1, 4, 8, 16, 32):
            try:
                P._raw([seq] * bs, bs, bf16)  # warm-up
                if P.torch is not None and not a.no_cuda:
                    P.torch.cuda.synchronize()
                t = time.time()
                n = max(32, bs * 4)
                P._raw([seq] * n, bs, bf16)
                if P.torch is not None and not a.no_cuda:
                    P.torch.cuda.synchronize()
                print(f"  {'bf16' if bf16 else 'fp32'} batch {bs:2d}: {1000 * (time.time() - t) / n:.1f} ms per pass")
            except RuntimeError as e:
                print(f"  {'bf16' if bf16 else 'fp32'} batch {bs:2d}: failed ({str(e)[:80]})")
                break
    print(f"\nwritten: {prefix}_singles.tsv, {prefix}_mux.tsv")


if __name__ == "__main__":
    main()
