#!/usr/bin/env python3
"""Step 16b: score the iSTART pairs with ProCapNet, each pair by a model that never saw its
chromosome.

Puffin-D (steps 07/07c) has only ~150-200 pairs on its held-out chromosomes, so how much
of the observed sense weakness sequence explains rests on chromosomes the model was
trained on. ProCapNet (K562) was trained 7 times, each time with a different group of
chromosomes held out as the test set (cv_params.fold_k.json). Here every pair is scored
by the model whose test set contains its chromosome ("held-out" prediction). For
comparison the same pair is also scored by the models that were TRAINED on its
chromosome (average over them: "train-fold" prediction). If the train-fold predictions
reproduce the observed asymmetry much better than the held-out ones, the models
memorised measured signal; the held-out predictions are the clean estimate.

Model: BPNet architecture as in ProCapNet (input 2114 nt, one-hot A,C,G,T, N = zeros;
output a profile over 1000 nt x 2 strands [plus, minus] whose softmax runs over both
strands together, and the log of the total count). Each window is scored as given and
as its reverse complement, mapped back and averaged (as ProCapNet's own evaluation):
probabilities averaged, log counts averaged. Predicted counts per base = probability x
(exp(log count) - 1). The window is centred between the two arms; pairs whose arms do
not fit in the central 1000 nt with the readout margin are marked.

Readouts per pair (s = sense arm on its strand, a = antisense arm on its strand):
  pred_DI       log2((counts within +/-hw of s + 1) / (same for a + 1))   [main; as observed DI]
  pred_DI_sum   the same with pseudocount 0.01                             [no pseudocount effect]
  pred_DI_point log2((max within +/-hw at s + 0.01) / (same at a + 0.01))  [peak height]
  pred_DI_wide  as pred_DI within +/-250 nt                                [wider window]
The column names match step 07c, so 07b_calibrate.py runs unchanged on the output.
Puffin-D's columns are kept with the prefix puffin_.

Checks printed: offset of the predicted maximum from each arm and which strand wins at
plus / minus arms (coordinates and strand order), forward vs reverse-complement
consistency, and Spearman of observed with predicted DI (held-out, train-fold, Puffin-D).

Usage:
  python3 16b_score_pairs_procapnet.py --pairs-dir ../07_puffin_pairs/results --lines K562,HCT116 \
      --genome hg38.fa --models-dir $SCRATCH/models/procapnet_K562/ENCFF976FHE \
      --cv-dir $SCRATCH/models/procapnet_K562/ENCFF033LJA --outdir results \
      [--signal-window 50] [--batch 64] [--no-train-folds] [--device auto] [--threads 0]
Writes <outdir>/procapnet_heldout_pairs_<line>.tsv and procapnet_trainfold_pairs_<line>.tsv.
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

IN_W, OUT_W = 2114, 1000
TRIM = (IN_W - OUT_W) // 2
WIDE = 250
OFFSET_R = 200
LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i
METRICS = ("pred_DI", "pred_DI_point", "pred_DI_wide", "pred_DI_sum")


def load_folds(models_dir, cv_dir):
    folds = []
    for p in sorted(glob.glob(os.path.join(cv_dir, "fold_*", "cv_params.fold_*.json"))):
        k = int(os.path.basename(p).split(".")[1].split("_")[1])
        with open(p) as fh:
            cv = json.load(fh)
        m = glob.glob(os.path.join(models_dir, f"fold_{k}", "*.state_dict.torch"))
        if len(m) != 1:
            sys.exit(f"STOP: expected one state_dict for fold {k} in {models_dir}, found {m}")
        folds.append({"k": k, "test": set(cv["test"]), "valid": set(cv["valid"]), "train": set(cv["train"]), "path": m[0]})
    if not folds:
        sys.exit(f"STOP: no cv_params files under {cv_dir}")
    allt = [c for f in folds for c in f["test"]]
    if len(allt) != len(set(allt)):
        sys.exit("STOP: a chromosome is in more than one test set")
    return folds


def torch_predictor(path, device, threads):
    """Load one ProCapNet state_dict; return f(X: (B,4,2114) float32) -> (logprob (B,2,1000), logcount (B,))."""
    import torch
    if threads:
        torch.set_num_threads(threads)

    class ProCapNet(torch.nn.Module):
        """Same layers and names as ProCapNet's BPNet_strand_merged_umap.Model."""

        def __init__(self, n_filters=512, n_layers=8, n_outputs=2, trimming=TRIM):
            super().__init__()
            self.n_layers, self.trimming, self.deconv_kernel_size = n_layers, trimming, 75
            self.iconv = torch.nn.Conv1d(4, n_filters, kernel_size=21, padding=10)
            self.rconvs = torch.nn.ModuleList([torch.nn.Conv1d(n_filters, n_filters, kernel_size=3, padding=2 ** i,
                                                               dilation=2 ** i) for i in range(1, n_layers + 1)])
            self.fconv = torch.nn.Conv1d(n_filters, n_outputs, kernel_size=self.deconv_kernel_size)
            self.linear = torch.nn.Linear(n_filters, 1)

        def forward(self, X):
            start, end = self.trimming, X.shape[2] - self.trimming
            X = torch.relu(self.iconv(X))
            for i in range(self.n_layers):
                X = X + torch.relu(self.rconvs[i](X))
            X = X[:, :, start - self.deconv_kernel_size // 2: end + self.deconv_kernel_size // 2]
            y_profile = self.fconv(X)
            y_counts = self.linear(torch.mean(X, dim=2)).reshape(X.shape[0], 1)
            return y_profile, y_counts

    try:
        sd = torch.load(path, map_location="cpu", weights_only=True)
    except TypeError:
        sd = torch.load(path, map_location="cpu")
    model = ProCapNet()
    model.load_state_dict(sd, strict=True)
    model.to(device).eval()

    def predict(X):
        with torch.no_grad():
            prof, cnt = model(torch.from_numpy(X).to(device))
            B = prof.shape[0]
            lp = torch.log_softmax(prof.reshape(B, -1), dim=-1).reshape(B, 2, -1)
        return lp.float().cpu().numpy(), cnt[:, 0].float().cpu().numpy()
    return predict


def fake_predictor(path, motif="TATAAA", bias=0.0):
    """Test stand-in: a start 3 nt after every motif on the plus strand and its mirror image on the
    minus strand; a fold-specific multiplier on the plus strand (bias) mimics memorisation."""
    mot = np.array([LUT[ord(c)] for c in motif])
    rc = 3 - mot[::-1]
    k = len(mot)

    def predict(X):
        B = X.shape[0]
        idx = np.where(X.sum(1) > 0, X.argmax(1), -1)
        out = np.full((B, 2, OUT_W), -6.0, dtype=np.float32)
        for b in range(B):
            s = idx[b]
            for j in range(TRIM - 10, TRIM + OUT_W + 10):
                if np.array_equal(s[j:j + k], mot) and 0 <= j + 3 - TRIM < OUT_W:
                    out[b, 0, j + 3 - TRIM] = 3.0 + bias
                if np.array_equal(s[j - k + 1:j + 1], rc) and 0 <= j - 3 - TRIM < OUT_W:
                    out[b, 1, j - 3 - TRIM] = 3.0
        flat = out.reshape(B, -1)
        lp = flat - np.log(np.exp(flat).sum(1, keepdims=True))
        return lp.reshape(B, 2, OUT_W), np.full(B, np.log(101.0), dtype=np.float32)
    return predict


def merged_counts(predict, X, batch):
    """RC-averaged predicted counts (B, 2, 1000) and the forward / RC-only versions."""
    outs = {"fwd": [], "rc": [], "avg": []}
    for i in range(0, len(X), batch):
        x = X[i:i + batch]
        lp_f, lc_f = predict(np.ascontiguousarray(x))
        lp_r, lc_r = predict(np.ascontiguousarray(x[:, ::-1, ::-1]))
        lp_r = lp_r[:, ::-1, ::-1]                                    # strands swapped, positions reversed
        p_f, p_r = np.exp(lp_f), np.exp(lp_r)
        tot_f, tot_r = np.expm1(lc_f).clip(0), np.expm1(lc_r).clip(0)
        tot_a = np.expm1((lc_f + lc_r) / 2).clip(0)
        outs["fwd"].append(p_f * tot_f[:, None, None])
        outs["rc"].append(p_r * tot_r[:, None, None])
        outs["avg"].append(((p_f + p_r) / 2) * tot_a[:, None, None])
    return {k: np.concatenate(v) for k, v in outs.items()}


def readouts(C, strand_sense, i_s, i_a, hw):
    """C: (n, 2, 1000) predicted counts; returns dict of arrays."""
    n = len(C)
    sense_tr = np.where(strand_sense == "+", 0, 1)
    anti_tr = 1 - sense_tr
    r = np.arange(n)

    def win(tr, i, h):
        lo, hi = np.clip(i - h, 0, OUT_W), np.clip(i + h + 1, 0, OUT_W)
        cs = np.concatenate([np.zeros((n, 1)), np.cumsum(C[r, tr], axis=1)], axis=1)
        return cs[r, hi] - cs[r, lo]

    def wmax(tr, i, h):
        return np.array([C[j, tr[j], max(0, i[j] - h):i[j] + h + 1].max() for j in range(n)])

    su_s, su_a = win(sense_tr, i_s, hw), win(anti_tr, i_a, hw)
    wd_s, wd_a = win(sense_tr, i_s, WIDE), win(anti_tr, i_a, WIDE)
    mx_s, mx_a = wmax(sense_tr, i_s, hw), wmax(anti_tr, i_a, hw)
    return {"pred_sense": su_s, "pred_antisense": su_a,
            "pred_DI": np.log2((su_s + 1) / (su_a + 1)), "pred_DI_sum": np.log2((su_s + 0.01) / (su_a + 0.01)),
            "pred_DI_point": np.log2((mx_s + 0.01) / (mx_a + 0.01)), "pred_DI_wide": np.log2((wd_s + 1) / (wd_a + 1))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-dir", required=True, help="directory with puffin_rc_pairs_<line>.tsv (step 07c)")
    ap.add_argument("--lines", default="K562,HCT116")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--models-dir", required=True)
    ap.add_argument("--cv-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--signal-window", type=int, default=50)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--no-train-folds", action="store_true", help="held-out predictions only (about 6x faster)")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--fake-model", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args()
    hw = a.signal_window
    os.makedirs(a.outdir, exist_ok=True)
    folds = load_folds(a.models_dir, a.cv_dir)
    print("folds (test chromosomes): " + "; ".join(f"{f['k']}: {','.join(sorted(f['test']))}" for f in folds))
    if a.fake_model:
        device = "cpu"
        get_pred = lambda f: fake_predictor(f["path"], bias=0.0 if f.get("role") == "test" else 0.7)
    else:
        import torch
        device = ("cuda" if torch.cuda.is_available() else "cpu") if a.device == "auto" else a.device
        get_pred = lambda f: torch_predictor(f["path"], device, a.threads)
        print(f"torch {torch.__version__}, device {device}" + (f" ({torch.cuda.get_device_name(0)})" if device == "cuda" else
                                                               f", {torch.get_num_threads()} threads"))
    from pyfaidx import Fasta
    fa = Fasta(a.genome, as_raw=True, sequence_always_upper=False)

    for L in [x for x in a.lines.split(",") if x]:
        t0 = time.time()
        f_in = os.path.join(a.pairs_dir, f"puffin_rc_pairs_{L}.tsv")
        if not os.path.isfile(f_in):
            sys.exit(f"STOP: missing {f_in}")
        d = pd.read_csv(f_in, sep="\t", keep_default_na=False, na_values=["", "NA", "nan", "NaN"],
                        dtype={"host": str, "klasa": str, "strand_sense": str, "chrom": str})
        ren = {c: "puffin_" + c for c in d.columns if c.startswith("pred_") or c in ("window_start", "rc_lag_sense",
               "calib_offset_sense", "calib_offset_antisense", "track_plus_wins_at_plus_arm", "track_minus_wins_at_minus_arm")}
        d = d.rename(columns=ren)
        n = len(d)
        center = ((d["pos_sense"].to_numpy(np.int64) + d["pos_antisense"].to_numpy(np.int64)) // 2)
        in_start = center - IN_W // 2
        out_start = in_start + TRIM
        i_s = d["pos_sense"].to_numpy(np.int64) - out_start
        i_a = d["pos_antisense"].to_numpy(np.int64) - out_start
        status = np.full(n, "ok", dtype=object)
        fold_of = {c: f["k"] for f in folds for c in f["test"]}
        chroms = d["chrom"].astype(str).to_numpy()
        status[~np.isin(chroms, list(fold_of))] = "chrom_not_in_any_fold"
        bad = (i_s < hw) | (i_s >= OUT_W - hw) | (i_a < hw) | (i_a >= OUT_W - hw)
        status[(status == "ok") & bad] = "arm_outside_window"
        X = np.zeros((n, 4, IN_W), dtype=np.float32)
        nfrac = np.ones(n)
        for c in np.unique(chroms[status == "ok"]):
            if c not in fa:
                status[(chroms == c) & (status == "ok")] = "no_such_contig"
                continue
            codes = LUT[np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)]
            for j in np.flatnonzero((chroms == c) & (status == "ok")):
                lo, hi = in_start[j], in_start[j] + IN_W
                seg = np.full(IN_W, -1, dtype=np.int8)
                a_, b_ = max(lo, 0), min(hi, len(codes))
                seg[a_ - lo:b_ - lo] = codes[a_:b_]
                ok = seg >= 0
                X[j, seg[ok], np.flatnonzero(ok)] = 1.0
                nfrac[j] = 1 - ok.mean()
        status[(status == "ok") & (nfrac > 0.5)] = "too_many_N"
        # local GC around the sense arm (as 07c) where 07c left it empty, so 07b keeps every scored pair
        if "gc_local" not in d.columns:
            d["gc_local"] = np.nan
        miss = np.flatnonzero(~np.isfinite(pd.to_numeric(d["gc_local"], errors="coerce").to_numpy(float)) & (status == "ok"))
        for j in miss:
            ps = int(d["pos_sense"].iat[j])
            loc = str(fa[chroms[j]][max(0, ps - 500):ps + 500]).upper()
            d.loc[d.index[j], "gc_local"] = (loc.count("G") + loc.count("C")) / max(1, len(loc))
        ok_idx = np.flatnonzero(status == "ok")
        print(f"\n######## {L}: {n:,} pairs; scoring {len(ok_idx):,}; " +
              ", ".join(f"{k} {v:,}" for k, v in pd.Series(status).value_counts().items() if k != "ok"))

        Ch = {k: np.zeros((n, 2, OUT_W), np.float32) for k in ("fwd", "rc", "avg")}
        Ct = np.zeros((n, 2, OUT_W), np.float32)
        nt = np.zeros(n)
        fold_used = np.full(n, -1)
        for f in folds:
            sel_h = ok_idx[np.isin(chroms[ok_idx], list(f["test"]))]
            sel_t = np.array([], dtype=int) if a.no_train_folds else ok_idx[np.isin(chroms[ok_idx], list(f["train"]))]
            if not len(sel_h) and not len(sel_t):
                continue
            t1 = time.time()
            if len(sel_h):
                f["role"] = "test"
                out = merged_counts(get_pred(f), X[sel_h], a.batch)
                for k in Ch:
                    Ch[k][sel_h] = out[k]
                fold_used[sel_h] = f["k"]
            if len(sel_t):
                f["role"] = "train"
                out = merged_counts(get_pred(f), X[sel_t], a.batch)
                Ct[sel_t] += out["avg"]
                nt[sel_t] += 1
            print(f"  fold {f['k']}: {len(sel_h):,} held-out, {len(sel_t):,} train-fold pairs in {time.time() - t1:.0f} s",
                  flush=True)
        Ct[nt > 0] /= nt[nt > 0][:, None, None]

        ss = d["strand_sense"].astype(str).to_numpy()
        res = {}
        for tag, C in (("heldout", Ch["avg"]), ("trainfold", Ct)):
            if tag == "trainfold" and a.no_train_folds:
                continue
            ro = readouts(C[ok_idx], ss[ok_idx], i_s[ok_idx], i_a[ok_idx], hw)
            D = d.copy()
            st = status.copy()
            if tag == "trainfold":
                st[(st == "ok") & (nt == 0)] = "no_train_fold"
            D["pred_status"] = st
            for k, v in ro.items():
                col = np.full(n, np.nan)
                col[ok_idx] = v
                D[k] = col
            D["pc_fold_heldout"] = fold_used
            D["pc_n_train_folds"] = nt.astype(int)
            D["pc_total_counts"] = np.where(status == "ok", C.sum(axis=(1, 2)), np.nan)
            if tag == "heldout":
                for side, C2 in (("fwd", Ch["fwd"]), ("rc", Ch["rc"])):
                    col = np.full(n, np.nan)
                    col[ok_idx] = readouts(C2[ok_idx], ss[ok_idx], i_s[ok_idx], i_a[ok_idx], hw)["pred_DI"]
                    D[f"pred_DI_{side}"] = col
                # coordinates and strand order: offset of the predicted maximum from each arm, winning strand
                for arm, ii in (("sense", i_s), ("antisense", i_a)):
                    tr = np.where((ss == "+") == (arm == "sense"), 0, 1)
                    off = np.full(n, np.nan)
                    for j in ok_idx:
                        lo = max(0, ii[j] - OFFSET_R)
                        off[j] = int(np.argmax(C[j, tr[j], lo:ii[j] + OFFSET_R + 1])) + lo - ii[j]
                    D[f"pc_offset_{arm}"] = off
                pp = np.where(ss == "+", i_s, i_a)
                pm = np.where(ss == "+", i_a, i_s)
                wins_p, wins_m = np.full(n, np.nan), np.full(n, np.nan)
                for j in ok_idx:
                    sl_p = slice(max(0, pp[j] - hw), pp[j] + hw + 1)
                    sl_m = slice(max(0, pm[j] - hw), pm[j] + hw + 1)
                    wins_p[j] = float(C[j, 0, sl_p].sum() > C[j, 1, sl_p].sum())
                    wins_m[j] = float(C[j, 1, sl_m].sum() > C[j, 0, sl_m].sum())
                D["pc_plus_wins_at_plus_arm"], D["pc_minus_wins_at_minus_arm"] = wins_p, wins_m
            out = os.path.join(a.outdir, f"procapnet_{tag}_pairs_{L}.tsv")
            D.to_csv(out, sep="\t", index=False)
            res[tag] = D
            print(f"  written {out}")

        H = res["heldout"]
        okh = H[H["pred_status"] == "ok"]
        print(f"\n  checks [{L}]")
        print(f"    offset of predicted maximum from the arm (within +/-{OFFSET_R} nt): sense median "
              f"{okh['pc_offset_sense'].median():.0f}, |offset| <= 5 nt for {100 * (okh['pc_offset_sense'].abs() <= 5).mean():.0f}%; "
              f"antisense median {okh['pc_offset_antisense'].median():.0f}, "
              f"|offset| <= 5 for {100 * (okh['pc_offset_antisense'].abs() <= 5).mean():.0f}%")
        print(f"    plus strand stronger at plus arms {100 * okh['pc_plus_wins_at_plus_arm'].mean():.0f}%, "
              f"minus strand stronger at minus arms {100 * okh['pc_minus_wins_at_minus_arm'].mean():.0f}%")
        print(f"    Spearman(DI forward pass, DI reverse-complement pass) = "
              f"{spearmanr(okh['pred_DI_fwd'], okh['pred_DI_rc'])[0]:.3f}")
        intr = okh["klasa"] == "intron genu"
        print("    Spearman(observed DI, predicted DI)       intergenic, both orientations | intronic")
        cmp_ = [("ProCapNet held-out", okh["pred_DI"])]
        if "trainfold" in res:
            T = res["trainfold"].loc[okh.index]
            cmp_.append(("ProCapNet train-fold", T["pred_DI"]))
        if "puffin_pred_DI" in okh.columns:
            cmp_.append(("Puffin-D (07c)", okh["puffin_pred_DI"]))
        for name, v in cmp_:
            m_ok = np.isfinite(v.to_numpy(float))
            x, y = okh["DI"].to_numpy(float), v.to_numpy(float)
            ci_ = (~intr.to_numpy()) & m_ok
            r_ig = spearmanr(np.r_[x[ci_], -x[ci_]], np.r_[y[ci_], -y[ci_]])[0]
            r_in = spearmanr(x[intr.to_numpy() & m_ok], y[intr.to_numpy() & m_ok])[0]
            print(f"      {name:22s} {r_ig:+.3f} | {r_in:+.3f}")
        print(f"  [{L}] done in {(time.time() - t0) / 60:.1f} min")
    print("\nnext: 07b_calibrate.py on procapnet_heldout_pairs_<line>.tsv and procapnet_trainfold_pairs_<line>.tsv")


if __name__ == "__main__":
    main()
