#!/usr/bin/env python3
"""Step 26b: gene-body DNA methylation and the strand deficit of measured start sites (K562).

Question. Step 17 found fewer strong PRO-cap start sites (>= 10 reads) on the gene strand of introns,
clearly so only in the more active genes. Inside transcribed genes, DNA methylation laid down after
H3K36me3 suppresses spurious initiation (Neri et al. 2017). Is the gene-strand deficit of strong sites
larger in genes whose bodies are more methylated, at the same level of host activity?

Inputs
  WGBS     ENCODE CpG methylation files of K562 (ENCSR765JPC): either bedMethyl (11 fields) or ENCODE
           bed9+ (14 fields; field 12 = reference context, rows kept when it starts with CG). Fields 10
           and 11 are read as coverage and percent methylated; the format is checked, not assumed.
           Several files (replicates) are pooled per CpG: methylated reads = round(coverage * percent /
           100), summed over strands and files. A CpG is used when its pooled coverage >= --min-cov.
  step 17  per-gene and per-block tables (17_obs_genes.tsv, 17_obs_blocks.tsv) and, as in step 17,
           the Puffin-D (13b) and ProCapNet (16d) gene tables.
  step 10  masks (eligible intronic positions per gene; intergenic positions) rebuilt with the same
           code, genome and annotation as step 17; the eligible nt per gene must match step 17.

Methylation of a gene = mean methylation of the CpGs (C at an eligible position) with coverage >=
--min-cov; genes with fewer than --min-cpg such CpGs are left out (counts printed). Intergenic 1-Mb
blocks: the same over intergenic positions.

Analyses (deficit = 100 * (1 - sense / antisense), summed over genes; CIs: 2,000 bootstrap resamples
of genes, stratified by the cells that define the groups; two-sided bootstrap p)
  1. QC: format; coordinate check against the genome (the C of each CpG must sit at a C followed by G
     on the plus strand), pooled and for a sample of rows of every file and strand (STOP below 95%);
     coverage; agreement between files; methylation by host-activity group.
  2. PRIMARY: active protein-coding genes (thirds A1-A3 of host activity, as in step 17), divided into
     thirds of gene-body methylation WITHIN each activity third (M1 low, M3 high), pooled over A1-A3;
     metric: measured sites >= 10 reads; contrast M3 - M1.
  3. The same for every activity group (silent, low, A1, A2, A3) and for the other metrics: measured
     sites >= 1 and >= 3 reads and summed reads; Puffin-D and ProCapNet weak, strong and summed signal
     (ProCapNet sees only local sequence: a sequence-only control).
  4. Negative control: intergenic 1-Mb blocks in thirds of methylation, plus strand vs minus strand
     (expected ratio 1 in every third).

Usage (CPU; files are read in chunks; peak memory estimated at about 10 GB for two ENCODE CpG files):
  python3 26b_methylation.py --meth ENCFF660IHA.bed.gz,ENCFF328NMN.bed.gz \
      --genes17 ../17_observed/results/17_obs_genes.tsv --blocks17 ../17_observed/results/17_obs_blocks.tsv \
      --chroms chr1,...,chr22 --genome hg38.fa --gtf gencode.v47.annotation.gtf.gz \
      --step10 ../10_shuffle/10_shuffle_test.py --step17 ../17_observed/17_observed_sites.py \
      --puffin-dir ../13_genomewide/results/shuffle --procapnet-dir ../16_procapnet/results/shuffle_pc \
      --out-prefix results/26b_meth
"""
import argparse
import gzip
import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ACT_ORDER = ["silent", "low", "A1", "A2", "A3"]
M_ORDER = ["M1", "M2", "M3"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read_meth(path, chroms, rng, chunksize=4_000_000):
    """One CpG methylation file, read in chunks -> {chrom: (cpg, cov, meth, plus)}, cpg = 0-based C of the CpG
    on the plus strand (plus-strand rows: start; minus-strand rows: start - 1), plus = row on the plus strand.
    Also returns a sample of up to 20,000 rows per chromosome and strand for the coordinate check."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline().rstrip("\n").split("\t")
    nf = len(first)
    if nf not in (11, 14):
        sys.exit(f"STOP: {path}: {nf} fields in line 1; expected 11 (bedMethyl) or 14 (ENCODE bed9+)")
    cols = [0, 1, 2, 5, 9, 10] + ([11] if nf == 14 else [])
    dt = {0: str, 1: np.int64, 2: np.int64, 5: str, 9: np.int64, 10: np.float64}
    if nf == 14:
        dt[11] = str
    names = ["chrom", "start", "end", "strand", "cov", "pct"] + (["ctx"] if nf == 14 else [])
    n0 = n_w1 = n_plus = n_minus = bad_strand = bad_pct = bad_cov = 0
    ctx_counts = {}
    acc = {c: [] for c in chroms}
    for x in pd.read_csv(path, sep="\t", header=None, usecols=cols, dtype=dt, chunksize=chunksize):
        x.columns = names
        n0 += len(x)
        if nf == 14:
            for k, v in x["ctx"].value_counts().items():
                ctx_counts[k] = ctx_counts.get(k, 0) + int(v)
            x = x[x["ctx"].str.startswith("CG", na=False)]
        w = (x["end"] - x["start"]).to_numpy()
        st = x["strand"].to_numpy()
        n_w1 += int((w == 1).sum())
        ip, im = st == "+", st == "-"
        n_plus += int(ip.sum())
        n_minus += int(im.sum())
        bad_strand += int((~(ip | im)).sum())
        bad_pct += int(((x["pct"] < 0) | (x["pct"] > 100)).sum())
        bad_cov += int((x["cov"] < 0).sum())
        x = x[x["chrom"].isin(chroms)]
        for c, d in x.groupby("chrom"):
            pl = (d["strand"] == "+").to_numpy()
            stt = d["start"].to_numpy(np.int64)
            cv = d["cov"].to_numpy(np.int64)
            acc[c].append((np.where(pl, stt, stt - 1), cv, np.rint(cv * d["pct"].to_numpy() / 100.0).astype(np.int64), pl))
    n_cg = n_plus + n_minus + bad_strand
    print(f"  {os.path.basename(path)}: {nf} fields, {n0:,} rows")
    if nf == 14:
        top = sorted(ctx_counts.items(), key=lambda kv: -kv[1])[:5]
        print("    reference context (field 12), top: " + ", ".join(f"{k} {v:,}" for k, v in top))
    print(f"    CpG rows {n_cg:,}; width 1: {100 * n_w1 / max(n_cg, 1):.2f}%; strands + {n_plus:,}, - {n_minus:,}, other {bad_strand:,}; "
          f"percent outside 0-100: {bad_pct:,}; negative coverage: {bad_cov:,}")
    if n_cg == 0 or (n_cg - n_w1) > 0.001 * n_cg or bad_strand or bad_pct or bad_cov:
        sys.exit(f"STOP: {path}: unexpected values (see the line above)")
    out, samp = {}, {}
    for c in chroms:
        if not acc[c]:
            continue
        cpg, cv, mt, pl = (np.concatenate(z) for z in zip(*acc[c]))
        acc[c] = None
        out[c] = (cpg, cv, mt, pl)
        for sd, m_ in (("+", pl), ("-", ~pl)):
            v = cpg[m_]
            samp[(c, sd)] = v if len(v) <= 20000 else rng.choice(v, 20000, replace=False)
    nr = sum(len(v[1]) for v in out.values())
    cov_all = np.concatenate([v[1] for v in out.values()]) if out else np.zeros(1)
    print(f"    on the analysed chromosomes: {nr:,} rows; coverage median {np.median(cov_all):.0f}, "
          f"rows with coverage >= 10: {100 * float((cov_all >= 10).mean()):.1f}%")
    return out, samp


def per_cpg(cpg, cov, meth):
    """Sum coverage and methylated reads over the two strands of each CpG -> sorted unique cpg, cov, meth."""
    u, inv = np.unique(cpg, return_inverse=True)
    return u, np.bincount(inv, weights=cov, minlength=len(u)), np.bincount(inv, weights=meth, minlength=len(u))


def stats_in(pos, beta, lo, hi, msk, off):
    """Mean methylation of the CpGs at positions lo <= pos < hi with msk[pos - off] true."""
    i0, i1 = np.searchsorted(pos, lo), np.searchsorted(pos, hi)
    if i1 <= i0:
        return 0, np.nan
    p = pos[i0:i1]
    keep = msk[p - off]
    b = beta[i0:i1][keep]
    return int(keep.sum()), (float(b.mean()) if len(b) else np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meth", required=True, help="comma-separated CpG methylation files, pooled")
    ap.add_argument("--genes17", required=True)
    ap.add_argument("--blocks17", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--step10", required=True)
    ap.add_argument("--step17", required=True)
    ap.add_argument("--puffin-dir", default="")
    ap.add_argument("--procapnet-dir", default="")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--min-cpg", type=int, default=20)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.time()
    chroms = [c for c in a.chroms.split(",") if c]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    s10 = load_module("s10", os.path.abspath(a.step10))
    s17 = load_module("s17", os.path.abspath(a.step17))
    rng = np.random.default_rng(a.seed)
    print(f"== 26b gene-body methylation and the strand deficit; min coverage {a.min_cov}, min CpGs {a.min_cpg}, "
          f"bootstrap {a.n_boot}, seed {a.seed}")

    # ---- 1. methylation files
    print("\n== 1. methylation files ==")
    files = [f for f in a.meth.split(",") if f]
    parts, samp = [], []
    for f in files:
        o, sm = read_meth(f, set(chroms), rng)
        parts.append(o)
        samp.append(sm)
    fchk = {(i, sd): [0, 0] for i in range(len(parts)) for sd in "+-"}
    agr = np.zeros(6)  # n, sx, sy, sxx, syy, sxy of per-CpG methylation, files 1 and 2, coverage >= min-cov in both
    byc, n_pool, s_pool = {}, 0, 0.0
    for c in chroms:
        have = [p[c] for p in parts if c in p]
        if not have:
            continue
        if len(parts) >= 2 and c in parts[0] and c in parts[1]:
            q1, q2 = per_cpg(*parts[0][c][:3]), per_cpg(*parts[1][c][:3])
            _, i1_, i2_ = np.intersect1d(q1[0], q2[0], assume_unique=True, return_indices=True)
            k_ = (q1[1][i1_] >= a.min_cov) & (q2[1][i2_] >= a.min_cov)
            bx, by = q1[2][i1_][k_] / q1[1][i1_][k_], q2[2][i2_][k_] / q2[1][i2_][k_]
            agr += [len(bx), bx.sum(), by.sum(), (bx * bx).sum(), (by * by).sum(), (bx * by).sum()]
            del q1, q2
        u, cv, mt = per_cpg(np.concatenate([h[0] for h in have]), np.concatenate([h[1] for h in have]),
                            np.concatenate([h[2] for h in have]))
        k_ = cv >= a.min_cov
        byc[c] = (u[k_], mt[k_] / cv[k_])
        n_pool += int(k_.sum())
        s_pool += float(byc[c][1].sum())
        for p in parts:
            p.pop(c, None)
    del parts
    if len(files) >= 2 and agr[0] > 2:
        n_, sx, sy, sxx, syy, sxy = agr
        r = (sxy - sx * sy / n_) / np.sqrt((sxx - sx * sx / n_) * (syy - sy * sy / n_))
        print(f"  agreement of files 1 and 2: {int(n_):,} CpGs with coverage >= {a.min_cov} in both; Pearson r of methylation {r:.3f}; "
              f"mean {sx / n_:.3f} and {sy / n_:.3f}")
    print(f"  pooled: {n_pool:,} CpGs with coverage >= {a.min_cov}; mean methylation {s_pool / max(n_pool, 1):.3f}")

    # ---- per chromosome: masks, coordinate check, gene and block methylation
    from pyfaidx import Fasta
    fa = Fasta(os.path.abspath(a.genome), as_raw=True, sequence_always_upper=False)
    genes, exons, tss = s10.parse_gtf(os.path.abspath(a.gtf), set(chroms))
    grows, brows = [], []
    chk_n = chk_ok = 0
    for c in chroms:
        t1 = time.time()
        raw = np.frombuffer(str(fa[c][:]).encode(), dtype=np.uint8)
        up = np.where(raw >= 97, raw - 32, raw).astype(np.uint8)
        del raw
        L = len(up)
        ok = up != ord("N")
        for i, dd in enumerate(samp):
            for sd in "+-":
                v = dd.get((c, sd))
                if v is None:
                    continue
                v = v[(v >= 0) & (v + 1 < L)]
                fchk[(i, sd)][0] += len(v)
                fchk[(i, sd)][1] += int(((up[v] == ord("C")) & (up[v + 1] == ord("G"))).sum())
        d = byc.pop(c, None)
        if d is None or not len(d[0]):
            print(f"  {c}: no CpGs")
            continue
        pos, beta = d
        inside = (pos >= 0) & (pos + 1 < L)
        pos, beta = pos[inside], beta[inside]
        isCG = (up[pos] == ord("C")) & (up[pos + 1] == ord("G"))
        chk_n += len(pos)
        chk_ok += int(isCG.sum())
        glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
        for r, s, e, sel in glist:
            n, mm = stats_in(pos, beta, s, e, sel, s)
            grows.append({"chrom": c, "gene_id": r.gene_id, "eligible_nt_26b": int(sel.sum()), "meth_cpg": n, "meth": mm})
        for b0 in range(0, L, s10.BLOCK):
            b1 = min(L, b0 + s10.BLOCK)
            if ig[b0:b1].any():
                n, mm = stats_in(pos, beta, b0, b1, ig[b0:b1], b0)
                brows.append({"chrom": c, "block_start": b0, "meth_cpg": n, "meth": mm})
        print(f"  {c}: {len(pos):,} CpGs, {100 * float(isCG.mean()):.2f}% at a CG of the genome; {len(glist):,} genes  "
              f"[{time.time() - t1:.0f} s]", flush=True)
        del up, ok
    frac = chk_ok / chk_n if chk_n else 0
    print(f"  coordinate check, pooled: {100 * frac:.2f}% of {chk_n:,} CpGs sit on a CG of the genome (C at the given position)")
    bad_f = []
    for (i, sd), (nn, kk) in sorted(fchk.items()):
        if nn:
            print(f"  coordinate check, file {i + 1} ({os.path.basename(files[i])}), strand {sd}: {100 * kk / nn:.2f}% of {nn:,} sampled rows")
            if kk / nn < 0.95:
                bad_f.append(f"file {i + 1} strand {sd}")
    if frac < 0.95 or bad_f:
        sys.exit("STOP: the CpG coordinates do not match the genome (convention or assembly)"
                 + (": " + ", ".join(bad_f) if bad_f else ""))
    GM, BM = pd.DataFrame(grows), pd.DataFrame(brows)

    # ---- join with step 17
    G = pd.read_csv(a.genes17, sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str, "gene_name": str})
    B = pd.read_csv(a.blocks17, sep="\t", keep_default_na=False, na_values=[""])
    G = G.merge(GM, on=["chrom", "gene_id"], how="inner")
    mism = int((G["eligible_nt"] != G["eligible_nt_26b"]).sum())
    print(f"\njoined with step 17: {len(G):,} genes; eligible nt differ for {mism:,}")
    if mism > 0.01 * len(G):
        sys.exit("STOP: the masks differ from step 17 (genome, annotation or step-10 code)")
    B = B.merge(BM, on=["chrom", "block_start"], how="inner")
    P_all = G[(G["gene_type"] == "protein_coding") & G["GC"].notna()].reset_index(drop=True)
    J = P_all.copy()
    MOD = []
    for pre, name, dd in (("pf", "Puffin-D", a.puffin_dir), ("pc", "ProCapNet", a.procapnet_dir)):
        if not dd:
            continue
        X = s17.load_model(dd, chroms, pre)
        n0 = len(J)
        J = J.merge(X, on=["chrom", "gene_id"], how="inner")
        print(f"  {name}: joined {len(J):,} of {n0:,} protein-coding genes")
        MOD += [(f"{name} {lb}", f"{pre}_{k}") for k, lb in s17.MODEL_KEYS]
    J = J.reset_index(drop=True)
    jh = J["host_tss_sense"].to_numpy(float)
    act = jh >= s17.ACTIVE_MIN
    lab = np.where(jh == 0, "silent", np.where(jh < s17.ACTIVE_MIN, "low", "active"))
    e = np.quantile(jh[act], [1 / 3, 2 / 3])
    lab = np.where(act, np.array(["A1", "A2", "A3"])[np.searchsorted(e, jh, side="right").clip(0, 2)], lab)
    J["activity"] = lab
    OBS = [("obs r1", "obs_all_r1"), ("obs r3", "obs_all_r3"), ("obs r10", "obs_all_r10"), ("obs sum", "obs_all_sum")]
    MET = OBS + MOD

    print("\n== 1b. methylation by host activity (protein-coding genes, joined) ==")
    for g_ in ACT_ORDER:
        x = J[J["activity"] == g_]
        ok_ = x["meth_cpg"] >= a.min_cpg
        print(f"  {g_:7s} genes {len(x):6,}; with >= {a.min_cpg} CpGs {int(ok_.sum()):6,}; median CpGs {x['meth_cpg'].median():7.0f}; "
              f"methylation median {x.loc[ok_, 'meth'].median():.3f} (IQR {x.loc[ok_, 'meth'].quantile(0.25):.3f}-{x.loc[ok_, 'meth'].quantile(0.75):.3f})")
    U = J[J["meth_cpg"] >= a.min_cpg].reset_index(drop=True)
    ua = U["activity"].isin(["A1", "A2", "A3"]).to_numpy()
    rho = spearmanr(U.loc[ua, "meth"], U.loc[ua, "host_tss_sense"])
    rho_gc = spearmanr(U["meth"], U["GC"])
    print(f"  Spearman, active genes: methylation vs host activity rho {rho.statistic:+.3f} (p {rho.pvalue:.2g}, n {int(ua.sum()):,}); "
          f"all genes: methylation vs GC rho {rho_gc.statistic:+.3f}")

    # methylation thirds within each activity group
    ml = np.full(len(U), "-", dtype=object)
    for g_ in ACT_ORDER:
        m_ = (U["activity"] == g_).to_numpy()
        if m_.sum() >= 6:
            ed = np.quantile(U.loc[m_, "meth"], [1 / 3, 2 / 3])
            ml[m_] = np.array(M_ORDER)[np.searchsorted(ed, U.loc[m_, "meth"].to_numpy(), side="right").clip(0, 2)]
    U["meth_third"] = ml

    rows = []

    def run(df, cell, mlab, metrics, analysis, subset, contrast=("M3", "M1")):
        """Deficit per methylation third; resampling stratified by `cell` (genes resampled within each cell)."""
        res = {}
        for b in M_ORDER:
            d = df[mlab == b]
            if len(d) < 2:
                continue
            cells = [np.flatnonzero((cell[mlab == b] == cc)) for cc in pd.unique(cell[mlab == b])]
            Ws = [s17.boot_W(len(ix), a.n_boot, rng) for ix in cells]
            for lbl, col in metrics:
                S, A = d[f"{col}_sense"].to_numpy(float), d[f"{col}_anti"].to_numpy(float)
                bs = sum(W @ S[ix] for W, ix in zip(Ws, cells))
                ba = sum(W @ A[ix] for W, ix in zip(Ws, cells))
                dr, pt = s17.deficit(bs, ba), float(s17.deficit(S.sum(), A.sum()))
                lo, hi = s17.ci(dr)
                res[(b, lbl)] = (pt, dr)
                rows.append({"analysis": analysis, "subset": subset, "bin": b, "metric": lbl, "n_units": len(d),
                             "sense": S.sum(), "anti": A.sum(), "deficit": pt, "lo": lo, "hi": hi,
                             "median_meth": float(d["meth"].median()), "p_boot": np.nan})
        bh, bl = contrast
        for lbl, _ in metrics:
            if (bh, lbl) in res and (bl, lbl) in res:
                v = res[(bh, lbl)][1] - res[(bl, lbl)][1]
                lo, hi = s17.ci(v)
                rows.append({"analysis": analysis, "subset": subset, "bin": f"{bh}-{bl}", "metric": lbl, "n_units": np.nan,
                             "sense": np.nan, "anti": np.nan, "deficit": res[(bh, lbl)][0] - res[(bl, lbl)][0],
                             "lo": lo, "hi": hi, "median_meth": np.nan, "p_boot": min(1.0, s17.pval(v))})

    # 2. primary: active genes, thirds within A1-A3, pooled
    A_ = U[U["activity"].isin(["A1", "A2", "A3"])].reset_index(drop=True)
    run(A_, A_["activity"].to_numpy(), A_["meth_third"].to_numpy(), MET, "meth_within_activity", "pc_active")
    # 3. every activity group on its own
    for g_ in ACT_ORDER:
        x = U[U["activity"] == g_].reset_index(drop=True)
        run(x, np.full(len(x), g_), x["meth_third"].to_numpy(), MET, "meth_within_activity", g_)
    # 4. intergenic blocks (plus vs minus), thirds of block methylation
    BB = B[B["meth_cpg"] >= a.min_cpg].reset_index(drop=True)
    if len(BB) >= 6:
        ed = np.quantile(BB["meth"], [1 / 3, 2 / 3])
        bl = np.array(M_ORDER)[np.searchsorted(ed, BB["meth"].to_numpy(), side="right").clip(0, 2)]
        run(BB, np.full(len(BB), "blocks"), bl, OBS, "meth_blocks", "intergenic")
        BB["meth_third"] = bl
    R = pd.DataFrame(rows)

    def show(analysis, subset, metrics, title):
        x = R[(R["analysis"] == analysis) & (R["subset"] == subset)]
        if not len(x):
            return
        print(f"\n-- {title} --")
        print("  " + f"{'bin':8s}{'n':>7s}{'meth':>8s}" + "".join(f"{m:>27s}" for m in metrics))
        for b in M_ORDER + ["M3-M1"]:
            xb = x[x["bin"] == b].set_index("metric")
            if not len(xb):
                continue
            n = xb["n_units"].iloc[0]
            md = xb["median_meth"].iloc[0]
            cells = [f"{xb.loc[m, 'deficit']:+.1f} ({xb.loc[m, 'lo']:+.1f}, {xb.loc[m, 'hi']:+.1f})" if m in xb.index else "" for m in metrics]
            print("  " + f"{b:8s}" + (f"{int(n):7,d}" if np.isfinite(n) else " " * 7) + (f"{md:8.3f}" if np.isfinite(md) else " " * 8)
                  + "".join(f"{cc:>27s}" for cc in cells))
            if b == "M3-M1":
                print("  " + " " * 23 + "".join(f"{'p=' + format(xb.loc[m, 'p_boot'], '.2g'):>27s}" if m in xb.index else " " * 27 for m in metrics))

    ol, mlb = [m for m, _ in OBS], [m for m, _ in MOD]
    print("\n== 2. PRIMARY: active protein-coding genes, methylation thirds within each activity third, pooled; "
          "sense deficit (%) with 95% CI ==")
    show("meth_within_activity", "pc_active", ol, "measured")
    show("meth_within_activity", "pc_active", mlb, "models")
    print("\n== 3. each activity group ==")
    for g_ in ACT_ORDER:
        show("meth_within_activity", g_, ol, f"{g_}: measured")
        show("meth_within_activity", g_, mlb, f"{g_}: models")
    print("\n== 4. negative control: intergenic 1-Mb blocks, plus vs minus strand, by methylation third ==")
    show("meth_blocks", "intergenic", ol, "intergenic blocks: measured")

    U.to_csv(a.out_prefix + "_genes.tsv", sep="\t", index=False)
    (BB if len(BB) else B).to_csv(a.out_prefix + "_blocks.tsv", sep="\t", index=False)
    R.to_csv(a.out_prefix + "_strata.tsv", sep="\t", index=False)
    print(f"\nwritten: {a.out_prefix}_genes.tsv, _blocks.tsv, _strata.tsv  [{(time.time() - t0) / 60:.1f} min]")


if __name__ == "__main__":
    main()
