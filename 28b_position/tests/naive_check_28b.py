"""Independent loop-based recount for step 28b on the synthetic data (run in the synthetic test folder).
Usage: python3 naive_check_28b.py <gtf> <genes28.tsv> <pairs.tsv> <arms.tsv> <28b prefix> [strong reads, default 10]"""
import gzip, re, sys, importlib.util, numpy as np, pandas as pd, pyBigWig
from pyfaidx import Fasta
gtf, genes28, pairs, arms, pre = sys.argv[1:6]
spec = importlib.util.spec_from_file_location("s10", "10_shuffle_test.py"); s10 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s10)
CH = ["chr1", "chr2"]; T = int(sys.argv[6]) if len(sys.argv) > 6 else 10
BINS = [1, 5, 10, 20, 40, 80]; LAB = ["1-5", "5-10", "10-20", "20-40", "40-80", ">=80"]
def blab(kb):
    k = 0
    for i, b in enumerate(BINS):
        if kb >= b: k = i
    return LAB[k]
fa = Fasta("synth.fa", as_raw=True, sequence_always_upper=False)
genes, exons, tss = s10.parse_gtf(gtf, set(CH))
tx = {}
for line in gzip.open(gtf, "rt"):
    f = line.split("\t")
    if f[2] == "transcript":
        g = re.search('gene_id "([^"]+)"', f[8]).group(1)
        tx.setdefault(g, []).append(int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1)
G28 = pd.read_csv(genes28, sep="\t").set_index("gene_id")
cnt = {}
for c in CH:
    seq = str(fa[c][:]); L = len(seq)
    ok = np.array([ch.upper() != "N" for ch in seq])
    P = np.zeros(L); M = np.zeros(L)
    for s, e, v in (pyBigWig.open("plus_all.bw").intervals(c) or []): P[s:e] += abs(v)
    for s, e, v in (pyBigWig.open("minus_all.bw").intervals(c) or []): M[s:e] += abs(v)
    glist, _ = s10.build_masks(c, L, ok, genes, exons, tss)
    for r, s, e, sel in glist:
        if r.gene_type != "protein_coding" or G28.loc[r.gene_id, "activity"] not in ("A1", "A2", "A3"):
            continue
        for nm, arr in (("sense", P if r.strand == "+" else M), ("anti", M if r.strand == "+" else P)):
            pos = [p for p in range(s, e) if sel[p - s] and arr[p] >= T]
            cl, cur = [], []
            for p in pos:
                if cur and p - cur[-1] > 50: cl.append(cur); cur = []
                cur.append(p)
            if cur: cl.append(cur)
            for x in cl:
                b = max(x, key=lambda p: (arr[p], -p))
                off = b - s if r.strand == "+" else e - 1 - b
                ups = [t for t in tx[r.gene_id] if (t <= b if r.strand == "+" else t >= b)]
                dt = (b - max(ups) if r.strand == "+" else min(ups) - b) if ups else off
                for col, lab in (("fifth", "P%d" % (min(int(off / (e - s) * 5), 4) + 1)), ("bin5", blab(off / 1e3)), ("bintss", blab(dt / 1e3))):
                    cnt[(col, lab, nm)] = cnt.get((col, lab, nm), 0) + 1
R = pd.read_csv(pre + "_strata.tsv", sep="\t")
bad = n = 0
for col in ("fifth", "bin5", "bintss"):
    d = R[(R.analysis == f"sites_strong_{col}") & (~R.stratum.str.contains("rest"))]
    for r in d.itertuples():
        for nm, v in (("sense", r.sense), ("anti", r.anti)):
            n += 1
            if int(v) != cnt.get((col, r.stratum, nm), 0):
                bad += 1; print("DIFF", col, r.stratum, nm, v, cnt.get((col, r.stratum, nm), 0))
Pp = pd.read_csv(pre + "_pairs.tsv", sep="\t")
for col in ("fifth", "dist_bin"):
    d = R[(R.analysis == f"pairs_{col}") & (~R.stratum.str.contains("rest"))]
    for r in d.itertuples():
        m = Pp[Pp[col] == r.stratum]
        exp = 100 * (1 - 2 ** (np.median(m.DI) - np.median(m.ctrl_DI)))
        n += 1
        if not np.isclose(exp, r.estimate):
            bad += 1; print("DIFF pairs", col, r.stratum, r.estimate, exp)
# positions of the pairs, recomputed
G = pd.read_csv(genes28, sep="\t")
G = G[G.gene_type == "protein_coding"]
M0 = pd.read_csv(pairs, sep="\t").merge(pd.read_csv(arms, sep="\t")[["chrom", "mid", "strand_sense"]], on=["chrom", "mid"])
for r in M0.itertuples():
    g = G[(G.chrom == r.chrom) & (G.strand == r.strand_sense) & (G.start <= r.mid) & (r.mid < G.end)].iloc[0]
    rel = (r.mid - g.start) / (g.end - g.start) if g.strand == "+" else (g.end - 1 - r.mid) / (g.end - g.start)
    row = Pp[(Pp.chrom == r.chrom) & (Pp.mid == r.mid)].iloc[0]
    n += 1
    if row.fifth != "P%d" % (min(int(rel * 5), 4) + 1):
        bad += 1; print("DIFF pair fifth", r.chrom, r.mid)
print(f"compared {n} values, {bad} differ")
