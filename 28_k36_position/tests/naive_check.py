"""Independent, loop-based recount of step-28 strata for the synthetic data; compares with the script's gene table."""
import gzip, importlib.util, numpy as np, pandas as pd, pyBigWig
from pyfaidx import Fasta
spec = importlib.util.spec_from_file_location("s10", "10_shuffle_test.py"); s10 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s10)
CH = ["chr1", "chr2"]
fa = Fasta("synth.fa", as_raw=True, sequence_always_upper=False)
genes, exons, tss = s10.parse_gtf("synth.gtf.gz", set(CH))
peaks = {c: [] for c in CH}
for line in gzip.open("peaks.bed.gz", "rt"):
    f = line.split("\t"); peaks[f[0]].append((int(f[1]), int(f[2])))
def sig(paths, c, L):
    a = np.zeros(L)
    for p in paths:
        bw = pyBigWig.open(p)
        for s, e, v in (bw.intervals(c) or []):
            a[s:e] += abs(v)
    return a
G = pd.read_csv("out28/28_SYN_genes.tsv", sep="\t")
bad = 0; n = 0
for c in CH:
    seq = str(fa[c][:]); L = len(seq)
    ok = np.array([ch.upper() != "N" for ch in seq])
    P = sig(["plus_rep1.bw", "plus_rep2.bw"], c, L); M = sig(["minus_rep1.bw", "minus_rep2.bw"], c, L)
    inpk = np.zeros(L, bool)
    for s, e in peaks[c]:
        inpk[s:e] = True
    glist, ig = s10.build_masks(c, L, ok, genes, exons, tss)
    for r, s, e, sel in glist:
        row = G[(G.chrom == c) & (G.gene_id == r.gene_id)].iloc[0]
        for T, tag in ((10, "r10"), (3, "r3")):
            for strand_name, arr in (("sense", P if r.strand == "+" else M), ("anti", M if r.strand == "+" else P)):
                pos = [p for p in range(s, e) if sel[p - s] and arr[p] >= T]
                clusters, cur = [], []
                for p in pos:
                    if cur and p - cur[-1] > 50:
                        clusters.append(cur); cur = []
                    cur.append(p)
                if cur:
                    clusters.append(cur)
                cnt = {k: 0 for k in ["all", "in", "out", "P1", "P2", "P3", "P4", "P5"]}
                for cl in clusters:
                    best = max(cl, key=lambda p: (arr[p], -p))
                    cnt["all"] += 1
                    cnt["in" if inpk[best] else "out"] += 1
                    rel = (best - s) / (e - s) if r.strand == "+" else (e - 1 - best) / (e - s)
                    cnt[f"P{min(int(rel * 5), 4) + 1}"] += 1
                for k, v in cnt.items():
                    n += 1
                    if row[f"{tag}_{k}_{strand_name}"] != v:
                        bad += 1
                        print("DIFF", c, r.gene_id, tag, k, strand_name, row[f"{tag}_{k}_{strand_name}"], v)
print(f"compared {n} counts, {bad} differ")
