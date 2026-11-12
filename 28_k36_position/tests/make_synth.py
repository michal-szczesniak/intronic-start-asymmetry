"""Synthetic genome, annotation, PRO-cap bigWigs (2 replicates per strand) and H3K36me3 peaks for testing step 28."""
import gzip, numpy as np, pyBigWig
rng = np.random.default_rng(7)
CH = {"chr1": 3_000_000, "chr2": 2_000_000}
# genome: random, some lowercase (repeats), a few N
with open("synth.fa", "w") as fh:
    for c, L in CH.items():
        s = rng.choice(list("ACGT"), L)
        low = rng.random(L) < 0.2
        s = np.where(low, np.char.lower(s), s)
        s[1000:1200] = "N"
        fh.write(f">{c}\n")
        seq = "".join(s)
        for i in range(0, L, 60):
            fh.write(seq[i:i + 60] + "\n")
genes = []  # (chrom, start0, end, strand, id, type, exons list)
gi = 0
for c, L in CH.items():
    pos = 20_000
    while pos < L - 60_000:
        if 1_000_000 <= pos < 1_250_000:          # a gene-free region for intergenic blocks
            pos = 1_250_000
            continue
        glen = int(rng.integers(8_000, 40_000))
        st = "+" if rng.random() < 0.5 else "-"
        tp = "protein_coding" if rng.random() < 0.8 else "lncRNA"
        s0, e = pos, pos + glen
        ex = [(s0, s0 + 300)]
        for k in range(int(rng.integers(1, 4))):
            x = int(rng.integers(s0 + 2000, e - 1000))
            ex.append((x, x + 150))
        ex.append((e - 500, e))
        ex = sorted(ex)
        genes.append((c, s0, e, st, f"G{gi:03d}", tp, ex))
        gi += 1
        pos = e + int(rng.integers(2_000, 30_000))
with gzip.open("synth.gtf.gz", "wt") as fh:
    for c, s0, e, st, gid, tp, ex in genes:
        at = f'gene_id "{gid}"; gene_type "{tp}"; gene_name "N{gid}";'
        fh.write(f"{c}\tS\tgene\t{s0 + 1}\t{e}\t.\t{st}\t.\t{at}\n")
        fh.write(f"{c}\tS\ttranscript\t{s0 + 1}\t{e}\t.\t{st}\t.\t{at} transcript_id \"T{gid}\";\n")
        for a, b in ex:
            fh.write(f"{c}\tS\texon\t{a + 1}\t{b}\t.\t{st}\t.\t{at} transcript_id \"T{gid}\";\n")
# signal: dense integer arrays per strand
P = {c: np.zeros(L, np.int64) for c, L in CH.items()}
M = {c: np.zeros(L, np.int64) for c, L in CH.items()}
peaks = []
for k, (c, s0, e, st, gid, tp, ex) in enumerate(genes):
    tssp = s0 if st == "+" else e - 1
    level = [0, 0, 2, 8, 30, 120][k % 6]                        # host activity classes
    S_, A_ = (P[c], M[c]) if st == "+" else (M[c], P[c])
    S_[tssp] += level
    A_[tssp + (-150 if st == "+" else 150)] += level // 3        # some antisense TSS reads
    if level >= 8:                                               # H3K36me3 over the 3' 60% of active genes
        if st == "+":
            peaks.append((c, s0 + int(0.4 * (e - s0)), e))
        else:
            peaks.append((c, s0, e - int(0.4 * (e - s0))))
    # intronic sites: antisense strong sites everywhere; sense strong sites rarer in the 3' part of active genes
    for strand_arr, is_sense in ((S_, True), (A_, False)):
        for _ in range(int(rng.integers(5, 25))):
            x = int(rng.integers(s0, e))
            rel = (x - s0) / (e - s0) if st == "+" else (e - 1 - x) / (e - s0)
            if is_sense and level >= 8 and rel > 0.4 and rng.random() < 0.6:
                continue
            v = int(rng.choice([1, 2, 3, 5, 10, 12, 25]))
            strand_arr[x] += v
            if rng.random() < 0.4:                               # a neighbour within 50 nt (merged site)
                strand_arr[min(x + int(rng.integers(5, 50)), len(strand_arr) - 1)] += int(rng.choice([1, 4, 11]))
# intergenic signal
for c, L in CH.items():
    for _ in range(400):
        x = int(rng.integers(0, L))
        (P if rng.random() < 0.5 else M)[c][x] += int(rng.choice([1, 3, 10, 15]))
    peaks.append((c, L - 30_000, L - 20_000))                    # an intergenic peak
# split into two replicates and write bigWigs (minus strand negative, as some pipelines do)
def write_bw(path, arrs, sign):
    bw = pyBigWig.open(path, "w")
    bw.addHeader([(c, CH[c]) for c in CH])
    for c in CH:
        x = arrs[c]
        nz = np.flatnonzero(x)
        if len(nz):
            bw.addEntries([c] * len(nz), nz.tolist(), ends=(nz + 1).tolist(), values=(sign * x[nz]).astype(float).tolist())
    bw.close()
for nm, D, sign in (("plus", P, 1), ("minus", M, -1)):
    r1 = {c: rng.binomial(D[c], 0.5) for c in CH}
    r2 = {c: D[c] - r1[c] for c in CH}
    write_bw(f"{nm}_rep1.bw", r1, sign)
    write_bw(f"{nm}_rep2.bw", r2, sign)
    write_bw(f"{nm}_all.bw", D, sign)
with gzip.open("peaks.bed.gz", "wt") as fh:
    for c, s, e in sorted(peaks):
        fh.write(f"{c}\t{s}\t{e}\tpk\t1000\t.\n")
print(len(genes), "genes;", sum(1 for g in genes if g[5] == "protein_coding"), "protein-coding;", len(peaks), "peaks")
