"""Synthetic CpG-island table (UCSC cpgIslandExt format) for the synthetic genome of step 28: islands next to the
5' end of some genes and elsewhere. Usage: python3 make_cgi.py synth.gtf.gz cgi.txt.gz"""
import gzip, sys, numpy as np
rng = np.random.default_rng(11)
rows = []
for line in gzip.open(sys.argv[1], "rt"):
    f = line.split("\t")
    if f[2] != "gene":
        continue
    s, e = int(f[3]) - 1, int(f[4])
    if rng.random() < 0.6:                                # an island over the promoter, reaching into the gene
        a = s - 300 if f[6] == "+" else e - 2500
        rows.append((f[0], max(0, a), max(0, a) + 2800))
for c, L in (("chr1", 3_000_000), ("chr2", 2_000_000)):
    for x in rng.integers(0, L - 2000, 15):
        rows.append((c, int(x), int(x) + int(rng.integers(300, 1500))))
with gzip.open(sys.argv[2], "wt") as fh:
    for c, a, b in sorted(rows):
        fh.write(f"1\t{c}\t{a}\t{b}\tCpG: 40\t{b - a}\t40\t200\t20\t60\t0.8\n")
print(len(rows), "islands")
