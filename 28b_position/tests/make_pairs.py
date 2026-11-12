"""Synthetic matched pairs and arm table for testing step 28b: pairs inside protein-coding genes of the synthetic
genome, DI lower (weaker sense start) near the 5' end. Hosts are given by gene_id, except every seventh pair, whose
host gets an identifier absent from the annotation (to test placement by position); with 'name', every host is
given by gene name. Usage: python3 make_pairs.py genes28.tsv pairs.tsv arms.tsv [name]"""
import sys, numpy as np, pandas as pd
rng = np.random.default_rng(3)
G = pd.read_csv(sys.argv[1], sep="\t")
G = G[(G.gene_type == "protein_coding") & (G.end - G.start > 5000)]
use_name = len(sys.argv) > 4 and sys.argv[4] == "name"
rows, arms = [], []
k = 0
for r in G.itertuples(index=False):
    for _ in range(int(rng.integers(1, 6))):
        mid = int(rng.integers(r.start + 1500, r.end - 500))
        rel = (mid - r.start) / (r.end - r.start) if r.strand == "+" else (r.end - 1 - mid) / (r.end - r.start)
        di = rng.normal(-1.0 if rel < 0.2 else -0.2, 1.0)
        k += 1
        host = r.gene_name if use_name else (f"ENSGRETIRED{k:05d}" if k % 7 == 0 else r.gene_id)
        rows.append({"line": "SYN", "chrom": r.chrom, "mid": mid, "host": host, "DI": round(di, 4),
                     "ctrl_chrom": "chr2", "ctrl_mid": 1_100_000, "ctrl_DI": round(rng.normal(0, 1.0), 4)})
        arms.append({"chrom": r.chrom, "mid": mid, "host": host, "strand_sense": r.strand,
                     "strand_antisense": "-" if r.strand == "+" else "+"})
M = pd.DataFrame(rows)
M.to_csv(sys.argv[2], sep="\t", index=False)
pd.DataFrame(arms).to_csv(sys.argv[3], sep="\t", index=False)
print(len(M), "pairs; expected overall bias", round(100 * (1 - 2 ** (np.median(M.DI) - np.median(M.ctrl_DI))), 4))
