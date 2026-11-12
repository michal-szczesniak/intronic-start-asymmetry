PURGE step 28b (10 Nov 2026), EXPLORATORY: position along the host gene of
(Q1) the weaker sense start of intronic start-site pairs (Results 1, step 23b pairs) and
(Q2) the gene-strand deficit of strong measured start sites (Results 4, step 28).
Specified after the step-28 results, which placed the deficit of strong sites in the 5' fifth of active genes
(opposite to the H3K36me3 prediction). No primary contrast; two-sided bootstrap P without correction.

Q1  bias 100 x (1 - 2^d), d = median DI of intronic pairs - median DI of their intergenic partners, by fifth of the
    host gene (5' -> 3') and by distance of the pair from the gene's 5' end (kb: 1-5, 5-10, 10-20, 20-40, 40-80,
    >= 80); contrast: 5' fifth (or 1-5 kb) minus the rest; CIs: 2,000 resamples of host genes. Check: over all
    pairs the estimate must equal step 23b (K562 33.0%, HCT116 20.6%), or the script stops. Host gene: by its
    gene_id in GENCODE v47 (any type) or, failing that, the only protein-coding gene containing the pair's
    midpoint on the strand of the sense start (iSTART arm table); stops if fewer than 98% of pairs are placed.
Q2  deficit of strong sites (K562 >= 10, HCT116 >= 24 reads, as in step 28) in active protein-coding genes by
    fifth, by distance from the gene's 5' end and from the nearest annotated TSS of the gene upstream of the
    site, and by fifth within thirds of gene length; >= 10 reads as secondary. CIs: 2,000 resamples of genes
    within activity thirds. Check: per-gene counts must equal step 28's tables, or the script stops.

Files
  28b_position.py              the analysis (loads the code of steps 10, 17 and 28)
  run_28b.sh                   md5 checks, synthetic test, the K562 and HCT116 runs, checksums, archive
  tests/make_synth_alt.py      adds alternative transcripts to step 28's synthetic annotation
  tests/make_pairs.py          synthetic matched pairs (hosts given by gene name)
  tests/naive_check_28b.py     independent loop-based recount of every stratum, compared with the script
Run: sbatch --mem=32G (about 15 min). Outputs in project_data/.../PURGE/28b_position/results; archive
scratch/.../PURGE/28b_results.tar.gz.
Tested on 10 Nov 2026: synthetic data (310 values equal to an independent recount; host genes found by gene_id
and by gene name; distances from alternative TSSs), stops on a changed step-28 table, a wrong step-23b value and
unplaceable pairs; run_28b.sh run end to end on a local copy of the directory tree with synthetic inputs.
Version b (10 Nov 2026, 17:10): version a stopped on the cluster because 127 of 1,610 K562 pairs had a host gene_id
not in the step-28 table (genes without eligible positions or identifiers outside GENCODE v47); hosts are now looked
up in the whole annotation, with placement by position and strand as fallback (tested with absent identifiers and
with hosts given by name).
Version c (10 Nov 2026, 17:15): version b stopped because 42 K562 pairs lie on chrX, whose genes had not been
read (the single-site analysis covers the 22 autosomes only); gene coordinates are now read for every chromosome
of the pairs, the number of such pairs is printed, and the bias is also given for autosomal pairs only.
