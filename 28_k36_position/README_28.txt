PURGE step 28 (10 Nov 2026): H3K36me3, position along the gene and HCT116
for the gene-strand deficit of strong measured start sites (follow-up of step 17)

Questions and primary analyses (specified before any data were seen; also recorded in the dated project log)
  A1  deficit of strong sites inside minus outside H3K36me3 peaks, same active genes; prediction > 0
  A2  deficit in the 3' fifth minus the 5' fifth of the gene; prediction > 0
  B   deficit of strong sites in active genes of HCT116; prediction > 0
  Primary: active protein-coding genes (host-activity thirds pooled, genes resampled within thirds), strong
  sites; K562: A1, A2; HCT116: A1, A2, B. Two-sided bootstrap (2,000 resamples of genes), 95% CIs, no
  correction for multiple testing. Strong = max(3, round(10 x N / 18,699,448)) reads, N = autosomal reads
  (10 reads in K562, as in step 17; the same reads per million in HCT116). Everything else is secondary.
  Controls: silent genes; intergenic 1-Mb blocks (plus vs minus strand).
  HCT116 bigWigs: the set of ENCODE's default analysis of ENCSR165QKR is primary (fallback, if neither or
  both sets are in it: set 266IRD, processed like the K562 data); the other set is a sensitivity analysis.

Files
  28_k36_position.py      the analysis (loads step-10 and step-17 code; reproduces step 17's K562 counts exactly
                          or stops)
  28_download.py          ENCODE metadata check, downloads with md5, choice of the primary HCT116 set
  run_28.sh               runs everything: md5 checks, synthetic test, downloads, three runs, checksums, archive
  tests/make_synth.py     synthetic genome, annotation, bigWigs (two replicates per strand) and peaks
  tests/naive_check.py    independent loop-based recount of every stratum, compared with the script's table
  MD5SUMS                 checksums of the files above

Run (inside srun with --mem=32G; about 1 h)
  bash run_28.sh
Outputs (project_data/.../PURGE/28_k36_position/results)
  28_synth_test.log, 28_download.log, 28_downloads.tsv, 28_hct116_sets.tsv,
  28_<line>.log, 28_<line>_genes.tsv, 28_<line>_blocks.tsv, 28_<line>_strata.tsv for K562, HCT116 and
  HCT116_sensitivity, 28_outputs.md5; archive scratch/.../PURGE/28_results.tar.gz (code and results)

Tested on 10 Nov 2026 on synthetic data: step 17 reproduced exactly, 3,424 stratum counts equal to an
independent recount, strands swapped automatically when the bigWigs are exchanged, stops on a changed step-17
table, on unknown arguments and on peaks with other chromosome names; run_28.sh run end to end on a local copy
of the directory tree with synthetic inputs and a mocked ENCODE.
Version b (10 Nov 2026, 15:10): 28_download.py retries ENCODE requests after server errors (the first run
stopped on a 504 Gateway Time-out of the portal); tested on a mocked ENCODE that fails twice, fails a download
once, answers 404 (stops at once) or stays down (stops after 6 attempts). The analysis script is unchanged.
