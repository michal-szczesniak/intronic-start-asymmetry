PURGE step 28d (10 Nov 2026), EXPLORATORY: CpG islands and the strand difference within 5 kb of the gene's start.
Question (specified before running): after removing strong sites and pairs within 2 kb of any CpG island (UCSC
hg38 cpgIslandExt), is the deficit of strong sites, and the bias of pairs, within 1-5 kb of the gene's 5' end still
above zero (95% CI)? If yes, CpG islands do not explain the concentration near the gene's start; if the estimate
falls to about the level beyond 5 kb, the concentration is tied to CpG islands. Margins of 1 and 5 kb, the subset
within 2 kb of an island and the difference < 5 kb minus >= 5 kb are also reported; no correction; exploratory.
Files
  28d_cpg.py              the analysis (step-28b sites and pairs, step-28 gene tables, step-17 bootstrap helpers)
  28d_download.py         downloads cpgIslandExt.txt.gz from UCSC (retries after server errors), records md5
  run_28d.sh              md5 checks, synthetic test, download, run, checksums, archive
  tests/make_cgi.py       how the synthetic island table was made
  tests/naive_check_28d.py  independent loop-based check (distances, counts and estimates for the 2-kb margin)
  tests/fixtures/         synthetic step-28 and 28b tables (labelled K562 and HCT116, but synthetic) and islands
Run: sbatch --mem=16G (a few minutes). Results in project_data/.../PURGE/28d_cpg/results; archive
scratch/.../PURGE/28d_results.tar.gz.
Tested on 10 Nov 2026: island distances equal brute force (300 random cases); synthetic tables: counts and
estimates equal an independent recount; the real step-28b tables with a random island file: checks pass; stops
below 20,000 islands; run_28d.sh run end to end on a local copy of the directory tree with a mocked download.
