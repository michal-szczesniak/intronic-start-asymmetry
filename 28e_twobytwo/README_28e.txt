PURGE step 28e (10 Nov 2026), EXPLORATORY: near the gene's start (< 5 kb from the 5' end) x near a CpG island
(<= 2 kb), for the pair bias (Results 1) and the deficit of strong sites (Results 4), per line and pooled over both
lines (genes resampled jointly). Specified in the dated project log before running. Computed from result tables only:
  28d_cpg/results: 28d_pairs.tsv, 28d_sites.tsv.gz, 28d_strata.tsv (md5 in 28d_outputs.md5; archive
  28d_results.tar.gz bc7e3f14391ae7ad13d948d6d723bf52)
  28_k36_position/results: 28_K562_genes.tsv (2286575131e7e76cb06812f65fd336ef), 28_HCT116_genes.tsv
  (0eb5f2cc15a26096d6355653760fbaae)
  17_observed/17_observed_sites.py (659dca46e47fe9a9f35df48cef694889)
Run (purge environment, about 1 min):
  python3 28e_twobytwo.py ../28d_cpg/results ../28_k36_position/results ../17_observed/17_observed_sites.py results/28e
  python3 tests/naive_check_28e.py ../28d_cpg/results ../28_k36_position/results results/28e
First run 10 Nov 2026 on copies of these tables (results/28e_cells.tsv, 28e.log, 28e_check.log). Tested before on
the synthetic tables of step 28d (24 point estimates equal to a direct computation).
