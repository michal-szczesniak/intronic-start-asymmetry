PURGE step 28f (10 Nov 2026), EXPLORATORY: does the near-start effect (strand difference < 5 kb from the gene's 5'
end minus >= 5 kb) grow with host activity? Specified in the dated project log before running. Pairs: thirds of
host_signal over the host genes of the pairs (per line); strong sites: activity thirds A1-A3 of step 28. Per line
and pooled, genes resampled jointly for both lines (2,000). Readout: interaction = effect in the most active third
minus effect in the least active third; secondary: the same > 2 kb from any CpG island.
Inputs (result tables only): 28d_cpg/results/28d_pairs.tsv, 28d_sites.tsv.gz, 28d_strata.tsv (archive
28d_results.tar.gz bc7e3f14391ae7ad13d948d6d723bf52); 28_k36_position/results/28_K562_genes.tsv
(2286575131e7e76cb06812f65fd336ef), 28_HCT116_genes.tsv (0eb5f2cc15a26096d6355653760fbaae);
17_observed/17_observed_sites.py (659dca46e47fe9a9f35df48cef694889).
Run (purge environment, about 2 min):
  python3 28f_activity.py ../28d_cpg/results ../28_k36_position/results ../17_observed/17_observed_sites.py results/28f
  python3 tests/naive_check_28f.py ../28d_cpg/results ../28_k36_position/results results/28f
First run 10 Nov 2026 on copies of these tables (results/28f_cells.tsv, 28f.log, 28f_check.log: 72 point estimates
equal to a direct computation). Tested before on the synthetic tables of step 28d (72 values equal).
