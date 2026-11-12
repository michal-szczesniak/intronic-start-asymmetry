PURGE step 28c (10 Nov 2026), EXPLORATORY: bias of pairs and deficit of strong sites beyond the first 5 kb of
the host gene, with 95% CIs (step 28b printed these estimates without CIs). Computed from the result tables only:
  28b_position/results: 28b_K562_pairs.tsv (8019e1b26beec25a1657264b137d3997), 28b_HCT116_pairs.tsv
  (951d27c8f32c37c62314fce81c04feec), 28b_K562_sites.tsv.gz (ac95f13d0aff5be4b1939e12621575c0),
  28b_HCT116_sites.tsv.gz (5fd83afc5e4c854fb3df008f4d7ec4db), 28b_*_strata.tsv (for the check)
  28_k36_position/results: 28_K562_genes.tsv (2286575131e7e76cb06812f65fd336ef), 28_HCT116_genes.tsv
  (0eb5f2cc15a26096d6355653760fbaae)
  17_observed/17_observed_sites.py (659dca46e47fe9a9f35df48cef694889; bootstrap helpers)
Run (in the purge environment, a few minutes):
  python3 28c_beyond5kb.py ../17_observed/17_observed_sites.py ../28b_position/results ../28_k36_position/results 28c_beyond5kb.tsv
The < 5 kb estimates are checked against step 28b (the script stops if they differ). First run: 10 Nov 2026 on a
copy of the result tables (output 28c_beyond5kb.tsv and 28c.log in this folder).
