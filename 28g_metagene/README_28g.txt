PURGE step 28g (12 Nov 2026; descriptive, for Figure 3F): strong measured start sites on the two strands of active
protein-coding genes in 1-kb bins of distance from the annotated 5' end of the gene (1-30 kb), with the strand ratio
in 1-kb bins and in 3-kb windows (95% CI, 2,000 resamples of genes within activity thirds, as in step 28b) and the
density per 1,000 active genes whose annotated span reaches the bin. Drawn after the analyses of distance of step
28b (exploratory); no test.

No new site calling: inputs are the site tables of step 28b (28b_K562_sites.tsv.gz, 28b_HCT116_sites.tsv.gz) and
of step 29b (29_GROcap_28b_sites.tsv.gz: the step-28b script on K562 GRO-cap), and the gene tables of steps 28 and
29b (activity and span of every gene). md5 of the inputs in run_28g.sh.
Checks in the script (STOP otherwise): number of active genes; all strong sites of active genes equal step 28
(K562 1,224 / 1,483; HCT116 927 / 1,092; GRO-cap 2,450 / 2,771, sense / antisense); sites below 5 kb equal the
"1-5" stratum of step 28b (232 / 431; 111 / 254; 412 / 667). All passed.

Files
  28g_metagene.py          the analysis
  tests/test_28g.py        synthetic genes and sites; counts, genes reaching each bin and bootstrap ratios recounted
                           by explicit loops (232 values, 0 differ); a wrong expected count must stop the script
  run_28g.sh               md5 of the inputs, synthetic test, the three runs, comparison with results/
  results/28g_{K562,HCT116,GROcap}_bins.tsv and .log

How the results were made: on 12 Nov 2026 the six input tables were copied from the cluster
(metagene_inputs_2026-11-12.tar.gz, md5 bfffaa8a323721bebc4b47327e29fc5c; every file matched the md5 lists of
steps 28b and 29b and run_28b.sh) and the script was run outside the cluster on these copies (Python 3.13, pandas 3.0.5,
numpy 2.5.3), after the synthetic test. run_28g.sh repeats this on the cluster and compares the tables.
