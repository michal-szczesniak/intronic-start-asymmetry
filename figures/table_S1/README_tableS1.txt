Supplemental Table S1 of the PURGE manuscript, built 11 Nov 2026

What this folder holds
  make_table_s1.py             builds the workbook (Supplemental Table S1 of the article) from the result tables;
                               every row names its result file and md5
  check_table_s1.py            checks the 116 numbers that the manuscript takes from the table against the workbook
                               cells, at the precision printed in the text (result: 116 agree, 0 differ)
  26c/26c_posthoc_rebuild.py   rebuild of step 26c (post hoc equal-weight methylation comparison); its outputs for
                               both methylation processings are in 26_methylation/results of the data archive

Sheets
  S1A Asymmetry                 depth bands and 1:1 matching; host-gene and single-pair resampling (23b, 02)
  S1B Promoter filter           distance from annotated starts (22c) and removal of pairs near them (22f)
  S1C Host activity, H3K36me3   host activity groups, trend, H3K36me3 within 1 kb (23b)
  S1D Share from sequence       shares under three calibrations, shuffles (23b); base composition (15)
  S1E Transposable elements     overlap, orientation, direction, removal of TE or repeat pairs (23c, 23b)
  S1F Genome-wide census        real vs shuffled introns, all and non-repeat positions (13c, 16)
  S1G Neighbours of silent genes three coverage rules (17b)
  S1H DNA methylation           step 26b strata (both processings) and step 26c contrasts, cells, rho
  S1I Population test           every class, set and expectation, both models (12, 19)
  The manuscript cites "Supplemental Table S1 (sheet S1H)" for the methylation analysis (Supplemental Fig. S4).

Provenance
  Inputs: 23_robustness/results/23d_figure_tables.tar.gz (md5 fb044b571501ff1deb44a5f00c7097dc) and
  26_methylation/results (26b_meth_primary_strata.tsv 489e9e17ed5f9b901c7df59b00e4525e,
  26b_meth_alternative_strata.tsv a8863ec55221e806e584f89e2c73eb49, gene tables 5e20e9b7996c19bbd06bdc34b56cb8a7
  and 08ba585e2473f4510a3bda3539e53708), as sent from the cluster on 11 Nov 2026.

  An earlier version of make_table_s1.py (md5 e4f2b94bef4dd621e195a64ade82e5f3, 23 Oct 2026) was rerun on the same
  23d archive and reproduced its workbook cell for cell (22,594 cells). The present version differs only in: README
  title; sheet S1H (DNA methylation) added; population test moved from S1H to S1I; the distance table of step 22c
  added to S1B; figure references in notes updated to the numbering of the manuscript. All values of the earlier
  sheets are unchanged (checked cell for cell).

  Step 26c. The original script (26c_weights_posthoc.py) and its tables were not kept; the primary log survived (sent 9 Nov 2026, md5 3a082807aeaa6052660f5be40d629757).
  26c_posthoc_rebuild.py recomputes the same quantities from the 26b gene tables (bootstrap 2,000, seed 1, genes
  resampled within each activity x methylation cell, cells drawn in the order A1-M1, A1-M2, ..., A3-M3). Its
  primary log is byte-identical to the surviving log (same md5), so every number of the original primary run is
  reproduced. The alternative processing gives the values cited in the Supplemental Fig. S4 legend (equal weight
  19.6, 95% CI -4.6 to 42.7).

Reproduce (in a folder holding the unpacked 23d archive as in23d/ and the 26 results as in26/26_methylation/results)
  python3 26c/26c_posthoc_rebuild.py --genes in26/26_methylation/results/26b_meth_primary_genes.tsv --set primary \
      --out-prefix in26/26_methylation/results/26c_posthoc --expect-md5 5e20e9b7996c19bbd06bdc34b56cb8a7 \
      --check-log <surviving 26c_posthoc_primary.log>
  python3 26c/26c_posthoc_rebuild.py --genes in26/26_methylation/results/26b_meth_alternative_genes.tsv --set alternative \
      --out-prefix in26/26_methylation/results/26c_posthoc --expect-md5 08ba585e2473f4510a3bda3539e53708
  python3 make_table_s1.py --data in23d --meth in26 --out Supplemental_Table_S1.xlsx
  python3 check_table_s1.py Supplemental_Table_S1.xlsx
  (Python 3, pandas, numpy, scipy, openpyxl; the xlsx file carries a creation time, so its md5 changes on every
  build while every cell stays the same.)
