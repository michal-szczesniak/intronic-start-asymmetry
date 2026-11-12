# Scripts

One line per script, from its own description. Steps are numbered in the order in which they were run; the
folders under `figures/` hold the figure scripts. Shell wrappers (`run_*.sh`) ran the steps on a SLURM cluster.

| Script | Description |
|---|---|
| `02_istart_baseline/02_istart_baseline.py` | Step 02: reproduce the iSTART directionality numbers before any new analysis. |
| `02_istart_baseline/02b_istart_contrasts.py` | Step 02b: does the intronic sense/antisense asymmetry depend on host activity? |
| `02_istart_baseline/02c_depth_matched.py` | Step 02c: is the intronic DI deficit an artefact of lower signal depth? |
| `02_istart_baseline/03_arm_coordinates.py` | Step 03: recover per-arm coordinates for the iSTART divergent pairs. |
| `02_istart_baseline/04_arm_qc.py` | Step 04: quality control of the recovered arm coordinates, before any modelling. |
| `05_env/05_setup_env.sh` | Step 05: build a self-contained environment for PURGE and install Puffin / Puffin-D. |
| `06_puffin_smoke/06_gpu_smoke_test.sh` | Step 06: check that Puffin-D runs on a GPU here, and measure how long one prediction takes. |
| `07_puffin_pairs/07_puffin_pairs.py` | Step 07: does Puffin-D reproduce the iSTART strand asymmetry from reference sequence alone? |
| `07_puffin_pairs/07b_calibrate.py` | Step 07b: how much of the observed strand asymmetry does Puffin-D reproduce? |
| `07_puffin_pairs/07c_puffin_pairs_rc.py` | Step 07c: step 07 again, with strand-symmetric (reverse-complement averaged) predictions. |
| `08_tracks/08_predict_tracks.py` | Step 08: strand-symmetric Puffin-D PRO-cap predictions along whole chromosomes. |
| `09_h1_introns/09_h1_intron_strand.py` | Step 09: H1 on whole introns - is predicted initiation weaker on the sense strand? |
| `09_h1_introns/09b_h1_diagnostics.py` | Step 09b: is the intronic sense-strand deficit of step 09 robust? |
| `10_shuffle/10_shuffle_test.py` | Step 10: does the sense-strand deficit survive a composition-preserving shuffle? |
| `11_gnomad/11a_check_gnomad_access.sh` | Step 11a: can this node reach gnomAD v4.1 genome VCFs, how fast, and with which tools? |
| `11_gnomad/11c_mux_test.py` | Step 11c: can many gnomAD variants be scored in one Puffin-D pass, and how fast? |
| `11_gnomad/11d_check_scores.py` | Step 11d check: what a scores table from 11d_score_variants.py contains. |
| `11_gnomad/11d_score_variants.py` | Step 11d: score gnomAD SNVs for their effect on predicted initiation (both strands). |
| `12_selection/12_selection_test.py` | Step 12: are variants that create a start on the gene's strand selected against? |
| `12_selection/12b_selection_robustness.py` | Step 12b: is the step-12 result robust to mutation-rate heterogeneity beyond the trinucleotide? |
| `12_selection/12c_selection_genomewide.py` | Step 12c: the selection test of steps 12/12b on many chromosomes, read one at a time. |
| `13_genomewide/13c_shuffle_combine.py` | Step 13c: genome-wide H1 map and shuffle test, from the per-chromosome step-10 tables of 13b. |
| `13_genomewide/13d_gc_strata.py` | Step 13d: does the sense deficit of potential start sites depend on regional composition? |
| `14_host/14_host_transcription.py` | Step 14: does the part of the sense weakness that sequence does not explain depend on |
| `15_composition/15_composition_calib.py` | Step 15: is the part of the sense weakness that Puffin-D does not explain carried by |
| `16_procapnet/16a_get_procapnet.py` | Step 16a: fetch the ProCapNet K562 models from ENCODE and check that they load. |
| `16_procapnet/16b_score_pairs_procapnet.py` | Step 16b: score the iSTART pairs with ProCapNet, each pair by a model that never saw its |
| `16_procapnet/16c_procapnet_tracks.py` | Step 16c: genome-wide ProCapNet tracks for one chromosome, real and shuffled, from the model |
| `16_procapnet/16d_procapnet_sites.py` | Step 16d: potential start sites in the ProCapNet tracks of 16c (real and shuffled), counted |
| `17_observed/17_observed_sites.py` | Step 17: observed K562 PRO-cap initiation in introns, sense vs antisense, as the arbiter |
| `17_observed/17b_neighbors.py` | Step 17b: is the observed sense deficit in silent genes real initiation, or transcription |
| `19_procapnet_gnomad/19a_score_variants_procapnet.py` | Step 19a: score the gnomAD SNVs of step 11d with ProCapNet, using the model that never saw |
| `19_procapnet_gnomad/19b_selection_classes.py` | Step 19b: step 12c with the class definitions as parameters, so that the variant classes can come |
| `19_procapnet_gnomad/19c_match_cutoffs.py` | Step 19c: change cut-offs for the ProCapNet variant classes, matched to Puffin-D's class sizes. |
| `19_procapnet_gnomad/patch_12c_to_19b.py` | Make 19b_selection_classes.py from 12c_selection_genomewide.py: the same test with the class |
| `20_pair_shuffle/20a_pairs_shuffle_procapnet.py` | Step 20a: where in the sequence of an iSTART pair does ProCapNet find the weaker sense start? |
| `22_open_items/22c_tss_distance.py` | Step 22c: distance from the iSTART pairs to the nearest annotated transcription start. |
| `22_open_items/22f_tss_filter.py` | Step 22f: write iSTART pair tables without pairs that sit at an annotated transcription start. |
| `23_robustness/23b_pairs.py` | Step 23b (pairs): the depth-matched asymmetry with confidence intervals that resample host genes, |
| `23_robustness/23b_shares.py` | Step 23b (shares): share of the intronic asymmetry predicted from sequence, on the same pairs as |
| `23_robustness/23c_transposons.py` | Step 23c: do transposable elements explain the weaker sense start of intronic pairs? |
| `26_methylation/26b_methylation.py` | Step 26b: gene-body DNA methylation and the strand deficit of measured start sites (K562). |
| `26_methylation/26c_posthoc_rebuild.py` | Step 26c (rebuild, 11 Nov 2026): post hoc weighting of the DNA-methylation comparison of step 26b. |
| `27a_start_motifs/27a_start_motifs.py` | Step 27a: sequence around the two starts of each pair - logos and pre-specified core-promoter elements. |
| `27a_start_motifs/run_27a.sh` | Step 27a (run on a second server, outside the cluster): sequence around the two starts of every depth-matched pair (PURGE step 23b) - base counts for |
| `28_k36_position/28_download.py` | Step 28: download the ENCODE files (H3K36me3 peaks of both lines, HCT116 PRO-cap bigWigs), check them against |
| `28_k36_position/28_k36_position.py` | Step 28: H3K36me3, position along the gene and a second cell line for the gene-strand deficit of strong |
| `28_k36_position/run_28.sh` | Step 28: H3K36me3, position along the gene and HCT116 for the gene-strand deficit of strong measured start sites. |
| `28_k36_position/tests/make_synth.py` | Synthetic genome, annotation, PRO-cap bigWigs (2 replicates per strand) and H3K36me3 peaks for testing step 28. |
| `28_k36_position/tests/naive_check.py` | Independent, loop-based recount of step-28 strata for the synthetic data; compares with the script's gene table. |
| `28b_position/28b_position.py` | Step 28b (EXPLORATORY, specified on 10 Nov 2026 after the step-28 results): where along the host gene lie |
| `28b_position/run_28b.sh` | Step 28b (EXPLORATORY): position along the host gene of the pair bias (Results 1) and of the strong-site deficit |
| `28b_position/tests/make_pairs.py` | Synthetic matched pairs and arm table for testing step 28b: pairs inside protein-coding genes of the synthetic |
| `28b_position/tests/make_synth_alt.py` | Copy synth.gtf.gz to synth2.gtf.gz, adding an alternative transcript that starts in the middle of every third |
| `28b_position/tests/naive_check_28b.py` | Independent loop-based recount for step 28b on the synthetic data (run in the synthetic test folder). |
| `28c_beyond5kb/28c_beyond5kb.py` | Step 28c (EXPLORATORY, 10 Nov 2026, after step 28b): bias and deficit beyond the first 5 kb of the gene. |
| `28d_cpg/28d_cpg.py` | Step 28d (EXPLORATORY; specified on 10 Nov 2026 after step 28c, before running): are CpG islands behind the |
| `28d_cpg/28d_download.py` | Step 28d: download the UCSC CpG-island table for hg38 (cpgIslandExt.txt.gz) with retries after server errors, |
| `28d_cpg/run_28d.sh` | Step 28d (EXPLORATORY): CpG islands and the strand difference within 5 kb of the gene's start. |
| `28d_cpg/tests/make_cgi.py` | Synthetic CpG-island table (UCSC cpgIslandExt format) for the synthetic genome of step 28: islands next to the |
| `28d_cpg/tests/naive_check_28d.py` | Independent loop-based check of step 28d: distances to islands, site counts and pair estimates for the subset |
| `28e_twobytwo/28e_twobytwo.py` | Step 28e (EXPLORATORY; specified on 10 Nov 2026 after step 28d, before running): near the gene's start versus near |
| `28e_twobytwo/tests/naive_check_28e.py` | Independent check of step 28e point estimates (direct pandas computation, no bootstrap). |
| `28f_activity/28f_activity.py` | Step 28f (EXPLORATORY; specified on 10 Nov 2026 after step 28e, before running): does the near-start effect grow |
| `28f_activity/tests/naive_check_28f.py` | Independent check of step 28f point estimates (direct computation). Usage: naive_check_28f.py <28d dir> <28 dir> <28f prefix> |
| `28g_metagene/28g_metagene.py` | Step 28g (12 Nov 2026; descriptive, for Figure 3F): strong measured start sites on the two strands of active |
| `28g_metagene/run_28g.sh` | Step 28g: strong measured sites in 1-kb bins from the 5' end of active genes (Figure 3F), from the site tables of |
| `28g_metagene/tests/test_28g.py` | Synthetic test of 28g_metagene.py: random genes and sites; counts per bin, genes reaching each bin and the |
| `29_grocap/29a_discover.sh` | Step 29a (version b: only file links named GSM..._ are kept; version a also took a link to an HHS web page) |
| `29_grocap/29a_inspect.py` | Step 29a: format check of the GRO-cap files (GEO GSM1480321 with TAP, GSM1480322 without TAP; hg19) and of |
| `29_grocap/29b_liftover.py` | Step 29b: lift single-base GRO-cap counts (bigWig, one value per base, as checked in step 29a) from hg19 to |
| `29_grocap/29b_summary.py` | Step 29b summary: the three pre-specified questions for K562 GRO-cap, next to the same estimates from K562 |
| `29_grocap/run_29b.sh` | Step 29b: K562 GRO-cap (GSM1480321, with TAP; hg19) lifted over to hg38, then the unchanged step-28 and step-28b |
| `29_grocap/tests/test_29a.py` | Synthetic test of 29a_inspect.py: writes small bigWig, bedGraph and chain files with known content into the |
| `29_grocap/tests/test_29b_liftover.py` | Synthetic test of 29b_liftover.py against an independent position-by-position mapper. |
| `figures/fig5/fig5.py` | Figure 5: the conclusions at three scales. Schematic; no data are plotted. |
| `figures/fig5/figstyle.py` | Shared style for the PURGE figures: sans serif 6.5-8 pt, bold 10-pt panel letters, |
| `figures/figS_motifs/figS_motifs.py` | Supplemental figure: sequence at the two starts of the depth-matched pairs and the pre-specified core-promoter |
| `figures/figS_motifs/figstyle.py` | Shared style for the PURGE figures: sans serif 6.5-8 pt, bold 10-pt panel letters, |
| `figures/figS_position/figS_position.py` | Supplemental figure: where along the gene the strand difference lies (steps 28, 28b, 28e, 28f and 29). |
| `figures/figS_position/figstyle.py` | Shared style for the PURGE figures: sans serif 6.5-8 pt, bold 10-pt panel letters, |
| `figures/main/check_figures.py` | Compares each 600-dpi figure in figures/ with the image of the same figure in the manuscript |
| `figures/main/fig1.py` | Figure 1: the question and the measured asymmetry. Data: source_data/fig1_asymmetry.csv (step 23b). |
| `figures/main/fig2.py` | Figure 2: host activity, the share predicted from sequence and where it lies. |
| `figures/main/fig3.py` | Figure 3: potential and measured start sites on the two strands of introns (version 3). |
| `figures/main/fig4.py` | Figure 4: no sign of selection against new sense start sites. |
| `figures/main/figS1.py` | Supplemental Figure S1: which model never saw which chromosome, and the share of the asymmetry predicted by |
| `figures/main/figS3.py` | Supplemental Figure S3: in silent genes, the strand difference of measured single-read sites follows the strand |
| `figures/main/figS4.py` | Supplemental Figure S4: gene-body DNA methylation and the strand deficit of strong |
| `figures/main/figS6.py` | Supplemental Figure S6: strand ratios of start sites by the keto skew of the gene strand. |
| `figures/main/figS7.py` | Supplemental Figure S7: the strand comparison of the population test under four settings of the expectation. |
| `figures/main/figS8.py` | Supplemental Figure S8: the two models class different variants. |
| `figures/main/figstyle.py` | Shared style for the PURGE figures: sans serif 6-8 pt, bold 10-pt panel letters, open axes, Okabe-Ito colours |
| `figures/main/graphical_abstract.py` | Graphical abstract (5:2, drawn at 127 x 50.8 mm). Qualitative summary; no data are plotted |
| `figures/main/make_source_data.py` | Source Data for the figures, read directly from the result tables of the analysis steps. |
| `figures/table_S1/26c/26c_posthoc_rebuild.py` | Step 26c (rebuild, 11 Nov 2026): post hoc weighting of the DNA-methylation comparison of step 26b. |
| `figures/table_S1/check_table_s1.py` | Checks every number that the manuscript (main text, Supplemental Results and Methods, legends) takes from |
| `figures/table_S1/make_table_s1.py` | Supplemental Table S1 (11 Nov 2026): robustness checks, the DNA methylation |
| `supplemental_data/build_supplemental_data.py` | Supplemental Data of the PURGE manuscript: every iSTART start-site pair of K562 and HCT116 |
