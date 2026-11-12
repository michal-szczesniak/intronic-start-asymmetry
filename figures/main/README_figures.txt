Figures 1-4, the graphical abstract and Supplemental Figs. S1, S3, S4, S6, S7 and S8 of the PURGE manuscript, with
the scripts and Source Data that draw them.

Files
  fig1.py fig2.py fig3.py fig4.py   Figures 1-4
  graphical_abstract.py             graphical abstract (PNG, PDF and SVG with live text)
  figS1.py figS3.py figS4.py figS6.py figS7.py figS8.py
                                    Supplemental Figs. S1, S3, S4, S6, S7 and S8
  figstyle.py                       shared style; registers the fonts in fonts/, stops if Liberation Sans is missing
                                    (instead of falling back silently to a wider font), prints the versions and font
                                    files, writes PDFs without a creation date
  fonts/                            Liberation Sans 2.1 (SIL Open Font License 1.1)
  make_source_data.py               writes source_data/ from the result tables; every CSV row names its result file
                                    and md5
  source_data/                      Source Data of all figures above
  check_figures.py                  compares newly drawn figures with the images of the figures in the manuscript
                                    (manuscript_images/, not included here)

Figure 3 (version 3, 12 Nov 2026)
  A  the DNA that each model reads, drawn to scale (schematic)
  B  potential start sites: 2 models x (introns, original DNA; introns, DNA shuffled; far from genes) x (weak,
     medium, strong); ratio and 95% CI in each cell; colour, log2 ratio (blue: fewer sites on the gene strand)
  C  measured K562 PRO-cap by host activity: all reads summed and sites with >= 1, 3 and 10 reads
  D  the two models by host activity, summed signal (top) and strong sites (bottom), with the measured curves of C
     as a faint reference
  E  strong measured sites by set of genes, with HCT116 PRO-cap (>= 24 reads) and K562 GRO-cap (>= 14 reads) in
     active protein-coding genes and intergenic DNA
  F  strong measured sites by distance from the 5' end of active genes in 1-kb bins (1-30 kb): mirrored densities
     for K562 PRO-cap (sites per 1,000 genes whose span reaches the bin) and the strand ratio in 3-kb windows with
     95% CI for the three data sets (step 28g)

Supplemental Fig. S4, data
  A  gene-body methylation of protein-coding genes with >= 20 CpGs by host activity (box: quartiles; whiskers:
     10th-90th percentile), ENCODE default processing   <- 26b_meth_primary_genes.tsv
  B  strand ratio (gene / opposite strand) of strong measured start sites (>= 10 reads) by methylation third, inside
     each activity third, pooled (pre-specified) and in intergenic 1-Mb blocks; 95% CIs  <- 26b_meth_primary_strata.tsv
  C  share of the opposite-strand strong sites that each activity third contributes to each methylation third
     <- 26b_meth_primary_strata.tsv (anti counts; 729, 528, 210 sites)
  D  high-minus-low methylation contrast of the deficit, pooled over activity thirds (pre-specified; step 26b) and
     with equal weight for each activity third (post hoc; step 26c), for both processings
     <- 26b_meth_{primary,alternative}_strata.tsv and 26c_posthoc_{primary,alternative}_contrasts.tsv

Inputs (paths relative to the project folder; all are in the data archive of the Zenodo record)
  steps 02-23: the result tables of these steps (the set of archive 23_robustness/results/23d_figure_tables.tar.gz,
               md5 fb044b571501ff1deb44a5f00c7097dc, as used here)
  step 26:     26_methylation/results/26b_meth_primary_genes.tsv       5e20e9b7996c19bbd06bdc34b56cb8a7
               26_methylation/results/26b_meth_primary_strata.tsv      489e9e17ed5f9b901c7df59b00e4525e
               26_methylation/results/26b_meth_alternative_strata.tsv  a8863ec55221e806e584f89e2c73eb49
               26_methylation/results/26c_posthoc_primary_contrasts.tsv      991170f570246c94e5f957973fac7430
               26_methylation/results/26c_posthoc_alternative_contrasts.tsv  a57ea151227757d7eb671ac71b1916f7
               (26c files: rebuild of 11 Nov 2026, see figures/table_S1/README_tableS1.txt)
  steps 28/29: figures/figS_position/data (copies of the strata tables of steps 28 and 29; see that folder)
               28_HCT116_strata.tsv       82f53f3dbcccd96709260750ccffaa2b  (Figure 3E, HCT116 PRO-cap)
               29_K562_GROcap_strata.tsv  15a7e68bbc6c418b55ef6188bcc13dc6  (Figure 3E, K562 GRO-cap)
               -> fig3_replication.csv (rows with metric "strong": >= 24 and 14 reads)
  step 28g:    28g_metagene/results/28g_{K562,HCT116,GROcap}_bins.tsv -> fig3_metagene.csv
  log-derived: fig2_folds.csv (fa398b7581cdec9f145f7d7361fe9313) and fig5_class_overlap.csv
               (dc393c9e8ed502090c70d226b49d53d2), as in 18_manuscript/source_data; they feed figS1_folds.csv and
               figS8_class_overlap.csv (marked "[log-derived; no result table]")

Reproduce (<data>: the folder into which the data archive of the Zenodo record was extracted, e.g. the top folder
of this repository; writes to rebuilt_sd/ and rebuilt/)
  python3 make_source_data.py --data <data> --meth <data> --old-sd <data>/18_manuscript/source_data \
      --position <data>/figures/figS_position/data --metagene <data>/28g_metagene/results --out rebuilt_sd
  cmp -s rebuilt_sd/figS4_contrasts.csv source_data/figS4_contrasts.csv && echo same    (likewise for each CSV)
  mkdir -p rebuilt
  for f in fig1 fig2 fig3 fig4 figS1 figS3 figS4 figS6 figS7 figS8; do python3 $f.py rebuilt_sd rebuilt/$(grep -o 'figures/Figure[A-Z0-9]*' $f.py | cut -d/ -f2); done
  python3 graphical_abstract.py rebuilt/GraphicalAbstract
  (Python 3 with pandas, numpy, matplotlib and Pillow)

Verification (12 Nov 2026)
  Every figure drawn by these scripts is identical, pixel for pixel, to the image of the figure in the manuscript
  (check_figures.py: 11 figures, 0 differing pixels each). Built with Python 3.13.18, matplotlib 3.11.2, FreeType
  2.14.3, pandas 3.0.5, numpy 2.5.3, Pillow 12.3.0, Liberation Sans 2.1 (fonts/) and DejaVu Sans Mono from
  matplotlib (DNA letters of the graphical abstract). With the same versions the PNG, PDF and SVG files are
  identical byte for byte from run to run; another FreeType or matplotlib version can change the PNG bytes and, in
  small ways, the drawing.
