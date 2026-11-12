Supplemental figure on where along the gene the strand difference lies (steps 28, 28b, 28e, 28f and 29),
PURGE manuscript, 10 Nov 2026: Supplemental Fig. S5.

Files
  figS_position.py        script (draws the figure and writes the Source Data; checks every input by md5)
  figstyle.py             shared style of the PURGE figures
  fonts/                  Liberation Sans 2.1 (SIL Open Font License 1.1)
  data/                   result tables of the steps, unchanged (in the data archive of the Zenodo record):
    28_K562_strata.tsv    e3b5e9eeb25d487a1d0b346eef5ca664  step 28  (28_results.tar.gz)
    28_HCT116_strata.tsv  82f53f3dbcccd96709260750ccffaa2b  step 28
    28b_K562_strata.tsv   8ecf4540c3f3bd108ae910a2ea006bb7  step 28b (28b_results.tar.gz)
    28b_HCT116_strata.tsv fbaed9ebac905eec7189d67a27ed2a11  step 28b
    28e_cells.tsv         3ae8b9c9c48a6b44ee253465ef2bd68f  step 28e (purge28e_2026-11-10.tar.gz)
    28f_cells.tsv         430a0d9ff54b9918fe5af451292fb925  step 28f (purge28f_2026-11-10.tar.gz)
    29_K562_GROcap_strata.tsv 15a7e68bbc6c418b55ef6188bcc13dc6  step 29b (29b_results.tar.gz)
    29_GROcap_28b_strata.tsv  7978f868ca2a43260e1aa7537df5430a  step 29b
  figures/Szczesniak_FigS_position_source_data.csv   Source Data, panels A-E (estimates, 95% CIs, P, n, source table)
  (the script also writes the figure as 600-dpi PNG, PDF and SVG; the images are not included here)

Panels
  A  deficit of strong measured sites in fifths of active protein-coding genes (steps 28 and 29, pre-specified)
  B  the same by distance from the annotated 5' end of the host gene (step 28b, exploratory; step 29, K562
     GRO-cap, pre-specified)
  C  depth-matched bias of the start-site pairs by distance from the 5' end (step 28b, exploratory)
  D  both lines pooled: within 5 kb of the 5' end or farther x within 2 kb of a CpG island or farther (step 28e)
  E  both lines pooled: excess within 5 kb by third of host activity (step 28f)

Rebuild:  python3 figS_position.py data rebuilt/Szczesniak_FigS_position
A rebuild gives byte-identical .png, .pdf, .svg and .csv (checked with matplotlib 3.11.2, FreeType 2.14.3).
The script stops if an input differs or Liberation Sans cannot be found.
Panel B shows the same sites as Figure 3F in coarser bins and over the whole length of the genes.
