Supplemental figure on core promoter elements at the two starts of the depth-matched pairs (step 27a),
PURGE manuscript, 10 Nov 2026: Supplemental Fig. S2.

Files
  figS_motifs.py                        script (draws the figure and writes the Source Data)
  figstyle.py                           shared style of the PURGE figures
  fonts/                                Liberation Sans 2.1 (SIL Open Font License 1.1)
  data/27a_pfm.tsv, data/27a_tests.tsv  step-27a result tables (checked by md5 in the script; in the data archive of the Zenodo record)
  figures/Szczesniak_FigS_motifs_logos.csv     Source Data, panel A (base counts and information content)
  figures/Szczesniak_FigS_motifs_elements.csv  Source Data, panel B (differences, CIs, P, Holm P, n)
  (the script also writes the figure as 600-dpi PNG, PDF and SVG; the images are not included here)

Rebuild:  python3 figS_motifs.py data rebuilt/Szczesniak_FigS_motifs
Inputs come from step 27a: script 27a_start_motifs.py (md5 cbe348734af56180b9bd2a48db278095), 27a_pfm.tsv
(bde79d95af096729398db01dac9b356e), 27a_tests.tsv (09407c3544bb97cfb84e2c13b687cd11). The script stops if
the inputs differ or Liberation Sans cannot be found. Built with matplotlib 3.11.2 and FreeType 2.14.3.
