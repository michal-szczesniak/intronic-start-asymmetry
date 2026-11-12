Figure 5 of the PURGE manuscript, 10 Nov 2026

Files
  fig5.py                        script that draws Figure 5 (schematic; no data are plotted)
  figstyle.py                    shared style of the PURGE figures
  fonts/LiberationSans-Regular.ttf, fonts/LiberationSans-Bold.ttf
                                 the font of the PURGE figures (Liberation 2.1, SIL Open Font License 1.1,
                                 fonts/LICENSE_LiberationFonts.txt); shipped because it is missing on many servers
  (the script writes the figure as 600-dpi PNG, PDF with embedded fonts and SVG with live text; the images
  are not included here)

Draw
  python3 fig5.py rebuilt/Szczesniak_Fig5
The script prints the matplotlib and FreeType versions and the font files it uses, stops if Liberation Sans
cannot be found (instead of switching silently to a wider font) and stops if any text leaves its band.
The figure in the manuscript was built with Python 3.13.18, matplotlib 3.11.2, FreeType 2.14.3, Liberation Sans 2.1 and
DejaVu Sans Mono from matplotlib (DNA letters). With the same versions the PNG and SVG are identical byte for
byte; the PDF bytes also depend on where the font file is found, and another FreeType version can change the
PNG bytes, in both cases without changing the figure.

Numbers quoted in the figure and where they come from:
  33% (K562), 21% (HCT116)       depth-matched bias of the sense start, Results 1, Fig. 1B
  30-42%                          share of the bias predicted by held-out sequence models, Results 2, Fig. 2B
  15-17%                          deficit of strong measured start sites in active genes, HCT116 and K562, Results 4 (step 28)
  46-56%                          the same deficit within 5 kb of the 5' end of the gene, Results 4, Supplemental Fig. S5 (step 28c)
  one hundredth                   upper limit relative to the excess of singletons at canonical splice sites, Results 5, Fig. 4C
  76,215 genomes                  gnomAD v4.1 genomes (Methods)
