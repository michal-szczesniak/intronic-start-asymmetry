# Intronic divergent start sites: strand asymmetry, sequence and selection

Code for the analyses and figures of:

> Szcześniak MW, Wanowska E. No evidence of genome-wide selection against sense transcription starts in human introns. Manuscript submitted for publication (2026).

The analyses compare the two starts of divergent start-site pairs inside introns of human protein-coding genes with pairs outside genes (PRO-cap, K562 and HCT116), ask how much of the weaker start in the gene's direction two sequence models (Puffin-D and ProCapNet) predict from the DNA alone, count potential and measured start sites on the two strands of introns genome-wide, and test for purifying selection against variants that create start sites (gnomAD v4.1).

## Contents

- `02_istart_baseline/` ... `29_grocap/`: the scripts of each analysis step (02–29), numbered in the order in which they were run; `tests/` hold the synthetic-data checks run before each real run, and `README_*.txt` describe the steps that have one. `SCRIPTS.md` gives a one-line description of every script.
- `figures/`: scripts and Source Data of every figure: `main/` (Figures 1–4, graphical abstract, Supplemental Figs. S1, S3, S4, S6–S8), `fig5/` (Figure 5), `figS_motifs/` (Supplemental Fig. S2), `figS_position/` (Supplemental Fig. S5) and `table_S1/` (Supplemental Table S1). Each folder has a README.
- `supplemental_data/`: the script that writes the Supplemental Data table (every start-site pair with its per-pair measures).
- `env/`: the conda environment of the analyses (`purge.yml`), the installed Python packages and the Python version.
- The result tables of every step are archived in Zenodo ([10.5281/zenodo.24161109](https://doi.org/10.5281/zenodo.24161109)).

## Data

All input data are public; accession numbers, versions and download dates are given in the Methods of the article. Scripts that download data record the md5 checksum of every file and stop if it differs. The start-site pairs of step 02 come from a catalogue of divergent start sites described in a separate manuscript; the pair tables used here are in the data archive (`02_istart_baseline/results`).

## Running

The steps ran on a SLURM cluster with the environment in `env/purge.yml` (Python 3.11; `conda env create -f env/purge.yml`); steps that score DNA with Puffin-D or ProCapNet need a GPU. Python scripts take their inputs and outputs as arguments (see the header of each script). In the shell wrappers (`run_*.sh` and a few `*.sh` steps), the absolute paths of the cluster are replaced by two variables: set `PURGE_PROJECT` (the project folder, which holds the step folders and their `results/`) and `PURGE_SCRATCH` (a scratch folder for the conda installation, the genome, downloaded files and temporary files) before running them; step 27a uses `PURGE27_DIR`.

The figures can be redrawn from the archived tables without the raw data: extract the data archive of the Zenodo record into this folder and follow `figures/main/README_figures.txt` (and the README of the other figure folders). The figures of the article were drawn with Python 3.13, matplotlib 3.11.2 and the Liberation Sans fonts included in `figures/*/fonts`.

## License

Code: MIT License (see `LICENSE`). Fonts in `figures/*/fonts`: SIL Open Font License 1.1 (`LICENSE_LiberationFonts.txt`). Data in the Zenodo archive: CC BY 4.0.

## Citation

Please cite the article above and the Zenodo archive ([10.5281/zenodo.24161109](https://doi.org/10.5281/zenodo.24161109)).
