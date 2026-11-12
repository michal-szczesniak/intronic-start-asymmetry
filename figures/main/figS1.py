#!/usr/bin/env python3
"""Supplemental Figure S1: which model never saw which chromosome, and the share of the asymmetry predicted by
models that did or did not see it, under three calibrations.
Data: source_data/figS1_folds.csv (step 16b log; Dudnyk et al. 2024) and fig2_shares.csv (step 23b).

Drawing unchanged from an earlier script of 23 Oct 2026."""
import sys

import numpy as np
import pandas as pd

from figstyle import GENE, GREY, HCT116, K562, LIGHT, MM, PALE, WIDTH, box, ci_point, fig_letter, plt, save

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS1"
H = 120 * MM
fig = plt.figure(figsize=(WIDTH, H))

# ------------------------------------------------------------------ A: chromosomes and folds
fo = pd.read_csv(f"{SD}/figS1_folds.csv").set_index("chrom")
fo = fo.loc[[f"chr{i}" for i in range(1, 23)]]
mb = (fo["length_bp"] / 1e6).to_numpy()
axA = fig.add_axes([0.06, 0.60, 0.92, 0.38])
fig_letter(fig, "A", 0.005, 0.995)
axA.set_xlim(-0.8, 23.9)
axA.set_ylim(-80, 290)
axA.axis("off")
unseen = [i for i, r in enumerate(fo["puffin_role"]) if r in ("test", "validation")]
for i, (c, r) in enumerate(fo.iterrows()):
    box(axA, i - 0.32, 0, 0.64, mb[i], fc=LIGHT, ec="none", r=0.2)
    axA.text(i, -10, c[3:], ha="center", va="top", fontsize=6.5)
    axA.text(i, mb[i] + 5, str(int(r["procapnet_test_fold"])), ha="center", va="bottom", fontsize=6.8,
             fontweight="bold", color=GENE)
for i in unseen:
    box(axA, i - 0.42, -2, 0.84, mb[i] + 4, fc="none", ec="black", lw=0.9, r=0.25, z=4)
axA.plot([min(unseen) - 0.42, max(unseen) + 0.42], [-42, -42], color="black", lw=0.8)
test = ", ".join(c[3:] for c in fo.index[fo["puffin_role"] == "test"])
val = ", ".join(c[3:] for c in fo.index[fo["puffin_role"] == "validation"])
axA.text(np.mean(unseen), -47, f"never seen by Puffin-D\n(chromosomes {test}: test; {val}: validation)",
         ha="center", va="top", fontsize=6.5, linespacing=1.1)
axA.text(10.4, 285, "Number above each chromosome: the ProCapNet model (one of seven, numbered 0–6)\n"
                    "trained without it. A held-out prediction for a start-site pair comes from that model;\n"
                    "a train-fold prediction is the mean of the six models that saw the chromosome.\n"
                    "The K562 and Caco-2 models use the same groups of chromosomes.",
         ha="left", va="top", fontsize=6.5, color=GENE, linespacing=1.2)
axA.plot([22.75, 22.75], [0, 100], color=GENE, lw=0.8)
for y_ in (0, 100):
    axA.plot([22.75, 22.9], [y_, y_], color=GENE, lw=0.6)
axA.text(23.0, 50, "100\nMb", ha="left", va="center", fontsize=6.0, color=GENE, linespacing=1.1)

# ------------------------------------------------------------------ B: shares under three calibrations
sh = pd.read_csv(f"{SD}/fig2_shares.csv")
sh = sh[sh["group"] == "intron_active"]
ROWS = [("puffin_D", "Puffin-D", "s"),
        ("pc_K562_heldout", "ProCapNet K562,\nheld-out", "o"),
        ("pc_K562_trainfold", "ProCapNet K562,\ntrain-fold", "o"),
        ("pc_Caco2_heldout", "ProCapNet Caco-2,\nheld-out", "o"),
        ("pc_Caco2_trainfold", "ProCapNet Caco-2,\ntrain-fold", "o")]
CAL = [("frac_linear", "calibrated on intergenic pairs\n(as in Figure 2B)"),
       ("frac_rank", "rank-based\n(no scale conversion)"),
       ("frac_depth", "calibrated within\nthirds of read depth")]
fig_letter(fig, "B", 0.005, 0.50)
x0, w, gap = 0.20, 0.245, 0.02
for k, (col_, title) in enumerate(CAL):
    ax = fig.add_axes([x0 + k * (w + gap), 0.09, w, 0.36])
    for i, (lab, name, mk) in enumerate(ROWS):
        yy = len(ROWS) - 1 - i
        if "trainfold" in lab:
            ax.axhspan(yy - 0.5, yy + 0.5, color=PALE, zorder=-1, lw=0)
        for line, col, dy in (("K562", K562, -0.15), ("HCT116", HCT116, 0.15)):
            r = sh[(sh["label"] == lab) & (sh["line"] == line)].iloc[0]
            ci_point(ax, 100 * r[col_], yy - dy, 100 * r[col_ + "_lo95"], 100 * r[col_ + "_hi95"], col, marker=mk,
                     ms=3.2, horizontal=True)
        if k == 0:
            ax.text(-38, yy, name, ha="right", va="center", fontsize=6.5, linespacing=1.1)
    ax.axvline(0, color=GREY, lw=0.6, zorder=0)
    ax.axvline(100, color=LIGHT, lw=0.6, ls=":", zorder=0)
    ax.set_xlim(-35, 125)
    ax.set_ylim(-0.5, len(ROWS) - 0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 50, 100])
    ax.set_title(title, fontsize=6.5, linespacing=1.1, color=GENE)
    if k == 1:
        ax.set_xlabel("share of the asymmetry predicted from sequence (%, 95% CI), active hosts")
    if k == 2:
        for line, col, yy in (("K562", K562, 0.16), ("HCT116", HCT116, 0.06)):
            ax.text(0.98, yy, line, transform=ax.transAxes, color=col, fontsize=6.8, fontweight="bold", ha="right")
save(fig, OUT)
