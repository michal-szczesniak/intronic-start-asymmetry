#!/usr/bin/env python3
"""Supplemental Figure S7: the strand comparison of the population test under four settings of the expectation.
Data: source_data/fig4_classes.csv and fig4_differences.csv (steps 12c and 19b).

Drawing unchanged from an earlier script of 23 Oct 2026."""
import sys

import pandas as pd

from figstyle import GENE, GREY, LIGHT, MM, PALE, PROCAP, PUFFIN, WIDTH, ci_point, fig_letter, plt, save

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS7"
H = 140 * MM
fig = plt.figure(figsize=(WIDTH, H))
cl = pd.read_csv(f"{SD}/fig4_classes.csv")
df = pd.read_csv(f"{SD}/fig4_differences.csv")
THR = {"Puffin-D": ("0.0439", "0.1432"), "ProCapNet": ("1.131", "2.658")}
NAME = {"creates": lambda m, s: f"creates_{s}_{THR[m][0]}", "creates_medium": lambda m, s: f"creates_{s}_{THR[m][1]}",
        "strong": lambda m, s: f"gain_{s}_strong", "moderate": lambda m, s: f"gain_{s}_mid",
        "removes": lambda m, s: f"destroys_{s}_{THR[m][0]}"}
DIFF = [("creates", "intron_protein_coding", "creates a start (weak site)"),
        ("creates_medium", "intron_protein_coding", "creates a start (medium site)"),
        ("strong", "intron_protein_coding", "strongly raises a start"),
        ("moderate", "intron_protein_coding", "moderately raises a start"),
        ("removes", "intron_protein_coding", "removes a start"),
        ("creates", "intergenic", "intergenic: creates a start\n(plus strand − minus strand)")]
SETTINGS = [("base", "three-base context"), ("noCpG", "without CpG sites"), ("gc5", "context × GC bin"),
            ("noCpG_gc5", "both")]
XL = (-0.016, 0.012)
STEP = len(SETTINGS) + 1.2


def panel(ax, model, col, mk, labels):
    can = cl[(cl["model"] == model) & (cl["class"] == "splice_canonical") & (cl["setting"] == "base")]["excess"].iloc[0]
    ax.axvspan(-0.01 * can, 0.01 * can, color=LIGHT, zorder=0, lw=0)
    ytop = len(DIFF) * STEP
    for i, (eff, st, lab) in enumerate(DIFF):
        y0 = ytop - i * STEP
        if i % 2 == 0:
            ax.axhspan(y0 - len(SETTINGS) - 0.1, y0 + 0.5, color=PALE, zorder=-1, lw=0)
        if labels:
            ax.text(XL[0] - 0.0122, y0 - (len(SETTINGS) - 1) / 2, lab, ha="right", va="center", fontsize=6.5,
                    linespacing=1.1)
        for j, (sett, slab) in enumerate(SETTINGS):
            yy = y0 - j
            r = df[(df["model"] == model) & (df["setting"] == sett) & (df["set"] == st)
                   & (df["sense_class"] == NAME[eff](model, "sense"))].iloc[0]
            lo, hi = max(r["diff_lo95"], XL[0]), min(r["diff_hi95"], XL[1])
            ci_point(ax, r["excess_diff"], yy, lo, hi, col, marker=mk, ms=3.0, horizontal=True, filled=(j == 0))
            if r["diff_lo95"] < XL[0]:
                ax.plot([XL[0] + 0.0004], [yy], marker="<", ms=2.8, color=col, zorder=4, clip_on=False)
            if r["diff_hi95"] > XL[1]:
                ax.plot([XL[1] - 0.0004], [yy], marker=">", ms=2.8, color=col, zorder=4, clip_on=False)
            if labels:
                ax.text(XL[0] - 0.0007, yy, slab, ha="right", va="center", fontsize=6.0, color=GENE)
    ax.axvline(0, color=GREY, lw=0.6, zorder=0)
    ax.set_xlim(*XL)
    ax.set_ylim(ytop - (len(DIFF) - 1) * STEP - len(SETTINGS) + 0.2, ytop + 0.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-0.01, 0, 0.01])
    ax.set_xlabel("excess of singletons,\ngene strand minus opposite strand (95% CI)")
    ax.set_title(f"classes from {model}", fontsize=7, color=col, fontweight="bold")


axA = fig.add_axes([0.375, 0.085, 0.27, 0.87])
fig_letter(fig, "A", 0.005, 0.995)
panel(axA, "ProCapNet", PROCAP, "o", True)
axB = fig.add_axes([0.715, 0.085, 0.27, 0.87])
fig_letter(fig, "B", 0.665, 0.995)
panel(axB, "Puffin-D", PUFFIN, "s", False)
save(fig, OUT)
