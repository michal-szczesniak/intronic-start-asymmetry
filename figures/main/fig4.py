#!/usr/bin/env python3
"""Figure 4: no sign of selection against new sense start sites.
Data: source_data/fig4_classes.csv and fig4_differences.csv (steps 12c and 19b; setting 'base').

Rebuilt on 11 Nov 2026 from an earlier script and the image of the figure in the manuscript (the script of that
image was not kept); the output is identical to that image pixel for pixel."""
import sys

import pandas as pd

from figstyle import (ANTI, GENE, GREY, LIGHT, MM, PROCAP, PUFFIN, SENSE, WIDTH, arrow, box, ci_point, fig_letter,
                      plt, save)

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Figure4"
H = 125 * MM
fig = plt.figure(figsize=(WIDTH, H))
cl = pd.read_csv(f"{SD}/fig4_classes.csv")
cl = cl[cl["setting"] == "base"]
df = pd.read_csv(f"{SD}/fig4_differences.csv")
df = df[df["setting"] == "base"]
THR = {"Puffin-D": ("0.0439", "0.1432"), "ProCapNet": ("1.131", "2.658")}
MODELS = (("ProCapNet", PROCAP, "o", 0.14), ("Puffin-D", PUFFIN, "s", -0.14))


def cls_name(model, effect, strand):
    w, m = THR[model]
    return {"creates": f"creates_{strand}_{w}", "creates_medium": f"creates_{strand}_{m}",
            "strong": f"gain_{strand}_strong", "moderate": f"gain_{strand}_mid", "removes": f"destroys_{strand}_{w}"}[effect]


# ------------------------------------------------------------------ A: logic of the test
axA = fig.add_axes([0.0, 0.84, 1.0, 0.16])
axA.set_xlim(0, 100)
axA.set_ylim(0, 10)
axA.axis("off")
fig_letter(fig, "A", 0.005, 0.995)
STEPS = ["Each variant in an intron is classified by its\npredicted effect on initiation on each strand\n"
         "(creates, raises or removes a start)",
         "Selection keeps harmful variants rare:\nmore of them are seen in only one\nof 76,215 genomes (singletons)",
         "Excess of singletons = observed share minus the\nshare expected from variants without effect\n"
         "(same mutation type and three-base context)"]
for i, t in enumerate(STEPS):
    x0 = 4 + 32.2 * i
    box(axA, x0, 1.0, 29.0, 8.0, fc="white", ec=GREY, lw=0.6, r=0.6)
    axA.text(x0 + 14.5, 5.0, t, ha="center", va="center", fontsize=6.5, linespacing=1.15)
    if i < 2:
        arrow(axA, x0 + 29.4, 5.0, x0 + 31.8, 5.0, GENE, lw=0.9, head=4)

# ------------------------------------------------------------------ B: excess of singletons by class
ROWS = [("creates", "sense", "intron_protein_coding", "creates a start, gene strand", SENSE),
        ("creates", "anti", "intron_protein_coding", "creates a start, opposite strand", ANTI),
        ("strong", "sense", "intron_protein_coding", "strongly raises, gene strand", SENSE),
        ("strong", "anti", "intron_protein_coding", "strongly raises, opposite strand", ANTI),
        ("moderate", "sense", "intron_protein_coding", "moderately raises, gene strand", SENSE),
        ("moderate", "anti", "intron_protein_coding", "moderately raises, opposite strand", ANTI),
        ("removes", "sense", "intron_protein_coding", "removes a start, gene strand", SENSE),
        ("removes", "anti", "intron_protein_coding", "removes a start, opposite strand", ANTI),
        ("creates", "sense", "intergenic", "intergenic: creates, plus strand", GREY),
        ("creates", "anti", "intergenic", "intergenic: creates, minus strand", GREY)]
axBc = fig.add_axes([0.25, 0.66, 0.24, 0.10])
fig_letter(fig, "B", 0.005, 0.80)
for j, (c, lab) in enumerate((("splice_canonical", "canonical splice sites"), ("splice_donor_3to6", "splice donor +3 to +6"))):
    yy = 1 - j
    r = cl[(cl["model"] == "ProCapNet") & (cl["class"] == c)].iloc[0]
    ci_point(axBc, r["excess"], yy, r["excess_lo95"], r["excess_hi95"], "black", marker="D", ms=3.0, horizontal=True)
    axBc.text(-0.006, yy, lab, ha="right", va="center", fontsize=6.4)
axBc.set_xlim(0, 0.16)
axBc.set_ylim(-0.6, 1.6)
axBc.set_yticks([])
axBc.spines["left"].set_visible(False)
axBc.set_xticks([0, 0.05, 0.10, 0.15])
axBc.tick_params(axis="x", labelsize=6)
axBc.text(0.16, 1.55, "positive controls", ha="right", va="bottom", fontsize=6.3, color=GREY)

axB = fig.add_axes([0.25, 0.09, 0.24, 0.52])
for i, (eff, strand, st, lab, col) in enumerate(ROWS):
    yy = len(ROWS) - 1 - i + (-0.5 if st == "intergenic" else 0)
    for model, mcol, mk, dy in MODELS:
        r = cl[(cl["model"] == model) & (cl["set"] == st) & (cl["class"] == cls_name(model, eff, strand))].iloc[0]
        ci_point(axB, r["excess"], yy + dy, r["excess_lo95"], r["excess_hi95"], col, marker=mk, ms=3.2,
                 horizontal=True, filled=(model == "ProCapNet"))
    axB.text(-0.0215, yy, lab, ha="right", va="center", fontsize=6.4)
axB.axvline(0, color=GREY, lw=0.6, zorder=0)
axB.set_xlim(-0.02, 0.02)
axB.set_ylim(-1.2, len(ROWS) - 0.4)
axB.set_yticks([])
axB.spines["left"].set_visible(False)
axB.set_xticks([-0.02, -0.01, 0, 0.01, 0.02])
axB.set_xlabel("excess of singletons (95% CI)")
hand = [plt.Line2D([], [], marker="o", ls="none", mfc=GREY, mec=GREY, ms=3.4, label="classes from ProCapNet"),
        plt.Line2D([], [], marker="s", ls="none", mfc="white", mec=GREY, ms=3.4, label="classes from Puffin-D")]
fig.legend(handles=hand, loc="lower left", bbox_to_anchor=(0.005, 0.0), ncol=1, handletextpad=0.2, frameon=False)

# ------------------------------------------------------------------ C: sense minus antisense
DIFF = [("creates", "intron_protein_coding", "creates a start (weak site)"),
        ("creates_medium", "intron_protein_coding", "creates a start (medium site)"),
        ("strong", "intron_protein_coding", "strongly raises a start"),
        ("moderate", "intron_protein_coding", "moderately raises a start"),
        ("removes", "intron_protein_coding", "removes a start"),
        ("creates", "intergenic", "intergenic: creates a start\n(plus strand − minus strand)")]
can = cl[(cl["model"] == "ProCapNet") & (cl["class"] == "splice_canonical")]["excess"].iloc[0]
axC = fig.add_axes([0.75, 0.09, 0.235, 0.60])
fig_letter(fig, "C", 0.53, 0.80)
XC0, XC1 = -0.015, 0.012
for i, (eff, st, lab) in enumerate(DIFF):
    yy = len(DIFF) - 1 - i
    for model, mcol, mk, dy in MODELS:
        r = df[(df["model"] == model) & (df["set"] == st) & (df["sense_class"] == cls_name(model, eff, "sense"))].iloc[0]
        lo, hi = max(r["diff_lo95"], XC0), min(r["diff_hi95"], XC1)
        ci_point(axC, r["excess_diff"], yy + dy, lo, hi, mcol, marker=mk, ms=3.2, horizontal=True)
        if r["diff_lo95"] < XC0:
            axC.plot([XC0 + 0.0003], [yy + dy], marker="<", ms=3.0, color=mcol, zorder=4, clip_on=False)
        if r["diff_hi95"] > XC1:
            axC.plot([XC1 - 0.0003], [yy + dy], marker=">", ms=3.0, color=mcol, zorder=4, clip_on=False)
    axC.text(XC0 - 0.0015, yy, lab, ha="right", va="center", fontsize=6.4, linespacing=1.1)
axC.axvspan(-0.01 * can, 0.01 * can, color=LIGHT, zorder=0, lw=0)
axC.axvline(0, color=GREY, lw=0.6, zorder=0)
axC.set_xlim(XC0, XC1)
axC.set_ylim(-0.6, len(DIFF) - 0.4)
axC.set_yticks([])
axC.spines["left"].set_visible(False)
axC.set_xticks([-0.01, 0, 0.01])
axC.set_xlabel("excess of singletons,\ngene strand minus opposite strand")
axC.text(XC1, len(DIFF) - 0.35, "gene-strand variants\nkept rarer →", ha="right", va="bottom", fontsize=6.0,
         color=GREY, linespacing=1.1)
axC.text(XC0, len(DIFF) - 0.35, "← opposite-strand\nvariants kept rarer", ha="left", va="bottom", fontsize=6.0,
         color=GREY, linespacing=1.1)
hand = [plt.Line2D([], [], marker="o", ls="none", mfc=PROCAP, mec=PROCAP, ms=3.4, label="classes from ProCapNet"),
        plt.Line2D([], [], marker="s", ls="none", mfc=PUFFIN, mec=PUFFIN, ms=3.4, label="classes from Puffin-D"),
        plt.Rectangle((0, 0), 1, 1, fc=LIGHT, ec="none", label="±1% of the excess at\ncanonical splice sites"),
        plt.Line2D([], [], marker=">", ls="none", mfc=GREY, mec=GREY, ms=3.0, label="interval continues")]
fig.legend(handles=hand, loc="upper left", bbox_to_anchor=(0.555, 0.80), ncol=1, handletextpad=0.4,
           borderaxespad=0.0, frameon=False)
save(fig, OUT)
