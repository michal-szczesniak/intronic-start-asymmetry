#!/usr/bin/env python3
"""Supplemental Figure S3: in silent genes, the strand difference of measured single-read sites follows the strand
reached by transcription from active neighbours. Data: source_data/figS3_neighbours.csv (step 17b).

Drawing unchanged from an earlier script of 23 Oct 2026."""
import sys

import numpy as np
import pandas as pd

from figstyle import (ANTI, GENE, GREY, LIGHT, MID, MM, PROCAP, PUFFIN, SENSE, WIDTH, arrow, box, ci_point,
                      fig_letter, log2_ratio_axis, plt, save)

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS3"
H = 112 * MM
fig = plt.figure(figsize=(WIDTH, H))
L2 = np.log2

# ------------------------------------------------------------------ A: the coverage rule, to scale (kb)
axA = fig.add_axes([0.0, 0.62, 1.0, 0.37])
fig_letter(fig, "A", 0.005, 0.995)
axA.set_xlim(-60, 120)
axA.set_ylim(-1.6, 8.6)
axA.axis("off")
COV = "#BDBDBD"


def strand_pair(y, label):
    for dy in (0.45, -0.45):
        axA.plot([-12, 95], [y + dy, y + dy], color=LIGHT, lw=0.6, zorder=0)
    axA.text(-13, y + 0.45, "+", ha="right", va="center", fontsize=6.5, color=GREY)
    axA.text(-13, y - 0.45, "−", ha="right", va="center", fontsize=6.5, color=GREY)
    axA.text(-59, y, label, ha="left", va="center", fontsize=6.5, linespacing=1.15)


def active_gene(y):
    box(axA, 0, y + 0.25, 30, 0.4, fc=GENE, z=3)                    # active gene on +
    arrow(axA, 0, y + 0.65, 0, y + 1.25, GENE, lw=0.8, style="-")
    arrow(axA, 0, y + 1.25, 4, y + 1.25, GENE, lw=0.8, head=4)
    box(axA, 0, y + 0.18, 50, 0.54, fc=COV, z=1)                    # own strand to 20 kb past the end
    box(axA, -5, y - 0.72, 5, 0.54, fc=COV, z=1)                    # opposite strand, 5 kb upstream
    axA.text(4.6, y + 1.25, "active gene", ha="left", va="center", fontsize=6.3, color=GENE)


y1, y2 = 5.9, 2.0
strand_pair(y1, "silent gene reached\non its gene strand")
active_gene(y1)
box(axA, 40, y1 + 0.25, 30, 0.4, fc="none", ec=SENSE, lw=0.9, z=3)
arrow(axA, 40, y1 + 0.65, 40, y1 + 1.25, SENSE, lw=0.8, style="-")
arrow(axA, 40, y1 + 1.25, 44, y1 + 1.25, SENSE, lw=0.8, head=4)
axA.text(44.6, y1 + 1.25, "silent gene", ha="left", va="center", fontsize=6.3, color=SENSE)

strand_pair(y2, "silent gene reached\non its opposite strand")
active_gene(y2)
box(axA, 40, y2 + 0.25 - 0.9, 30, 0.4, fc="none", ec=ANTI, lw=0.9, z=3)
arrow(axA, 70, y2 - 0.65, 70, y2 - 1.25, ANTI, lw=0.8, style="-")
arrow(axA, 70, y2 - 1.25, 66, y2 - 1.25, ANTI, lw=0.8, head=4)
axA.text(65.4, y2 - 1.25, "silent gene", ha="right", va="center", fontsize=6.3, color=ANTI)

box(axA, -12, -1.05, 4, 0.54, fc=COV, z=1)
axA.text(-7, -0.78, "covered by the active gene: its own strand from the start to 20 kb past the end,\n"
                    "and the opposite strand for 5 kb upstream of the start", ha="left", va="center", fontsize=6.3,
         color=GENE, linespacing=1.15)
axA.plot([75, 95], [-0.9, -0.9], color=GENE, lw=0.8)
axA.text(85, -0.75, "20 kb", ha="center", va="bottom", fontsize=6.0, color=GENE)

# ------------------------------------------------------------------ B, C: strand ratios by class
nb = pd.read_csv(f"{SD}/figS3_neighbours.csv")
CLS = [("isolated", "not reached"), ("sense_only", "reached on the\ngene strand only"),
       ("antisense_only", "reached on the\nopposite strand only"), ("both", "reached on\nboth strands")]


def class_axis(ax, metrics, xl, ticks, xlabel, letter, lx, labels):
    for i, (b, lab) in enumerate(CLS):
        yy = len(CLS) - 1 - i
        for metric, col, mk, filled, dy in metrics:
            r = nb[(nb["bin"] == b) & (nb["metric"] == metric)].iloc[0]
            lo, hi = max(L2(r["ratio_lo95"]), xl[0]), min(L2(r["ratio_hi95"]), xl[1])
            ci_point(ax, L2(r["ratio"]), yy + dy, lo, hi, col, marker=mk, ms=3.2, horizontal=True, filled=filled)
            if L2(r["ratio_lo95"]) < xl[0]:
                ax.plot([xl[0] + 0.03 * (xl[1] - xl[0])], [yy + dy], marker="<", ms=3.0, color=col, clip_on=False)
            if L2(r["ratio_hi95"]) > xl[1]:
                ax.plot([xl[1] - 0.03 * (xl[1] - xl[0])], [yy + dy], marker=">", ms=3.0, color=col, clip_on=False)
        if labels:
            n = int(nb[nb["bin"] == b]["n_units"].iloc[0])
            ax.text(xl[0] - 0.04 * (xl[1] - xl[0]), yy, f"{lab}\n({n:,} genes)", ha="right", va="center",
                    fontsize=6.4, linespacing=1.1)
    ax.axvline(0, color=GREY, lw=0.6, zorder=0)
    ax.set_xlim(*xl)
    ax.set_ylim(-0.6, len(CLS) - 0.4)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    log2_ratio_axis(ax, ticks=ticks, axis="x")
    ax.set_xlabel(xlabel)
    fig_letter(fig, letter, lx, 0.53)


axB = fig.add_axes([0.20, 0.10, 0.33, 0.40])
class_axis(axB, [("obs r1", "#BDBDBD", "D", True, 0.22), ("obs r3", MID, "D", True, 0.0),
                 ("obs r10", "black", "D", True, -0.22)],
           (L2(0.25), L2(4)), (0.25, 0.5, 1, 2, 4), "measured start sites,\ngene strand / opposite strand", "B",
           0.005, True)
hand = [plt.Line2D([], [], marker="D", ls="none", mfc=c, mec=c, ms=3.2, label=l)
        for c, l in (("#BDBDBD", "≥1 read"), (MID, "≥3 reads"), ("black", "≥10 reads"))]
axB.legend(handles=hand, loc="upper right", handletextpad=0.2, borderaxespad=0.1)

axC = fig.add_axes([0.66, 0.10, 0.30, 0.40])
class_axis(axC, [("Puffin-D sum", PUFFIN, "s", True, 0.14), ("ProCapNet sum", PROCAP, "o", True, -0.14)],
           (L2(0.5), L2(2.5)), (0.5, 0.75, 1, 1.5, 2), "summed predicted signal,\ngene strand / opposite strand",
           "C", 0.57, False)
hand = [plt.Line2D([], [], marker="s", ls="none", mfc=PUFFIN, mec=PUFFIN, ms=3.2, label="Puffin-D"),
        plt.Line2D([], [], marker="o", ls="none", mfc=PROCAP, mec=PROCAP, ms=3.2, label="ProCapNet")]
axC.legend(handles=hand, loc="lower right", handletextpad=0.2, borderaxespad=0.1)
save(fig, OUT)
