#!/usr/bin/env python3
"""Graphical abstract (5:2, drawn at 127 x 50.8 mm). Qualitative summary; no data are plotted
except the range 30-42% (share of the asymmetry predicted by ProCapNet models that never saw the chromosome,
Figure 2B). Panel 3 shows the two tests of selection: start sites on both strands across introns (Figure 3) and
single-base variants in people (Figure 4).
Saved also as SVG with live text, for editing in a vector program.

Rebuilt on 11 Nov 2026 from an earlier script and the image of the graphical abstract (the script of that image
was not kept); the output is identical to that image pixel for pixel."""
import sys

import numpy as np

from figstyle import ANTI, GENE, GREY, MM, PROCAP, SENSE, arrow, box, plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "figures/GraphicalAbstract"
W, H = 127 * MM, 50.8 * MM
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 127)
ax.set_ylim(0, 50.8)
ax.axis("off")
HEAD, BODY = 12, 12
BG = "#F3F3F3"
COLS = [(1.0, 39.5), (43.8, 39.5), (86.5, 39.5)]
for x0, w in COLS:
    box(ax, x0, 1.0, w, 48.8, fc=BG, ec="none", r=2.0, z=0)


def head(x0, w, txt):
    ax.text(x0 + w / 2, 47.5, txt, ha="center", va="top", fontsize=HEAD, fontweight="bold", linespacing=1.0,
            color="#222222")


def foot(x0, w, txt):
    ax.text(x0 + w / 2, 3.0, txt, ha="center", va="bottom", fontsize=BODY, linespacing=1.0, color="#222222")


# ------------------------------------------------------------------ 1: weaker sense starts inside genes
x0, w = COLS[0]
head(x0, w, "Sense starts\nare weaker")
yg = 24.5
ax.plot([x0 + 2.5, x0 + w - 2.5], [yg, yg], color=GENE, lw=1.2, solid_capstyle="round")
for a_, b_ in ((x0 + 2.5, x0 + 6.0), (x0 + w - 6.0, x0 + w - 2.5)):
    box(ax, a_, yg - 1.4, b_ - a_, 2.8, fc=GENE, r=0.4, z=2)
arrow(ax, x0 + 2.5, yg + 1.4, x0 + 2.5, yg + 3.6, GENE, lw=1.1, style="-")
arrow(ax, x0 + 2.5, yg + 3.6, x0 + 7.0, yg + 3.6, GENE, lw=1.1, head=5)
ax.text(x0 + 7.6, yg + 3.6, "gene", va="center", fontsize=9, color=GENE)
xm = x0 + w / 2
box(ax, xm - 3.2, yg - 0.9, 6.4, 1.8, fc="white", ec=GREY, lw=0.6, r=0.6, z=3)
arrow(ax, xm + 1.2, yg + 0.9, xm + 1.2, yg + 3.2, SENSE, lw=1.0, style="-", z=4)
arrow(ax, xm + 1.2, yg + 3.2, xm + 6.0, yg + 3.2, SENSE, lw=1.0, head=4, z=4)
ax.text(xm + 6.6, yg + 3.2, "sense", va="center", fontsize=9, color=SENSE)
arrow(ax, xm - 1.2, yg - 0.9, xm - 1.2, yg - 4.6, ANTI, lw=3.0, style="-", z=4)
arrow(ax, xm - 1.2, yg - 4.6, xm - 9.0, yg - 4.6, ANTI, lw=3.0, head=7, z=4)
ax.text(xm - 5.0, yg - 6.6, "antisense", va="top", ha="center", fontsize=9, color=ANTI)
foot(x0, w, "inside genes, by a\nfifth to a third")

# ------------------------------------------------------------------ 2: partly encoded in the element
x0, w = COLS[1]
head(x0, w, "Partly predicted\nfrom DNA")
yg = 29.0
rng = np.random.default_rng(5)
TONES = ["#4D4D4D", "#8C8C8C", "#B3B3B3", "#6E6E6E"]
for xx in np.linspace(x0 + 4, x0 + 16, 20):
    ax.plot([xx, xx], [yg - 2.0, yg + 2.0], color=TONES[rng.integers(0, 4)], lw=1.5, solid_capstyle="butt")
arrow(ax, x0 + 17.5, yg, x0 + 21.5, yg, GENE, lw=1.1, head=5)
box(ax, x0 + 22.5, yg - 2.6, 13.5, 5.2, fc=PROCAP, r=0.9, z=2)
ax.text(x0 + 29.25, yg, "model", ha="center", va="center", fontsize=10, color="white", fontweight="bold")
ax.text(x0 + w / 2, 19.5, "30–42%", ha="center", va="center", fontsize=16, fontweight="bold", color=PROCAP)
foot(x0, w, "of the difference\nis predicted\nfrom the DNA")

# ------------------------------------------------------------------ 3: no genome-wide purging
x0, w = COLS[2]
head(x0, w, "No sign\nof selection")
yg = 28.0
for yy, col, d in ((yg + 4.0, SENSE, 1), (yg - 0.5, ANTI, -1)):
    ax.plot([x0 + 3.5, x0 + w - 3.5], [yy, yy], color=GENE, lw=0.9)
    for xx in np.linspace(x0 + 8, x0 + w - 9, 4):
        arrow(ax, xx, yy, xx, yy + 1.8 * d, col, lw=1.0, style="-")
        arrow(ax, xx, yy + 1.8 * d, xx + 2.6 * d, yy + 1.8 * d, col, lw=1.0, head=4)
# a stretch of intron with two single-base variants (the human-variant test of Figure 4)
SEQ, ALT = "GACTTAGCATCGA", {3: "C", 9: "G"}
YS, YV, BW, BH = 15.4, 20.05, 2.3, 3.1
for i, (xx, b) in enumerate(zip(np.linspace(x0 + 5, x0 + w - 5, len(SEQ)), SEQ)):
    if i in ALT:
        box(ax, xx - BW / 2, YS - BH / 2, BW, BH, fc="white", ec=GENE, lw=0.7, r=0.4, z=3)
        ax.text(xx, YS, b, ha="center", va="center", family="monospace", fontsize=8.5, fontweight="bold", color=GENE,
                zorder=4)
        ax.plot([xx, xx], [YS + BH / 2, YV - BH / 2], color=GENE, lw=0.7, zorder=2)
        box(ax, xx - BW / 2, YV - BH / 2, BW, BH, fc=GENE, r=0.4, z=3)
        ax.text(xx, YV, ALT[i], ha="center", va="center", family="monospace", fontsize=8.5, fontweight="bold",
                color="white", zorder=4)
    else:
        ax.text(xx, YS, b, ha="center", va="center", family="monospace", fontsize=8.5, color=GREY)
foot(x0, w, "in all introns or\nin human variants")

fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf", metadata={"CreationDate": None})
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["svg.hashsalt"] = "purge"
fig.savefig(OUT + ".svg", metadata={"Date": None})
plt.close(fig)
