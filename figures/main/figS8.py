#!/usr/bin/env python3
"""Supplemental Figure S8: the two models class different variants.
Data: source_data/figS8_class_overlap.csv (step 19c log; intronic variants of protein-coding genes, chr8 + chr9).

Drawing unchanged from an earlier script of 23 Oct 2026."""
import sys

import pandas as pd

from figstyle import GENE, MM, PROCAP, PUFFIN, WIDTH, box, plt, save

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS8"
H = 75 * MM
fig = plt.figure(figsize=(WIDTH, H))
ov = pd.read_csv(f"{SD}/figS8_class_overlap.csv").set_index("class")
ROWS = [("creates sense start (weak)", "creates a start, gene strand (weak site)"),
        ("creates antisense start (weak)", "creates a start, opposite strand (weak site)"),
        ("creates sense start (medium)", "creates a start, gene strand (medium site)"),
        ("removes sense start (weak)", "removes a start, gene strand (weak site)"),
        ("strong gain, sense", "strongly raises a start, gene strand"),
        ("strong gain, antisense", "strongly raises a start, opposite strand"),
        ("moderate gain, sense", "moderately raises a start, gene strand"),
        ("strong loss, sense", "strongly lowers a start, gene strand")]
BOTH = GENE


def pct(a, b):
    p = 100 * a / b
    return f"{p:.1f}%" if p < 1 else f"{p:.0f}%"


ax = fig.add_axes([0.30, 0.15, 0.50, 0.74])
for i, (c, lab) in enumerate(ROWS):
    yy = len(ROWS) - 1 - i
    nC, nP, nB = (int(ov.loc[c, k]) for k in ("procapnet", "puffin_d", "both"))
    box(ax, 0, yy - 0.3, nC - nB, 0.6, fc=PROCAP, z=2)
    box(ax, nC - nB, yy - 0.3, nB, 0.6, fc=BOTH, z=2)
    box(ax, nC, yy - 0.3, nP - nB, 0.6, fc=PUFFIN, z=2)
    ax.text(-3000, yy, lab, ha="right", va="center", fontsize=6.5)
    ax.text(nC + nP - nB + 3000, yy, f"{nB:,} in both ({pct(nB, nP)} of Puffin-D,\n{pct(nB, nC)} of ProCapNet)",
            ha="left", va="center", fontsize=6.0, color=GENE, linespacing=1.05)
ax.set_ylim(-0.6, len(ROWS) - 0.4)
ax.set_xlim(0, 190000)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks([0, 50000, 100000, 150000])
ax.set_xticklabels(["0", "50,000", "100,000", "150,000"])
ax.set_xlabel("intronic variants of protein-coding genes on chromosomes 8 and 9")
hand = [plt.Rectangle((0, 0), 1, 1, fc=PROCAP, label="ProCapNet only"),
        plt.Rectangle((0, 0), 1, 1, fc=BOTH, label="both models"),
        plt.Rectangle((0, 0), 1, 1, fc=PUFFIN, label="Puffin-D only")]
ax.legend(handles=hand, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, handletextpad=0.4, columnspacing=1.2,
          borderaxespad=0.2)
save(fig, OUT)
