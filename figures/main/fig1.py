#!/usr/bin/env python3
"""Figure 1: the question and the measured asymmetry. Data: source_data/fig1_asymmetry.csv (step 23b).

Rebuilt on 11 Nov 2026 from an earlier script and the image of the figure in the manuscript (the script of that
image was not kept); the output is identical to that image pixel for pixel."""
import sys

import numpy as np
import pandas as pd

from figstyle import (ANTI, GENE, GREY, HCT116, K562, LIGHT, MM, PALE, SENSE, WIDTH, arrow, box, ci_point,
                      fig_letter, plt, save)

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Figure1"

H = 112 * MM
fig = plt.figure(figsize=(WIDTH, H))

# ------------------------------------------------------------------ A: a divergent pair inside an intron
axA = fig.add_axes([0.0, 0.565, 1.0, 0.435])
axA.set_xlim(0, 100)
axA.set_ylim(0, 30)
axA.axis("off")
fig_letter(fig, "A", 0.005, 0.995)
yg = 24.0
axA.plot([5, 95], [yg, yg], color=GENE, lw=0.9, zorder=1)
for x0, x1 in ((5, 8.5), (22, 25), (68, 71), (89, 95)):
    box(axA, x0, yg - 1.3, x1 - x0, 2.6, fc=GENE, z=2)
arrow(axA, 5, yg + 1.3, 5, yg + 3.4, GENE, lw=0.9, style="-")
arrow(axA, 5, yg + 3.4, 9.0, yg + 3.4, GENE, lw=0.9, head=5)
axA.text(10, yg + 3.4, "host gene, transcribed left to right", va="center", fontsize=7, color=GENE)
xp = 46.5
axA.text(xp, yg + 2.0, "intron", ha="center", fontsize=6.5, color=GREY, style="italic")
arrow(axA, xp + 0.35, yg + 0.15, xp + 0.35, yg + 1.5, SENSE, lw=0.9, style="-")
arrow(axA, xp - 0.35, yg - 0.15, xp - 0.35, yg - 1.5, ANTI, lw=0.9, style="-")
axA.plot([xp - 1.2, 16], [yg - 1.9, 18.6], color=LIGHT, lw=0.7, ls="--", zorder=0)
axA.plot([xp + 1.2, 78], [yg - 1.9, 18.6], color=LIGHT, lw=0.7, ls="--", zorder=0)

yt, yb = 13.0, 9.6
for yy, l, r in ((yt, "5′", "3′"), (yb, "3′", "5′")):
    axA.plot([17, 77], [yy, yy], color=GENE, lw=0.8)
    axA.text(16.3, yy, l, ha="right", va="center", fontsize=6.5, color=GENE)
    axA.text(77.7, yy, r, ha="left", va="center", fontsize=6.5, color=GENE)
axA.text(17.3, yt + 0.7, "gene strand", ha="left", fontsize=6.5, color=GREY)
box(axA, 43.6, yb - 0.5, 7.0, yt - yb + 1.0, fc=PALE, ec=GREY, lw=0.5, r=0.5, z=0)
axA.text(51.2, yb - 0.9, "element, a few hundred nt", ha="left", va="top", fontsize=6.5, color=GREY)
# sense start with its reads (S) and RNA
xs_, xa_ = 49.2, 45.0
for k, h in enumerate([2.4, 1.8, 1.1, 0.6]):
    box(axA, xs_ + 0.15 + 0.45 * k, yt, 0.36, h, fc=SENSE, z=3)
axA.text(xs_ + 2.9, yt + 1.0, "S", color=SENSE, fontsize=8, fontweight="bold", ha="center")
arrow(axA, xs_, yt, xs_, yt + 4.2, SENSE, lw=1.4, style="-")
arrow(axA, xs_, yt + 4.2, 55.0, yt + 4.2, SENSE, lw=1.4, head=6)
xx = np.linspace(56, 76, 200)
axA.plot(xx, yt + 4.2 + 0.4 * np.sin(1.1 * (xx - 56)), color=SENSE, lw=0.9)
axA.text(50.5, yt + 6.0, "sense start", fontsize=7, color=SENSE, fontweight="bold", ha="left")
axA.text(79.0, yt + 4.2, "can use the gene’s own splice\nand polyadenylation signals", fontsize=6.5, va="center",
         ha="left", color=GENE, linespacing=1.15)
# antisense start with its reads (A) and RNA
for k, h in enumerate([0.7, 1.3, 2.1, 2.8, 2.0]):
    box(axA, xa_ - 0.5 - 0.45 * k, yb - h, 0.36, h, fc=ANTI, z=3)
axA.text(xa_ - 3.6, yb - 2.2, "A", color=ANTI, fontsize=8, fontweight="bold", ha="center")
arrow(axA, xa_, yb, xa_, yb - 4.2, ANTI, lw=1.4, style="-")
arrow(axA, xa_, yb - 4.2, 39.0, yb - 4.2, ANTI, lw=1.4, head=6)
xx = np.linspace(28, 38, 120)
axA.plot(xx, yb - 4.2 + 0.4 * np.sin(1.1 * (xx - 28)), color=ANTI, lw=0.9)
axA.text(43.5, yb - 6.4, "antisense start", fontsize=7, color=ANTI, fontweight="bold", ha="right")
axA.text(27.0, yb - 4.2, "short-lived, degraded\nin the nucleus", fontsize=6.5, va="center", ha="right", color=GENE,
         linespacing=1.15)
axA.text(80.5, yb - 2.8, "asymmetry of a pair\n= log2((S + 1) / (A + 1));\nreference: pairs outside\nprotein-coding genes, with\nthe same "
         "read depth", fontsize=6.5, va="center", ha="left", color="black", linespacing=1.15)

# ------------------------------------------------------------------ B: the measured asymmetry
d = pd.read_csv(f"{SD}/fig1_asymmetry.csv")
axB = fig.add_axes([0.085, 0.11, 0.40, 0.40])
fig_letter(fig, "B", 0.005, 0.545)
for line, col, dx in (("K562", K562, -0.13), ("HCT116", HCT116, 0.13)):
    b = d[(d["line"] == line) & (d["comparison"] == "depth band")].sort_values("stratum")
    for _, r in b.iterrows():
        ci_point(axB, r["stratum"] + 1 + dx, r["sense_weaker_pct"], r["pct_lo95"], r["pct_hi95"], col, ms=3.4)
    m = d[(d["line"] == line) & (d["comparison"] == "matched")].iloc[0]
    ci_point(axB, 6.7 + 2 * dx, m["sense_weaker_pct"], m["pct_lo95"], m["pct_hi95"], col, marker="D", ms=4.2, lw=1.2)
axB.axhline(0, color=GREY, lw=0.6, zorder=0)
axB.axvline(6.0, color=LIGHT, lw=0.6, ls=":", zorder=0)
axB.set_xticks([1, 2, 3, 4, 5, 6.7])
axB.set_xticklabels(["1", "2", "3", "4", "5", "all,\nmatched"])
axB.set_xlim(0.4, 7.3)
axB.set_ylim(-40, 100)
axB.set_yticks([-40, -20, 0, 20, 40, 60, 80, 100])
axB.set_xlabel("read depth of the pair, fifths (low → high)")
axB.set_ylabel("sense start weaker than in\nintergenic pairs (%, 95% CI)")
for line, col, yy in (("K562", K562, 0.97), ("HCT116", HCT116, 0.89)):
    m = d[(d["line"] == line) & (d["comparison"] == "matched")].iloc[0]
    axB.text(0.02, yy, f"{line}: {m['sense_weaker_pct']:.0f}% ({m['pct_lo95']:.0f}–{m['pct_hi95']:.0f}%), "
                       f"{int(m['n_pairs']):,} matched pairs", transform=axB.transAxes, ha="left", va="top",
             fontsize=6.5, color=col, fontweight="bold")

# ------------------------------------------------------------------ C: three explanations and their predictions
axC = fig.add_axes([0.555, 0.03, 0.44, 0.49])
axC.set_xlim(0, 10)
axC.set_ylim(0, 10)
axC.axis("off")
fig_letter(fig, "C", 0.535, 0.545)
cols = ["weaker sense\nstarts at active\nelements (Fig. 2)", "fewer potential\nsense starts in all\nintrons (Fig. 3)",
        "new sense starts\nkept rare in\npeople (Fig. 4)"]
rows = ["selection removes\nsense starts", "mutation shifts\nbase composition", "property of the\nelements"]
pred = [["yes", "yes", "yes"], ["maybe", "yes, from\ncomposition", "no"], ["yes", "not\nrequired", "not\nrequired"]]
x0, cw, rh, ytop = 3.0, 2.33, 1.62, 7.45
for j, c in enumerate(cols):
    axC.text(x0 + cw * (j + 0.5), 9.05, c, ha="center", va="center", fontsize=6.5, linespacing=1.1, color=GENE)
for i, r in enumerate(rows):
    yy = ytop - rh * i
    axC.text(0.0, yy, r, va="center", fontsize=6.8, linespacing=1.1)
    for j in range(3):
        v = pred[i][j]
        kind = v.split(",")[0]
        fc, ec, ls = {"yes": (LIGHT, "none", "-"), "no": ("white", GREY, "-"), "not\nrequired": ("white", GREY, "-"),
                      "maybe": ("white", GREY, "--")}.get(kind, (LIGHT, "none", "-"))
        box(axC, x0 + cw * j + 0.09, yy - rh / 2 + 0.12, cw - 0.18, rh - 0.24, fc=fc, ec=ec, lw=0.6, r=0.12, ls=ls)
        axC.text(x0 + cw * (j + 0.5), yy, v, ha="center", va="center", fontsize=6.8, linespacing=1.05)
# bottom row: what the tests found
axC.plot([0, 10], [3.21, 3.21], color=GENE, lw=0.6)
yy = 2.34
axC.text(0.0, yy, "found in this\nstudy", va="center", fontsize=6.8, linespacing=1.1, fontweight="bold")
for j, v in enumerate(["yes", "little or\nnone", "no"]):
    box(axC, x0 + cw * j + 0.09, yy - rh / 2 + 0.12, cw - 0.18, rh - 0.24, fc="white", ec="black", lw=1.0, r=0.12)
    axC.text(x0 + cw * (j + 0.5), yy, v, ha="center", va="center", fontsize=6.8, linespacing=1.05, fontweight="bold")
axC.text(0.0, 0.55, "Filled: predicted; open: not predicted or not required; dashed: possible.\n"
         "Bottom row: what the tests of Figures 2–4 found.", fontsize=6.3, color=GREY)

save(fig, OUT)
