#!/usr/bin/env python3
"""Supplemental Figure S4: gene-body DNA methylation and the strand deficit of strong
measured start sites in K562 (steps 26b and 26c).
Data: source_data/figS4_meth_by_activity.csv (A), figS4_strand_ratio.csv (B, C) and figS4_contrasts.csv (D).

Rebuilt on 11 Nov 2026 from the image of an earlier figure, whose script was not kept; the rebuilt figure is
identical to that image pixel for pixel (checked at the size of the manuscript image). The positions of the key of
panel D are the measured positions of that image (axes fraction).
"""
import sys

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

from figstyle import GENE, GREY, LIGHT, MID, MM, PALE, WIDTH, ci_point, fig_letter, log2_ratio_axis, plt, save

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS4"
fig = plt.figure(figsize=(WIDTH, 112 * MM))
L2 = np.log2
DARK = "#4D4D4D"
BINS = ["M1", "M2", "M3"]                      # methylation thirds, formed inside each activity third
st = pd.read_csv(f"{SD}/figS4_strand_ratio.csv")


def cell(sub, b):
    return st[(st["subset"] == sub) & (st["bin"] == b)].iloc[0]


# ------------------------------------------------------------------ A: methylation by host activity
bx = pd.read_csv(f"{SD}/figS4_meth_by_activity.csv")
axA = fig.add_axes([0.075, 0.60, 0.25, 0.35])
fig_letter(fig, "A", 0.005, 0.99)
LAB = {"silent": "silent\n", "low": "weakly\nactive", "A1": "active,\nleast", "A2": "active,\nmiddle", "A3": "active,\nmost"}
stats = [dict(med=r["median"], q1=r["q1"], q3=r["q3"], whislo=r["p10"], whishi=r["p90"],
              label=f"{LAB[r['group']]}\n({int(r['n_genes']):,})") for _, r in bx.iterrows()]
axA.bxp(stats, positions=range(len(stats)), widths=0.6, showfliers=False, showcaps=False, patch_artist=True,
        boxprops=dict(facecolor=LIGHT, edgecolor=GENE, lw=0.6), medianprops=dict(color="black", lw=1.2),
        whiskerprops=dict(color=GENE, lw=0.7, zorder=1.5))  # whiskers: 10th-90th percentile, behind the boxes
axA.axvspan(1.5, 4.5, color=PALE, zorder=0, lw=0)
axA.tick_params(axis="x", labelsize=6.0)
axA.set_xlim(-0.6, 4.6)
axA.set_ylim(0, 1)
axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axA.set_yticklabels(["0.00", "0.25", "0.50", "0.75", "1.00"])
axA.set_ylabel("gene-body methylation\n(fraction of methylated reads)")
axA.text(1.0, 1.02, "genes in brackets", transform=axA.transAxes, ha="right", va="bottom", fontsize=5.8, color=GENE)

# ------------------------------------------------------------------ B: strand ratio of strong sites by methylation third
axB = fig.add_axes([0.43, 0.60, 0.555, 0.35])
fig_letter(fig, "B", 0.36, 0.99)
GROUPS = [("A1", "active genes,\nleast active third"), ("A2", "active genes,\nmiddle third"),
          ("A3", "active genes,\nmost active third"), ("pc_active", "all active genes\n(pre-specified)"),
          ("intergenic", "intergenic DNA\n(plus / minus strand)")]
tr = blended_transform_factory(axB.transData, axB.transAxes)
for g, (sub, lab) in enumerate(GROUPS):
    col = GREY if sub == "intergenic" else "black"
    xs_, ys_ = [], []
    for j, b in enumerate(BINS):
        r = cell(sub, b)
        ci_point(axB, 4 * g + j, L2(r["ratio"]), L2(r["ratio_lo95"]), L2(r["ratio_hi95"]), col, marker="D", ms=3.2)
        xs_.append(4 * g + j)
        ys_.append(L2(r["ratio"]))
    axB.plot(xs_, ys_, color=col, lw=0.7, zorder=2)
    axB.text(4 * g + 1, -0.16, lab, transform=tr, ha="center", va="top", fontsize=6.2)
axB.text(-0.8, -0.075, "methylation:", transform=tr, ha="right", va="top", fontsize=6.0)
axB.axhline(0, color=GREY, lw=0.6, zorder=0)
axB.axvspan(11.4, 14.6, color=PALE, zorder=-1, lw=0)
axB.set_xticks([4 * g + j for g in range(len(GROUPS)) for j in range(3)])
axB.set_xticklabels(["low", "mid", "high"] * len(GROUPS))
axB.tick_params(axis="x", labelsize=6.0)
axB.set_xlim(-0.8, 18.8)
log2_ratio_axis(axB, ticks=(0.25, 0.5, 0.75, 1, 1.5, 2))
axB.set_ylim(L2(0.23), L2(2.1))
axB.set_ylabel("strong start sites (≥10 reads),\ngene strand / opposite strand")

# ------------------------------------------------------------------ C: where the opposite-strand strong sites come from
axC = fig.add_axes([0.13, 0.10, 0.30, 0.25])
fig_letter(fig, "C", 0.005, 0.44)
COLS = {"A1": LIGHT, "A2": MID, "A3": DARK}
for i, m in enumerate(BINS):
    y = len(BINS) - 1 - i
    tot = sum(cell(a, m)["anti"] for a in ("A1", "A2", "A3"))
    left = 0.0
    for a in ("A1", "A2", "A3"):
        w = 100 * cell(a, m)["anti"] / tot
        axC.barh(y, w - 0.4, left=left + 0.2, height=0.62, color=COLS[a], lw=0)  # white gaps of 0.4 percentage points
        axC.text(left + w / 2, y, f"{w:.0f}%", ha="center", va="center", fontsize=6.0,
                 color="white" if a == "A3" else "black")
        left += w
    axC.text(102, y, f"{int(tot):,} sites", ha="left", va="center", fontsize=6.0, color=GENE)
axC.set_yticks([2, 1, 0])
axC.set_yticklabels(["low methylation", "mid methylation", "high methylation"])
axC.tick_params(axis="y", length=0)
axC.spines["left"].set_visible(False)
axC.set_xlim(0, 100)
axC.set_ylim(-0.6, 2.9)
axC.set_xticks([0, 25, 50, 75, 100])
axC.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
axC.set_xlabel("share of the strong sites on the opposite strand")
for x_, a, lab in ((0, "A1", "least active third"), (37, "A2", "middle third"), (66, "A3", "most active third")):
    axC.add_patch(Rectangle((x_, 2.66 - 0.11), 4, 0.22, color=COLS[a], lw=0, clip_on=False))
    axC.text(x_ + 5, 2.66, lab, ha="left", va="center", fontsize=6.0)

# ------------------------------------------------------------------ D: the high-minus-low contrast, two weightings
ct = pd.read_csv(f"{SD}/figS4_contrasts.csv")
axD = fig.add_axes([0.64, 0.10, 0.345, 0.25])
fig_letter(fig, "D", 0.535, 0.44)
for y, comp in ((1, "pooled over activity thirds (pre-specified)"), (0, "activity thirds weighted equally (post hoc)")):
    for proc, dy, filled in (("primary", 0.14, True), ("alternative", -0.14, False)):
        r = ct[(ct["processing"] == proc) & (ct["comparison"] == comp)].iloc[0]
        ci_point(axD, r["value"], y + dy, r["lo95"], r["hi95"], "black", marker="D", ms=3.2, horizontal=True,
                 filled=filled)
axD.axvline(0, color=GREY, lw=0.6, zorder=0)
axD.set_yticks([1, 0])
axD.set_yticklabels(["pooled over activity thirds\n(pre-specified)", "activity thirds weighted\nequally (post hoc)"],
                    fontsize=6.0)
axD.tick_params(axis="y", length=0)
axD.spines["left"].set_visible(False)
axD.set_xlim(-25, 55)
axD.set_xticks([-20, 0, 20, 40])
axD.set_ylim(-0.6, 1.9)
axD.set_xlabel("deficit, high minus low methylation\n(percentage points)")
for (mx, tx), filled, lab in (((0.0646, 0.10567), True, "ENCODE default processing"),
                              ((0.606, 0.6467), False, "earlier processing")):
    axD.plot([mx], [0.96445], marker="D", ms=3.2, mfc="black" if filled else "white", mec="black", mew=0.8,
             ls="none", transform=axD.transAxes, clip_on=False)
    axD.text(tx, 0.96445, lab, ha="left", va="center_baseline", fontsize=6.0, transform=axD.transAxes)
save(fig, OUT)
