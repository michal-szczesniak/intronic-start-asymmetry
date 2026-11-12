#!/usr/bin/env python3
"""Supplemental Figure S6: strand ratios of start sites by the keto skew of the gene strand.
Data: source_data/figS6_keto.csv (steps 13d, 16 and 17).

Drawing unchanged from an earlier script of 23 Oct 2026."""
import sys

import numpy as np
import pandas as pd

from figstyle import GENE, GREY, MM, PROCAP, PUFFIN, WIDTH, ci_point, fig_letter, log2_ratio_axis, plt, save

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/FigureS6"
H = 70 * MM
fig = plt.figure(figsize=(WIDTH, H))
L2 = np.log2
k = pd.read_csv(f"{SD}/figS6_keto.csv")
Q = ["Q1", "Q2", "Q3", "Q4", "Q5"]
med = k[(k["model"] == "measured") & k["bin"].isin(Q)].set_index("bin").loc[Q, "keto_median"]
XT = [f"{i + 1}\n({m:+.3f})".replace("-", "−") for i, m in enumerate(med)]


def series(ax, model, metric, col, mk, dx, label):
    d = k[(k["model"] == model) & (k["metric"] == metric)].set_index("bin").loc[Q]
    for i, q in enumerate(Q):
        r = d.loc[q]
        ci_point(ax, i + dx, L2(r["ratio"]), L2(r["ratio_lo95"]), L2(r["ratio_hi95"]), col, marker=mk, ms=3.4,
                 label=label if i == 0 else None)
    ax.plot(np.arange(5) + dx, L2(d["ratio"].to_numpy(float)), color=col, lw=0.8, zorder=2)


def frame(ax, ylab):
    ax.axhline(0, color=GREY, lw=0.6, zorder=0)
    ax.set_xticks(range(5))
    ax.set_xticklabels(XT, fontsize=6.2)
    ax.set_xlim(-0.5, 4.5)
    log2_ratio_axis(ax, ticks=(0.5, 0.75, 1, 1.5))
    ax.set_ylim(L2(0.5), L2(1.6))
    ax.set_ylabel(ylab)
    ax.set_xlabel("fifth of genes by keto skew of the gene strand\n(median (G + T − A − C)/(A + C + G + T))")


axA = fig.add_axes([0.095, 0.25, 0.36, 0.68])
fig_letter(fig, "A", 0.005, 0.99)
series(axA, "measured", "obs r10", "black", "D", 0.0, "measured, ≥10 reads")
series(axA, "Puffin-D", "top0.0001", PUFFIN, "s", -0.18, "Puffin-D, strong (top 0.01%)")
series(axA, "ProCapNet", "top0.0001", PROCAP, "o", 0.18, "ProCapNet, strong (top 0.01%)")
frame(axA, "strong start sites,\ngene strand / opposite strand")
axA.legend(loc="lower left", handletextpad=0.2, borderaxespad=0.2)

axB = fig.add_axes([0.615, 0.25, 0.36, 0.68])
fig_letter(fig, "B", 0.50, 0.99)
series(axB, "Puffin-D", "top0.001", PUFFIN, "s", -0.12, "Puffin-D, weak (top 0.1%)")
series(axB, "ProCapNet", "top0.001", PROCAP, "o", 0.12, "ProCapNet, weak (top 0.1%)")
frame(axB, "weak potential start sites,\ngene strand / opposite strand")
axB.legend(loc="lower left", handletextpad=0.2, borderaxespad=0.2)
n = int(k[(k["model"] == "measured") & (k["bin"] == "Q1")]["n_units"].iloc[0])
axB.text(1.0, 1.0, f"about {n:,} protein-coding genes per fifth", transform=axB.transAxes, ha="right", va="bottom",
         fontsize=6.2, color=GENE)
save(fig, OUT)
