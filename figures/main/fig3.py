#!/usr/bin/env python3
"""Figure 3: potential and measured start sites on the two strands of introns (version 3).
Data: source_data/fig3_predicted_sites.csv (step 13c), fig3_by_activity.csv and fig3_strong_sites.csv (step 17),
fig3_replication.csv (steps 28 and 29) and fig3_metagene.csv (step 28g, from the site tables of steps 28b and 29b).

Version 3 (12 Nov 2026): each panel asks one question and uses the form that suits it. Panel A is unchanged; no estimate of version 1 changed.
  B  potential start sites: annotated matrix (2 models x 3 kinds of DNA x 3 strengths); colour, log2 ratio on a
     diverging scale whose two hues are those of the two strands in A (blue: fewer sites on the gene strand)
  C  measured sites in K562 at increasing read thresholds, and all reads summed
  D  the two models by gene activity, summed signal and strong sites, with the measured curve of C as a faint
     reference (one axis per facet)
  E  strong measured sites by set of genes with HCT116 PRO-cap and K562 GRO-cap, and a column of values
  F  strong measured sites along active genes in 1-kb bins from the 5' end: mirrored densities (K562 PRO-cap) and
     the strand ratio with its 95% CI for the three data sets (step 28g)"""
import sys

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

from figstyle import (ANTI, GENE, GREY, HCT116, K562, LIGHT, MM, PROCAP, PUFFIN, SENSE, WIDTH, arrow, box, ci_point,
                      log2_ratio_axis, plt, save)

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Figure3"
HMM = 212
fig = plt.figure(figsize=(WIDTH, HMM * MM))
L2 = np.log2
SHADE = "#F3F3F3"
INK2 = "#333333"


def ax_mm(x, top, w, h):
    """Axes from millimetres: left edge, distance of the top edge from the top of the figure, width, height."""
    return fig.add_axes([x / 180, 1 - (top + h) / HMM, w / 180, h / HMM])


def letter(s, x, top):
    fig.text(x / 180, 1 - top / HMM, s, fontsize=10, fontweight="bold", va="top", ha="left")


def subtitle(ax, s, bold=False):
    ax.set_title(s, fontsize=6.5, loc="left", pad=3, color="black" if bold else INK2,
                 fontweight="bold" if bold else "normal")


# ------------------------------------------------------------------ A: what each model reads, to scale
axA = ax_mm(5.4, 1.0, 171, 26.25)
letter("A", 0.9, 1.0)
axA.set_xlim(-2, 162)
axA.set_ylim(-5.5, 6.2)
axA.axis("off")
y0 = 0.0
axA.plot([0, 160], [y0, y0], color=LIGHT, lw=0.8, zorder=0)
# host gene (+), 32-135 kb
axA.plot([32, 135], [y0, y0], color=GENE, lw=1.0, zorder=1)
for a_, b_ in ((32, 33.5), (52, 53.2), (78, 79.2), (104, 105.2), (132, 135)):
    box(axA, a_, y0 - 0.55, b_ - a_, 1.1, fc=GENE, z=2)
arrow(axA, 32, y0 + 0.55, 32, y0 + 1.7, GENE, lw=0.9, style="-")
arrow(axA, 32, y0 + 1.7, 36, y0 + 1.7, GENE, lw=0.9, head=4)
axA.text(36.8, y0 + 1.7, "host gene", va="center", fontsize=6.5, color=GENE)
# neighbouring gene (-), 4-27 kb, start at 27 kb
axA.plot([4, 27], [y0, y0], color=GENE, lw=1.0, zorder=1)
for a_, b_ in ((4, 6), (14, 15), (25.8, 27)):
    box(axA, a_, y0 - 0.55, b_ - a_, 1.1, fc=GENE, z=2)
arrow(axA, 27, y0 - 0.55, 27, y0 - 1.7, GENE, lw=0.9, style="-")
arrow(axA, 27, y0 - 1.7, 23, y0 - 1.7, GENE, lw=0.9, head=4)
axA.text(22.2, y0 - 1.7, "neighbouring gene", va="center", ha="right", fontsize=6.5, color=GENE)
# trace of transcription: sparse reads along each transcribed gene, on its own strand
rng = np.random.default_rng(3)
for x_ in np.sort(rng.uniform(34, 134, 46)):
    axA.plot([x_, x_], [y0 + 0.65, y0 + 1.05], color=SENSE, lw=0.6)
for x_ in np.sort(rng.uniform(5, 26, 12)):
    axA.plot([x_, x_], [y0 - 0.65, y0 - 1.05], color=ANTI, lw=0.6)
# scored position and the two input windows
xs = 75.0
axA.plot([xs, xs], [y0 - 0.9, y0 + 0.9], color="black", lw=1.0, zorder=4)
box(axA, xs - 50, 3.0, 100, 1.0, fc=PUFFIN, alpha=0.35, z=1)
axA.text(xs, 4.35, "Puffin-D sees 100 kb of DNA: here, both promoters", ha="center", va="bottom",
         fontsize=6.5, color=GENE)
box(axA, xs - 1.05, -4.2, 2.1, 1.0, fc=PROCAP, alpha=0.8, z=1)
axA.text(xs + 2.2, -3.7, "ProCapNet sees 2.1 kb of DNA", ha="left", va="center", fontsize=6.5, color=GENE)
axA.text(xs, -1.2, "scored position", ha="center", va="top", fontsize=6.0, color=GENE)
axA.plot([140, 150], [-4.6, -4.6], color=GENE, lw=1.0)
axA.text(145, -4.4, "10 kb", ha="center", va="bottom", fontsize=6.0, color=GENE)


# ------------------------------------------------------------------ B: potential start sites, annotated matrix
ps = pd.read_csv(f"{SD}/fig3_predicted_sites.csv")
ps = ps[ps["seq"] == "all"]
letter("B", 0.9, 36)
axB = ax_mm(30, 45, 39, 52)
MET = [("top0.001", "weak\n(top 0.1%)"), ("top0.0003", "medium\n(top 0.03%)"), ("top0.0001", "strong\n(top 0.01%)")]
KINDS = [("intron_protein_coding", "real", "introns,\noriginal DNA"), ("intron_protein_coding", "shuf", "introns,\nDNA shuffled"),
         ("intergenic_plus_vs_minus", "real", "far from genes,\nplus / minus")]
cmap = LinearSegmentedColormap.from_list("strand", ["#16679E", "#7FB9E0", "#F2F2F2", "#F2BE6E", "#A85F00"])
norm = TwoSlopeNorm(vmin=L2(0.6), vcenter=0, vmax=L2(1 / 0.6))
mat = np.full((7, 3), np.nan)
ylab = []
for mi, (model, col, mk) in enumerate((("Puffin-D", PUFFIN, "s"), ("ProCapNet", PROCAP, "o"))):
    for ki, (st_, pre, lab) in enumerate(KINDS):
        k = mi * 4 + ki
        for c, (m, _) in enumerate(MET):
            r = ps[(ps["model"] == model) & (ps["metric"] == m) & (ps["set"] == st_)].iloc[0]
            v, lo, hi = r[f"{pre}_ratio"], r[f"{pre}_ratio_lo95"], r[f"{pre}_ratio_hi95"]
            mat[k, c] = L2(v)
            axB.text(c, k, f"{v:.2f}\n({lo:.2f}–{hi:.2f})", ha="center", va="center", fontsize=5.5, linespacing=1.1,
                     color="white" if abs(L2(v)) > 0.33 else "#1A1A1A")
    hk = -0.75 if mi == 0 else 3.25
    axB.plot([-0.62], [hk], marker=mk, ms=3.6, color=col, ls="none", clip_on=False)
    axB.text(-0.52, hk, model, ha="left", va="center", fontsize=7, fontweight="bold")
cmap.set_bad("white")
im = axB.imshow(mat, cmap=cmap, norm=norm, aspect="auto")
for k in (0.5, 1.5, 4.5, 5.5):
    axB.axhline(k, color="white", lw=0.8)
for c in (0.5, 1.5):
    axB.axvline(c, color="white", lw=0.8)
axB.set_xticks(range(3))
axB.set_xticklabels([m[1] for m in MET])
axB.set_yticks([0, 1, 2, 4, 5, 6])
axB.set_yticklabels([k[2] for k in KINDS] * 2, fontsize=6, linespacing=1.05)
axB.tick_params(length=0, pad=2)
for s in axB.spines.values():
    s.set_visible(False)
axB.set_xlabel("strength of the potential start site")
axB.set_ylim(6.5, -1.15)
cax = ax_mm(70.5, 56, 2.0, 34)
cb = fig.colorbar(im, cax=cax, ticks=L2([0.6, 0.8, 1, 1.25, 1 / 0.6]))
cb.ax.set_yticklabels(["0.6×", "0.8×", "1×", "1.25×", "1.67×"], fontsize=5.5)
cb.outline.set_visible(False)
cb.ax.tick_params(length=0, pad=1.5)
cb.ax.set_title("gene /\nopposite", fontsize=5.5, pad=5, loc="left")

# ------------------------------------------------------------------ shared drawing of the activity groups
ba = pd.read_csv(f"{SD}/fig3_by_activity.csv")
BINS = [("silent", "0"), ("low", "1–4"), ("A1", "low"), ("A2", "mid"), ("A3", "high")]


def activity_axis(ax, ylim, ticks, xlabels=True, shade_label=True):
    ax.axhline(0, color=GREY, lw=0.6, zorder=1)
    ax.axvspan(1.5, 4.5, color=SHADE, zorder=0, lw=0)
    ax.set_xticks(range(5))
    ax.set_xticklabels([b[1] for b in BINS] if xlabels else [])
    ax.set_xlim(-0.5, 4.5)
    log2_ratio_axis(ax, ticks=ticks)
    ax.set_ylim(L2(ylim[0]), L2(ylim[1]))
    if shade_label:
        ax.text(4.42, L2(ylim[0]) + 0.02 * (L2(ylim[1]) - L2(ylim[0])), "active genes", ha="right", va="bottom",
                fontsize=6.0, color=GREY)


def values(metric):
    r = ba[ba["metric"] == metric].set_index("bin").loc[[b for b, _ in BINS]]
    return L2(r["ratio"].to_numpy()), L2(r["ratio_lo95"].to_numpy()), L2(r["ratio_hi95"].to_numpy())


def series(ax, metric, col, mk, label=None, ms=3.2, filled=True, ls="-", ci=True, alpha=1.0, lw=0.9):
    y, lo, hi = values(metric)
    if ci:
        for i in range(5):
            ax.plot([i, i], [lo[i], hi[i]], color=col, lw=0.8, zorder=3, alpha=alpha, solid_capstyle="butt")
    ax.plot(range(5), y, color=col, lw=lw, ls=ls, zorder=3, alpha=alpha)
    ax.plot(range(5), y, marker=mk, ms=ms, mfc=col if filled else "white", mec=col, mew=0.8, ls="none", zorder=4,
            alpha=alpha)
    return Line2D([], [], color=col, lw=lw, ls=ls, marker=mk, ms=ms, mfc=col if filled else "white", mec=col,
                  mew=0.8, alpha=alpha, label=label)


# ------------------------------------------------------------------ C: measured sites at increasing read thresholds
letter("C", 80, 36)
axC = ax_mm(92, 45, 42, 52)
THR = (("obs sum", "all reads, summed", "#6B6B6B", False, (0, (2.2, 1.4))), ("obs r1", "sites with ≥1 read", "#ABABAB", True, "-"),
       ("obs r3", "sites with ≥3 reads", "#5E5E5E", True, "-"), ("obs r10", "sites with ≥10 reads", "black", True, "-"))
hC = [series(axC, m, c, "D", lab, filled=f, ls=ls) for m, lab, c, f, ls in THR]
activity_axis(axC, (0.5, 18), (0.5, 1, 2, 4, 8, 16))
axC.set_ylabel("measured PRO-cap, K562,\ngene strand / opposite strand")
axC.set_xlabel("activity of the gene\n(reads at its starts)")
axC.legend(handles=hC, loc="upper left", handlelength=2.2, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25,
           fontsize=6.0)
subtitle(axC, "the stricter the threshold, the smaller the excess")

# ------------------------------------------------------------------ D: the models next to the measured data
letter("D", 136, 36)
axD1 = ax_mm(150, 45, 28, 23)
axD2 = ax_mm(150, 76, 28, 21)
for ax, (meas, pmet, cmet), ylim, ticks, ttl in (
        (axD1, ("obs sum", "Puffin-D sum", "ProCapNet sum"), (0.5, 26), (0.5, 1, 2, 4, 8, 16), "summed signal"),
        (axD2, ("obs r10", "Puffin-D strong", "ProCapNet strong"), (0.5, 1.6), (0.5, 0.75, 1, 1.5), "strong sites")):
    hm = series(ax, meas, "black", "D", "measured (from C)", ms=2.6, ci=False, alpha=0.3, lw=0.8)
    hp = series(ax, pmet, PUFFIN, "s", "Puffin-D", ms=3.0)
    hc = series(ax, cmet, PROCAP, "o", "ProCapNet", ms=3.0)
    activity_axis(ax, ylim, ticks, xlabels=ax is axD2, shade_label=False)
    ax.set_ylabel("gene / opposite\nstrand")
    subtitle(ax, ttl, bold=True)
axD2.set_xlabel("activity of the gene")
axD1.legend(handles=[hp, hc, hm], loc="upper left", handlelength=1.6, handletextpad=0.3, borderaxespad=0.1,
            labelspacing=0.12, fontsize=5.8)

# ------------------------------------------------------------------ E: strong measured sites by set of genes, three data sets
st = pd.read_csv(f"{SD}/fig3_strong_sites.csv")
rp = pd.read_csv(f"{SD}/fig3_replication.csv")
DS = (("K562 PRO-cap", K562, True), ("HCT116 PRO-cap", HCT116, True), ("K562 GRO-cap", K562, False))
letter("E", 0.9, 113)
axE = ax_mm(30, 122, 28, 66)
ROWS = [("pc_all", None, "protein-coding genes"), ("pc_active", "active", "active"), ("pc_silent", None, "silent"),
        ("lncRNA", None, "lncRNA genes"), ("intergenic", "intergenic", "intergenic,\nplus / minus")]
tE = blended_transform_factory(axE.transAxes, axE.transData)
for i, (b, rsub, lab) in enumerate(ROWS):
    yy = len(ROWS) - 1 - i
    if i % 2 == 0:
        axE.axhspan(yy - 0.5, yy + 0.5, xmin=-1.05, xmax=1.95, color="#F2F2F2", lw=0, zorder=0, clip_on=False)
    r = st[st["bin"] == b].iloc[0]
    unit = "1-Mb blocks" if b == "intergenic" else "genes"
    pts = [(DS[0], (r["ratio"], r["ratio_lo95"], r["ratio_hi95"]))]
    if rsub:
        for ds in DS[1:]:
            q = rp[(rp["data_set"] == ds[0]) & (rp["subset"] == rsub)].iloc[0]
            pts.append((ds, (q["ratio"], q["ratio_lo95"], q["ratio_hi95"])))
    offs = [0.0] if len(pts) == 1 else [0.25, 0.0, -0.25]
    for ((ds, col, filled), (v, lo, hi)), dy in zip(pts, offs):
        ci_point(axE, L2(v), yy + dy, L2(lo), L2(hi), col, marker="D", ms=3.0, horizontal=True, filled=filled)
        axE.text(1.06, yy + dy, f"{v:.2f} ({lo:.2f}–{hi:.2f})", transform=tE, ha="left", va="center", fontsize=5.6,
                 color=INK2)
    axE.text(-0.05, yy, f"{lab}\n({int(r['n_units']):,} {unit if b != 'intergenic' else 'blocks'})", transform=tE,
             ha="right", va="center", fontsize=6.0, linespacing=1.05)
axE.text(1.06, len(ROWS) - 0.45, "ratio (95% CI)", transform=tE, ha="left", va="bottom", fontsize=5.6, color=INK2)
axE.axvline(0, color=GREY, lw=0.6, zorder=1)
axE.set_yticks([])
axE.spines["left"].set_visible(False)
axE.set_ylim(-0.5, len(ROWS) - 0.5)
log2_ratio_axis(axE, ticks=(0.5, 0.75, 1, 1.5), axis="x")
axE.set_xlim(L2(0.5), L2(1.6))
axE.set_xlabel("strong measured sites,\ngene strand / opposite strand")
# one key for E and F, at the bottom of the figure
hE = [Line2D([], [], marker="D", ms=3.2, color=c, mfc=c if f else "white", mew=0.8, lw=1.1,
             ls="-" if f else (0, (3, 1.6)), label=f"{d}, strong sites: ≥{t} reads")
      for (d, c, f), t in zip(DS, (10, 24, 14))]
fig.legend(handles=hE, loc="lower center", bbox_to_anchor=(0.5, 2.0 / HMM), ncol=3, handlelength=2.4,
           handletextpad=0.4, columnspacing=2.0, borderaxespad=0, fontsize=6.0, frameon=False)

# ------------------------------------------------------------------ F: along active genes, 1-kb bins (step 28g)
mg = pd.read_csv(f"{SD}/fig3_metagene.csv")
letter("F", 94, 113)
axF1 = ax_mm(110, 122, 68, 33)
axF2 = ax_mm(110, 158, 68, 30)
k = mg[(mg["data_set"] == "K562 PRO-cap") & mg["bin_hi_kb"].notna()]
xc = (k["bin_lo_kb"] + k["bin_hi_kb"]).to_numpy() / 2
axF1.bar(xc, k["sense_per1000"], width=0.86, color=SENSE, lw=0, zorder=2)
axF1.bar(xc, -k["anti_per1000"], width=0.86, color=ANTI, lw=0, zorder=2)
axF1.axhline(0, color="black", lw=0.6, zorder=3)
top = float(np.ceil(max(k["sense_per1000"].max(), k["anti_per1000"].max()) / 5) * 5)
axF1.set_ylim(-top, top)
step = 5 if top <= 15 else 10
tk = np.arange(-(top // step) * step, top + 0.1, step)
axF1.set_yticks(tk)
axF1.set_yticklabels([f"{abs(t):g}" for t in tk])
axF1.set_ylabel("strong sites per 1,000\ngenes reaching the bin")
axF1.set_xlim(0, 30)
axF1.set_xticklabels([])
axF1.axvspan(0, 5, color="black", alpha=0.05, lw=0, zorder=0)
axF1.text(2.5, top * 0.93, "first\n5 kb", ha="center", va="top", fontsize=6.0, color=INK2, linespacing=1.05)
axF1.legend(handles=[Patch(color=SENSE, label="gene strand (sense)"), Patch(color=ANTI, label="opposite strand (antisense)")],
            loc="upper right", handlelength=1.0, handleheight=0.8, handletextpad=0.4, borderaxespad=0.1, fontsize=6.0)
subtitle(axF1, "K562 PRO-cap, strong sites (≥10 reads), active genes")
COL = "ratio3" if "ratio3" in mg.columns else "ratio"
for (ds, col, filled) in DS:
    d = mg[(mg["data_set"] == ds) & mg["bin_hi_kb"].notna()]
    x = (d["bin_lo_kb"] + d["bin_hi_kb"]).to_numpy() / 2
    lo_c, hi_c = ("lo3", "hi3") if COL == "ratio3" else ("lo", "hi")
    axF2.fill_between(x, L2(d[lo_c]), L2(d[hi_c]), color=col, alpha=0.12 if filled else 0.0, lw=0, zorder=1)
    if not filled:
        axF2.plot(x, L2(d[lo_c]), color=col, lw=0.5, ls=(0, (1.5, 1.5)), zorder=2)
        axF2.plot(x, L2(d[hi_c]), color=col, lw=0.5, ls=(0, (1.5, 1.5)), zorder=2)
    axF2.plot(x, L2(d[COL]), color=col if ds != "K562 GRO-cap" else col, lw=1.1, ls="-" if filled else (0, (3, 1.6)),
              zorder=3, label=ds)
axF2.axhline(0, color=GREY, lw=0.6, zorder=1)
axF2.axvspan(0, 5, color="black", alpha=0.05, lw=0, zorder=0)
log2_ratio_axis(axF2, ticks=(0.25, 0.5, 1, 2))
axF2.set_ylim(L2(0.25), L2(2.7))
axF2.set_xlim(0, 30)
axF2.set_xlabel("distance from the annotated 5′ end of the gene (kb)")
axF2.set_ylabel("gene / opposite\nstrand (95% CI)")
save(fig, OUT)
