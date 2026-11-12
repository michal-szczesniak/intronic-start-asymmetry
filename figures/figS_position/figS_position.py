#!/usr/bin/env python3
"""Supplemental figure: where along the gene the strand difference lies (steps 28, 28b, 28e, 28f and 29).
Drawn from the result tables of these steps; no other input.

Inputs (checked by md5)
  28_K562_strata.tsv, 28_HCT116_strata.tsv    step 28 (28_k36_position.py c020e287552668cec7a74a786f2aa5b4)
  28b_K562_strata.tsv, 28b_HCT116_strata.tsv  step 28b (28b_position.py 08504653e2a0ea5742c10d0bd714461a)
  28e_cells.tsv                               step 28e (28e_twobytwo.py d053a8fead86cacf3de182ed2af6fe07)
  28f_cells.tsv                               step 28f (28f_activity.py 3d020b9596c91944684deb57af56a65b)
  29_K562_GROcap_strata.tsv, 29_GROcap_28b_strata.tsv   step 29b: the step-28 and step-28b scripts on K562 GRO-cap
                                              (GSM1480321, lifted over from hg19; 29b_liftover.py 68eb7bd5b143bd3f60686c1f3300b81f)
Version b (11 Nov 2026): K562 GRO-cap added to panels A and B (open blue diamonds; pre-specified in step 29).
Panels
  A  Deficit of strong measured sites, 100 x (1 - sense / antisense), in fifths of active protein-coding genes
     (5' to 3'), K562 and HCT116 PRO-cap and K562 GRO-cap (steps 28 and 29; pre-specified contrast 3' fifth minus
     5' fifth printed).
  B  The same deficit by distance of the site from the annotated 5' end of the host gene (step 28b, exploratory, for
     PRO-cap; step 29, pre-specified, for GRO-cap).
  C  Depth-matched bias of the start-site pairs, 100 x (1 - 2^d), by distance of the pair from the 5' end of the
     host gene (step 28b; exploratory).
  D  Both cell lines pooled (step 28e; exploratory): strong sites and pairs within 5 kb of the 5' end or farther,
     within 2 kb of a CpG island (filled) or farther (open).
  E  Both cell lines pooled (step 28f; exploratory): excess within 5 kb (within minus beyond, percentage points) in
     thirds of host activity.
Bars: 95% bootstrap CIs over genes (host genes for pairs); intervals beyond the axis end in an arrowhead.

Usage: python3 figS_position.py [data_dir] [output_prefix]
  defaults: data_dir = data, output_prefix = figures/Szczesniak_FigS_position
Writes <prefix>.png (600 dpi), .pdf, .svg and the Source Data <prefix>_source_data.csv.
Needs figstyle.py and the fonts/ folder (Liberation Sans) in the same folder as this script.
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D

from figstyle import GREY, HCT116, K562, LIGHT, MM, OBS, fig_letter, plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Szczesniak_FigS_position"
MD5 = {"28_K562_strata.tsv": "e3b5e9eeb25d487a1d0b346eef5ca664",
       "28_HCT116_strata.tsv": "82f53f3dbcccd96709260750ccffaa2b",
       "28b_K562_strata.tsv": "8ecf4540c3f3bd108ae910a2ea006bb7",
       "28b_HCT116_strata.tsv": "fbaed9ebac905eec7189d67a27ed2a11",
       "28e_cells.tsv": "3ae8b9c9c48a6b44ee253465ef2bd68f",
       "28f_cells.tsv": "430a0d9ff54b9918fe5af451292fb925",
       "29_K562_GROcap_strata.tsv": "15a7e68bbc6c418b55ef6188bcc13dc6",
       "29_GROcap_28b_strata.tsv": "7978f868ca2a43260e1aa7537df5430a"}


def md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


for f_, m_ in MD5.items():
    if md5(os.path.join(DATA, f_)) != m_:
        sys.exit(f"STOP: md5 of {f_} differs from the step result")
fdir = os.path.join(HERE, "fonts")
if os.path.isdir(fdir):
    for f_ in sorted(os.listdir(fdir)):
        if f_.lower().endswith(".ttf"):
            fm.fontManager.addfont(os.path.join(fdir, f_))
try:
    fm.findfont(FontProperties(family="Liberation Sans"), fallback_to_default=False)
except ValueError:
    sys.exit("STOP: Liberation Sans not found; keep the fonts/ folder next to the script")
plt.rcParams.update({"svg.hashsalt": "PURGE-FigS-position"})

LINES = [("K562", K562), ("HCT116", HCT116)]
DX = {"K562": -0.13, "HCT116": 0.13}
S28 = {ln: pd.read_csv(os.path.join(DATA, f"28_{ln}_strata.tsv"), sep="\t") for ln, _ in LINES}
S28B = {ln: pd.read_csv(os.path.join(DATA, f"28b_{ln}_strata.tsv"), sep="\t") for ln, _ in LINES}
C28E = pd.read_csv(os.path.join(DATA, "28e_cells.tsv"), sep="\t")
C28F = pd.read_csv(os.path.join(DATA, "28f_cells.tsv"), sep="\t")
G28 = pd.read_csv(os.path.join(DATA, "29_K562_GROcap_strata.tsv"), sep="\t")
G28B = pd.read_csv(os.path.join(DATA, "29_GROcap_28b_strata.tsv"), sep="\t")
# strong sites: three data sets (label, colour, filled, step-28 table, step-28b table, x offset, source names)
SITES = [("K562 PRO-cap", K562, True, S28["K562"], S28B["K562"], -0.2, "28_K562_strata.tsv", "28b_K562_strata.tsv"),
         ("K562 GRO-cap", K562, False, G28, G28B, 0.0, "29_K562_GROcap_strata.tsv", "29_GROcap_28b_strata.tsv"),
         ("HCT116 PRO-cap", HCT116, True, S28["HCT116"], S28B["HCT116"], 0.2, "28_HCT116_strata.tsv",
          "28b_HCT116_strata.tsv")]
rows = []


def one(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    if m.sum() != 1:
        sys.exit(f"STOP: {kw}: {int(m.sum())} rows instead of 1")
    return df[m].iloc[0]


def point(ax, x, est, lo, hi, color, ylim, filled=True, ms=3.6):
    """Diamond (measured data) with a 95% CI bar, clipped to the axis with an arrowhead."""
    y0, y1 = ylim
    lo_c, hi_c = max(lo, y0), min(hi, y1)
    ax.plot([x, x], [lo_c, hi_c], color=color, lw=0.9, zorder=3, solid_capstyle="butt")
    for v, lim, mk in ((lo, y0, "v"), (hi, y1, "^")):
        if (v < lim) if mk == "v" else (v > lim):
            ax.plot([x], [lim], marker=mk, ms=3.0, color=color, mew=0, zorder=4, clip_on=False)
    ax.plot([x], [est], marker="D", ms=ms, mfc=color if filled else "white", mec=color, mew=0.8, ls="none", zorder=5)


def frame(ax, ylim, ylabel=None):
    ax.axhline(0, color=GREY, lw=0.6, zorder=1)
    ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel)


def title(ax, s):
    ax.set_title(s, fontsize=7, loc="left", pad=4)


fig = plt.figure(figsize=(180 * MM, 128 * MM))
TOP, H1, H2, BOT = 0.905, 0.29, 0.29, 0.115
W3, X3 = 0.25, [0.075, 0.405, 0.735]
YL_SITE = (-60, 75)


def mn(v):
    """Number with a typographic minus sign."""
    return f"{v:.1f}".replace("-", "\u2212")
YL_PAIR = (-100, 100)

# ---- A: fifths of the gene (steps 28 and 29, pre-specified)
axa = fig.add_axes([X3[0], TOP - H1, W3, H1])
for lab, col, filled, t28, _, dx, src, _ in SITES:
    for i, p in enumerate(["P1", "P2", "P3", "P4", "P5"]):
        r = one(t28, subset="active", metric="strong", stratum=p)
        point(axa, i + dx, r["deficit"], r["lo"], r["hi"], col, YL_SITE, filled=filled, ms=3.2)
        rows.append({"panel": "A", "scope": lab, "what": "strong sites, active protein-coding genes", "group": p,
                     "n_genes": int(r["n_units"]), "sense": int(r["sense"]), "antisense": int(r["anti"]),
                     "estimate": r["deficit"], "lo95": r["lo"], "hi95": r["hi"], "p_boot": np.nan, "source": src})
    c = one(t28, subset="active", metric="strong", stratum="P5-P1")
    rows.append({"panel": "A", "scope": lab, "what": "contrast 3' fifth minus 5' fifth (pre-specified)",
                 "group": "P5-P1", "n_genes": int(c["n_units"]), "sense": np.nan, "antisense": np.nan,
                 "estimate": c["deficit"], "lo95": c["lo"], "hi95": c["hi"], "p_boot": c["p_boot"], "source": src})
frame(axa, YL_SITE, "deficit of strong sense sites, %")
axa.set_xticks(range(5))
axa.set_xticklabels(["1st\n(5′)", "2nd", "3rd", "4th", "5th\n(3′)"])
axa.set_xlabel("fifth of the gene")
axa.set_xlim(-0.5, 4.5)
title(axa, "Strong sites, fifths of the gene")
axa.text(0.02, 0.20, "3′ fifth − 5′ fifth (pre-specified):", transform=axa.transAxes, ha="left", va="bottom",
         fontsize=5.5, color=GREY)
for k, (lab, col, filled, t28, _, _, _, _) in enumerate(SITES):
    c = one(t28, subset="active", metric="strong", stratum="P5-P1")
    pv = f"P = {c['p_boot']:.2g}" if c["p_boot"] >= 0.001 else "P < 0.001"
    axa.text(0.02, 0.14 - 0.062 * k, f"{lab} {mn(c['deficit'])} ({mn(c['lo'])} to {mn(c['hi'])}), {pv}",
             transform=axa.transAxes, ha="left", va="bottom", fontsize=5.5, color=col,
             fontstyle="normal" if filled else "italic")

# ---- B and C: distance from the 5' end (step 28b)
BINS = ["1-5", "5-10", "10-20", "20-40", "40-80", ">=80"]
BLAB = ["1–5", "5–10", "10–20", "20–40", "40–80", "≥80"]
for ax_x, an, ylim, ylab, ttl, pnl in ((X3[1], "sites_strong_bin5", YL_SITE, "deficit of strong sense sites, %",
                                        "Strong sites, distance from the 5′ end", "B"),
                                       (X3[2], "pairs_dist_bin", YL_PAIR, "bias of the start-site pairs, %",
                                        "Start-site pairs, distance from the 5′ end", "C")):
    ax = fig.add_axes([ax_x, TOP - H1, W3, H1])
    if pnl == "B":
        series = [(lab, col, filled, t28b, dx, src) for lab, col, filled, _, t28b, dx, _, src in SITES]
    else:
        series = [(ln, col, True, S28B[ln], DX[ln], f"28b_{ln}_strata.tsv") for ln, col in LINES]
    for lab, col, filled, t28b, dx, src in series:
        for i, b in enumerate(BINS):
            r = one(t28b, analysis=an, stratum=b)
            point(ax, i + dx, r["estimate"], r["lo"], r["hi"], col, ylim, filled=filled, ms=3.2 if pnl == "B" else 3.6)
            rows.append({"panel": pnl, "scope": lab, "what": an, "group": f"{b} kb",
                         "n_genes": int(r["units"]), "sense": r["sense"], "antisense": r["anti"],
                         "n_pairs": int(r["n"]) if an.startswith("pairs") else np.nan,
                         "estimate": r["estimate"], "lo95": r["lo"], "hi95": r["hi"], "p_boot": np.nan,
                         "source": src})
        if pnl == "B":
            r = one(t28b, analysis=an, stratum="1-5-rest")
            rows.append({"panel": "B", "scope": lab, "what": "contrast 1-5 kb minus the rest (pre-specified for GRO-cap)"
                         if not filled else "contrast 1-5 kb minus the rest", "group": "1-5-rest",
                         "n_genes": int(r["units"]), "sense": np.nan, "antisense": np.nan, "estimate": r["estimate"],
                         "lo95": r["lo"], "hi95": r["hi"], "p_boot": r["p_boot"], "source": src})
    frame(ax, ylim, ylab)
    ax.axvspan(-0.5, 0.5, color=LIGHT, alpha=0.5, lw=0, zorder=0)
    ax.set_xticks(range(len(BINS)))
    ax.set_xticklabels(BLAB)
    ax.set_xlim(-0.5, len(BINS) - 0.5)
    ax.set_xlabel("distance from the 5′ end of the host gene, kb")
    title(ax, ttl)
    if pnl == "B":
        ax.legend(handles=[Line2D([], [], marker="D", ms=3.4, mfc=c_ if f_ else "white", mec=c_, mew=0.8, lw=0,
                                  label=l_) for l_, c_, f_, *_ in SITES],
                  loc="upper right", handletextpad=0.3, borderaxespad=0.2)
    elif pnl == "C":
        ax.legend(handles=[Line2D([], [], marker="D", ms=3.6, color=c_, lw=0, label=f"{l_} PRO-cap")
                           for l_, c_ in LINES], loc="upper right", handletextpad=0.3, borderaxespad=0.2)

# ---- D: pooled 2 x 2 (step 28e)
axd = fig.add_axes([0.075, BOT, 0.50, H2])
CELLS = [("near start, away", "≤ 5 kb from\nthe 5′ end", False), ("near start, island", None, True),
         ("far start, away", "> 5 kb from\nthe 5′ end", False), ("far start, island", None, True)]
GROUPS = [("sites", "strong sites (deficit)", 0.0), ("pairs", "start-site pairs (bias)", 2.6)]
YL_D = (-20, 100)
for what, glab, x0 in GROUPS:
    for k, (cell, _, filled) in enumerate(CELLS):
        r = one(C28E, what=what, scope="pooled", cell=cell)
        x = x0 + (k // 2) * 1.0 + (0.17 if filled else -0.17)
        point(axd, x, r["estimate"], r["lo"], r["hi"], OBS, YL_D, filled=filled)
        rows.append({"panel": "D", "scope": "pooled", "what": what, "group": cell, "n_genes": np.nan,
                     "sense": np.nan, "antisense": np.nan, "n_pairs": int(r["n"]) if what == "pairs" else np.nan,
                     "n_sites": int(r["n"]) if what == "sites" else np.nan, "estimate": r["estimate"],
                     "lo95": r["lo"], "hi95": r["hi"], "p_boot": np.nan, "source": "28e_cells.tsv"})
    sub = C28E[(C28E["what"] == what) & (C28E["scope"] == "pooled")]
    for lab in [c_ for c_ in sub["cell"] if str(c_).startswith(("R1", "R2", "R3", "R4"))]:
        r = one(C28E, what=what, scope="pooled", cell=lab)
        rows.append({"panel": "D", "scope": "pooled", "what": what, "group": lab, "n_genes": np.nan,
                     "sense": np.nan, "antisense": np.nan, "estimate": r["estimate"], "lo95": r["lo"],
                     "hi95": r["hi"], "p_boot": r["p_boot"], "source": "28e_cells.tsv"})
    axd.text(x0 + 0.5, YL_D[1] * 1.02, glab, ha="center", va="bottom", fontsize=6.5)
frame(axd, YL_D, "deficit or bias, %")
ticks, labs = [], []
for _, _, x0 in GROUPS:
    ticks += [x0, x0 + 1.0]
    labs += ["≤ 5 kb", "> 5 kb"]
axd.set_xticks(ticks)
axd.set_xticklabels(labs)
axd.set_xlim(-0.55, 4.15)
axd.axvline(1.8, ymax=0.8, color=LIGHT, lw=0.6, zorder=0)
axd.set_xlabel("distance from the 5′ end of the host gene (both cell lines pooled)")
axd.legend(handles=[Line2D([], [], marker="D", ms=3.6, mfc="white", mec=OBS, mew=0.8, lw=0,
                           label="> 2 kb from any CpG island"),
                    Line2D([], [], marker="D", ms=3.6, mfc=OBS, mec=OBS, lw=0, label="≤ 2 kb from a CpG island")],
           loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.3, columnspacing=1.5,
           borderaxespad=0.2)

# ---- E: excess within 5 kb by activity third (step 28f)
axe = fig.add_axes([0.705, BOT, 0.28, H2])
YL_E = (-20, 80)
SER = [("sites", "strong sites", OBS, -0.13), ("pairs", "start-site pairs", GREY, 0.13)]
for what, lab, col, dx in SER:
    for i, t in enumerate(["T1", "T2", "T3"]):
        r = one(C28F, what=what, subset="all", scope="pooled", third=t, cell="near minus far")
        point(axe, i + dx, r["estimate"], r["lo"], r["hi"], col, YL_E)
        rows.append({"panel": "E", "scope": "pooled", "what": what, "group": t, "n_genes": np.nan, "sense": np.nan,
                     "antisense": np.nan, "estimate": r["estimate"], "lo95": r["lo"], "hi95": r["hi"],
                     "p_boot": r["p_boot"], "source": "28f_cells.tsv"})
    r = one(C28F, what=what, subset="all", scope="pooled", third="T3-T1", cell="interaction")
    rows.append({"panel": "E", "scope": "pooled", "what": what, "group": "T3-T1 (interaction)", "n_genes": np.nan,
                 "sense": np.nan, "antisense": np.nan, "estimate": r["estimate"], "lo95": r["lo"], "hi95": r["hi"],
                 "p_boot": r["p_boot"], "source": "28f_cells.tsv"})
frame(axe, YL_E, "excess within 5 kb, percentage points")
axe.set_xticks(range(3))
axe.set_xticklabels(["least", "middle", "most"])
axe.set_xlim(-0.5, 2.5)
axe.set_xlabel("third of host activity (pooled)")
axe.legend(handles=[Line2D([], [], marker="D", ms=3.6, color=c_, lw=0, label=l_) for _, l_, c_, _ in SER],
           loc="upper left", handletextpad=0.3, borderaxespad=0.2)

for s_, x_, y_ in (("A", 0.008, 0.985), ("B", 0.338, 0.985), ("C", 0.668, 0.985), ("D", 0.008, 0.48),
                   ("E", 0.638, 0.48)):
    fig_letter(fig, s_, x_, y_)

# fit check: every text inside the figure
fig.canvas.draw()
rend = fig.canvas.get_renderer()
for t in fig.findobj(match=lambda o: hasattr(o, "get_window_extent") and o.__class__.__name__ in ("Text", "Annotation")):
    if not t.get_text().strip() or not t.get_visible():
        continue
    bb = t.get_window_extent(rend)
    W, H = fig.bbox.width, fig.bbox.height
    if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
        sys.exit(f"STOP: text outside the figure: {t.get_text()!r}")

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf", metadata={"CreationDate": None, "ModDate": None})
fig.savefig(OUT + ".svg", metadata={"Date": None})
plt.close(fig)
cols = ["panel", "scope", "what", "group", "n_genes", "n_pairs", "n_sites", "sense", "antisense", "estimate", "lo95",
        "hi95", "p_boot", "source"]
pd.DataFrame(rows).reindex(columns=cols).to_csv(OUT + "_source_data.csv", index=False, float_format="%.6g")
print(f"wrote {OUT}.png/.pdf/.svg and {OUT}_source_data.csv ({len(rows)} rows)")
