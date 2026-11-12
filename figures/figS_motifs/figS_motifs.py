#!/usr/bin/env python3
"""Supplemental figure: sequence at the two starts of the depth-matched pairs and the pre-specified core-promoter
elements (step 27a). Drawn from the step-27a result tables; no other input.

Inputs (checked by md5; from 27a_results.tar.gz of step 27a, script 27a_start_motifs.py md5 cbe348734af56180b9bd2a48db278095)
  27a_pfm.tsv    base counts per group and position (-40..+40; 0 = first transcribed base)   bde79d95af096729398db01dac9b356e
  27a_tests.tsv  per feature: means per group, differences, bootstrap CIs and p, Holm p      09407c3544bb97cfb84e2c13b687cd11
Groups: IS / IA = sense / antisense start of intronic pairs; GP / GM = plus / minus start of their depth-matched
intergenic partners (drawn together, as the label 'sense' is arbitrary between genes).
Panel A: sequence logos, information content in bits (2 + sum p log2 p, no small-sample correction).
Panel B: for each element, (sense - antisense, intronic) - (plus - minus, intergenic) of the share of starts carrying
the element, observed minus composition-expected (20 dinucleotide shuffles of +/-60 nt per start), in percentage
points, with 95% bootstrap CI over host genes; filled symbols: the five pre-specified elements (Holm-corrected p
over these five within each cell line printed at right); open: TATA box with one mismatch (sensitivity).

Usage: python3 figS_motifs.py [data_dir] [output_prefix]
  defaults: data_dir = data, output_prefix = figures/Szczesniak_FigS_motifs
Writes <prefix>.png (600 dpi), .pdf, .svg and the Source Data <prefix>_logos.csv and <prefix>_elements.csv.
Needs figstyle.py and the fonts/ folder (Liberation Sans) in the same folder as this script.
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd
from matplotlib import font_manager as fm
from matplotlib import transforms
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath

from figstyle import GREY, HCT116, K562, LIGHT, MM, ci_point, fig_letter, plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else "data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Szczesniak_FigS_motifs"
MD5 = {"27a_pfm.tsv": "bde79d95af096729398db01dac9b356e", "27a_tests.tsv": "09407c3544bb97cfb84e2c13b687cd11"}


def md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


for f_, m_ in MD5.items():
    if md5(os.path.join(DATA, f_)) != m_:
        sys.exit(f"STOP: md5 of {f_} differs from the step-27a result")
fdir = os.path.join(HERE, "fonts")
if os.path.isdir(fdir):
    for f_ in sorted(os.listdir(fdir)):
        if f_.lower().endswith(".ttf"):
            fm.fontManager.addfont(os.path.join(fdir, f_))
try:
    fm.findfont(FontProperties(family="Liberation Sans"), fallback_to_default=False)
except ValueError:
    sys.exit("STOP: Liberation Sans not found; keep the fonts/ folder next to the script")
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "svg.fonttype": "none", "svg.hashsalt": "PURGE-FigS-motifs"})

PFM = pd.read_csv(os.path.join(DATA, "27a_pfm.tsv"), sep="\t")
TST = pd.read_csv(os.path.join(DATA, "27a_tests.tsv"), sep="\t")
for col in ("line", "group", "pos", "A", "C", "G", "T", "n"):
    if col not in PFM.columns:
        sys.exit(f"STOP: column {col} missing in 27a_pfm.tsv")
LINES = [("K562", K562), ("HCT116", HCT116)]
ROWS = [("IS", "intronic pairs, sense start"), ("IA", "intronic pairs, antisense start"),
        ("GP+GM", "intergenic partners, both starts")]
BASE_COL = {"A": "#1A9850", "C": "#2166AC", "G": "#F0A30A", "T": "#C2272D"}
LOGO_FONT = FontProperties(family="Liberation Sans", weight="bold")
WINDOWS = [("TATA", -36, -17), ("Inr", -6, 5), ("DPE", 24, 36)]   # positions a pre-specified element may cover


def counts(line, groups):
    d = PFM[(PFM["line"] == line) & (PFM["group"].isin(groups))]
    g = d.groupby("pos")[["A", "C", "G", "T", "n"]].sum().sort_index()
    return g


def letter(ax, ch, x, y0, h, w=0.86):
    if h <= 0:
        return
    tp = TextPath((0, 0), ch, size=1, prop=LOGO_FONT)
    bb = tp.get_extents()
    tr = transforms.Affine2D().translate(-bb.x0, -bb.y0).scale(w / bb.width, h / bb.height).translate(x - w / 2, y0)
    ax.add_patch(PathPatch(tr.transform_path(tp), fc=BASE_COL[ch], ec="none", lw=0, zorder=3))


logo_rows = []
fig = plt.figure(figsize=(180 * MM, 140 * MM))
L0, L1, COLW, GAPX = 0.065, 0.985, 0.43, 0.05
for j, (line, color) in enumerate(LINES):
    for i, (grp, title) in enumerate(ROWS):
        ax = fig.add_axes([L0 + j * (COLW + GAPX), 0.835 - i * 0.16, COLW, 0.118])
        g = counts(line, grp.split("+"))
        P = g[["A", "C", "G", "T"]].to_numpy(float) / g["n"].to_numpy(float)[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            ic = 2 + np.nansum(np.where(P > 0, P * np.log2(P), 0), axis=1)
        for pos, p, info in zip(g.index, P, ic):
            y = 0.0
            for b, pb in sorted(zip("ACGT", p), key=lambda t: t[1]):
                letter(ax, b, pos, y, pb * info)
                y += pb * info
            logo_rows.append({"line": line, "group": grp, "pos": int(pos), **{b: int(g.loc[pos, b]) for b in "ACGT"},
                              "n": int(g.loc[pos, "n"]), "information_bits": round(float(info), 5)})
        for nm, a_, b_ in WINDOWS:
            ax.axvspan(a_ - 0.5, b_ + 0.5, color=LIGHT, alpha=0.45, lw=0, zorder=0)
            if i == 0:
                ax.text((a_ + b_) / 2, 1.04, nm, ha="center", va="bottom", fontsize=6.5, color=GREY,
                        transform=ax.get_xaxis_transform())
        ax.axvline(0, color=GREY, lw=0.4, ls=":", zorder=1)
        ax.set_xlim(-40.6, 40.6)
        ax.set_ylim(0, 1.0)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_xticks(range(-40, 41, 10))
        n_ = int(g["n"].iloc[0])
        ax.set_title(f"{line}: {title} (n = {n_:,})", fontsize=7.5, pad=10 if i == 0 else 3, loc="left", color=color)
        if j == 0:
            ax.set_ylabel("bits")
        else:
            ax.set_yticklabels([])
        if i == len(ROWS) - 1:
            ax.set_xlabel("position relative to the start (0 = first transcribed base)")
        else:
            ax.set_xticklabels([])
fig_letter(fig, "A", 0.008, 0.985)

# ---- panel B: forest of the intronic-minus-intergenic difference
ELEM = [("TATA", "TATA box", True), ("Inr", "Initiator (Inr)", True), ("DPE", "Downstream promoter element", True),
        ("YY1f", "YY1 motif, forward", True), ("YY1r", "YY1 motif, reverse", True),
        ("TATA1", "TATA box, one mismatch", False)]
axb = fig.add_axes([0.30, 0.095, 0.50, 0.29])
el_rows = []
for k, (feat, name, prim) in enumerate(ELEM):
    y0 = len(ELEM) - 1 - k
    for j, (line, color) in enumerate(LINES):
        r = TST[(TST["line"] == line) & (TST["feature"] == feat) & (TST["measure"] == "adj")]
        if len(r) != 1:
            sys.exit(f"STOP: {line} {feat} adj: {len(r)} rows in 27a_tests.tsv")
        r = r.iloc[0]
        y = y0 + (0.16 if j == 0 else -0.16)
        ci_point(axb, 100 * r["Diff"], y, 100 * r["Diff_lo95"], 100 * r["Diff_hi95"], color, horizontal=True, filled=prim,
                 ms=3.4)
        ph = r["Diff_p_holm"]
        txt = (f"{line}: Holm P = {ph:.2g}" if prim else f"{line}: P = {r['Diff_p']:.2g}")
        axb.text(1.02, y, txt, transform=axb.get_yaxis_transform(), ha="left", va="center", fontsize=6.5, color=color)
        el_rows.append({"line": line, "feature": feat, "name": name, "pre_specified": prim, "n_pairs": int(r["n_pairs"]),
                        "n_host_genes": int(r["n_host_genes"]), "diff_pp": 100 * r["Diff"], "lo95_pp": 100 * r["Diff_lo95"],
                        "hi95_pp": 100 * r["Diff_hi95"], "p": r["Diff_p"], "p_holm": ph if prim else np.nan})
axb.axvline(0, color=GREY, lw=0.6, zorder=1)
axb.axhline(0.5, color=LIGHT, lw=0.6, zorder=0)
axb.set_yticks(range(len(ELEM)))
axb.set_yticklabels([nm for _, nm, _ in ELEM][::-1])
axb.set_ylim(-0.6, len(ELEM) - 0.4)
axb.set_xlabel("(sense − antisense, intronic) − (plus − minus, intergenic),\npercentage points of starts, observed − expected")
lim = max(abs(v) for row in el_rows for v in (row["lo95_pp"], row["hi95_pp"])) * 1.1
axb.set_xlim(-lim, lim)
fig_letter(fig, "B", 0.008, 0.405)

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf", metadata={"CreationDate": None, "ModDate": None})
fig.savefig(OUT + ".svg", metadata={"Date": None})
plt.close(fig)
pd.DataFrame(logo_rows).to_csv(OUT + "_logos.csv", index=False)
pd.DataFrame(el_rows).to_csv(OUT + "_elements.csv", index=False, float_format="%.6g")
print(f"wrote {OUT}.png/.pdf/.svg, {OUT}_logos.csv, {OUT}_elements.csv")
