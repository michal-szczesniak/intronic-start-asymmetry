#!/usr/bin/env python3
"""Figure 2: host activity, the share predicted from sequence and where it lies.
Data: source_data/fig2_host_activity.csv, fig2_host_trend.csv (step 23b), fig2_shares.csv (steps 23b and 15).

Rebuilt on 11 Nov 2026 from an earlier script and the image of the figure in the manuscript (the script of that
image was not kept); the output is identical to that image pixel for pixel."""
import sys

import pandas as pd

from figstyle import (ANTI, GENE, GREY, HCT116, K562, MM, SENSE, WIDTH, arrow, box, ci_point, fig_letter,
                      plt, save)

SD = sys.argv[1] if len(sys.argv) > 1 else "source_data"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/Figure2"
H = 128 * MM
fig = plt.figure(figsize=(WIDTH, H))
LINES = (("K562", K562, -0.13), ("HCT116", HCT116, 0.13))

# ------------------------------------------------------------------ A: by host activity
ha = pd.read_csv(f"{SD}/fig2_host_activity.csv")
tr = pd.read_csv(f"{SD}/fig2_host_trend.csv")
mt = pd.read_csv(f"{SD}/fig1_asymmetry.csv")
axA = fig.add_axes([0.085, 0.565, 0.37, 0.37])
fig_letter(fig, "A", 0.005, 0.995)
GROUPS = [("silent", "0"), ("weak", "1–4"), ("active_T1", "low"), ("active_T2", "mid"), ("active_T3", "high")]
YMIN, YMAX = -100, 100
for line, col, dx in LINES:
    m = mt[(mt["line"] == line) & (mt["comparison"] == "matched")].iloc[0]
    axA.axhline(m["sense_weaker_pct"], color=col, lw=0.7, ls="--", zorder=0)
    for i, (g, _) in enumerate(GROUPS):
        r = ha[(ha["line"] == line) & (ha["group"] == g)].iloc[0]
        lo = max(r["pct_lo95"], YMIN)
        ci_point(axA, i + dx, r["sense_weaker_pct"], lo, min(r["pct_hi95"], YMAX), col, ms=3.4)
        if r["pct_lo95"] < YMIN:
            axA.plot([i + dx], [YMIN + 2], marker="v", ms=3.2, color=col, zorder=4)
axA.axhline(0, color=GREY, lw=0.6, zorder=0)
axA.axvspan(1.5, 4.5, color="#F5F5F5", zorder=-1, lw=0)
axA.set_xticks(range(len(GROUPS)))
axA.set_xticklabels([g[1] for g in GROUPS])
axA.set_xlim(-0.5, len(GROUPS) - 0.5)
axA.set_ylim(YMIN, YMAX)
axA.set_xlabel("activity of the host gene (reads at its start: silent, weak,\nand thirds of pairs in active hosts)")
axA.set_ylabel("sense start weaker than in matched\nintergenic pairs (%, 95% CI)")
axA.text(3.0, YMIN + 6, "active hosts", ha="center", va="bottom", fontsize=6.3, color=GREY)
txt = "; ".join(f"{r['line']} {r['spearman_rho_hostsignal_vs_withinpair_diff']:+.2f} "
                f"({r['rho_lo95']:+.2f} to {r['rho_hi95']:+.2f})" for _, r in tr.iterrows())
txt = txt.replace("-", "−")
axA.text(0.0, 1.03, f"Spearman ρ with host activity: {txt}", transform=axA.transAxes, fontsize=6.2, color=GENE,
         ha="left", va="bottom")
for line, col, yy in (("K562", K562, 0.97), ("HCT116", HCT116, 0.90)):
    axA.text(0.99, yy, line, transform=axA.transAxes, color=col, fontsize=6.8, fontweight="bold", ha="right", va="top")

# ------------------------------------------------------------------ B: share predicted from sequence
sh = pd.read_csv(f"{SD}/fig2_shares.csv")
sh = sh[sh["group"] == "intron_active"]
ROWS = [("composition_300nt_mono", "base composition only\n(within 300 nt)", "D"),
        ("puffin_D", "Puffin-D", "s"),
        ("pc_K562_heldout", "ProCapNet trained on K562,\nnever saw the chromosome", "o"),
        ("pc_Caco2_heldout", "ProCapNet trained on Caco-2,\nnever saw the chromosome", "o")]
axB = fig.add_axes([0.70, 0.565, 0.28, 0.37])
fig_letter(fig, "B", 0.50, 0.995)
for i, (lab, name, mk) in enumerate(ROWS):
    yy = len(ROWS) - 1 - i
    for line, col, dy in LINES:
        r = sh[(sh["label"] == lab) & (sh["line"] == line)].iloc[0]
        ci_point(axB, 100 * r["frac_linear"], yy - dy, 100 * r["frac_linear_lo95"], 100 * r["frac_linear_hi95"], col,
                 marker=mk, ms=3.6 if mk != "D" else 3.2, horizontal=True)
    axB.text(-58, yy, name, ha="right", va="center", fontsize=6.6, linespacing=1.1)
axB.axvline(0, color=GREY, lw=0.6, zorder=0)
axB.set_yticks([])
axB.spines["left"].set_visible(False)
axB.set_ylim(-0.6, len(ROWS) - 0.4)
axB.set_xlim(-55, 105)
axB.set_xticks([-50, 0, 50, 100])
axB.set_xlabel("share of the asymmetry\npredicted from sequence (%, 95% CI)")

# ------------------------------------------------------------------ C: where in the sequence (pair-level shuffles)
REG = [("shuf_none", "original sequence"), ("shuf_window", "whole window shuffled"),
       ("shuf_element", "element shuffled"), ("shuf_flanks", "flanks shuffled,\nelement intact")]
axC = fig.add_axes([0.0, 0.11, 0.60, 0.34])
fig_letter(fig, "C", 0.005, 0.47)
axC.set_xlim(0, 100)
axC.set_ylim(-0.75, len(REG) - 0.05)
axC.axis("off")
W0, W1 = 25.0, 98.0
E0, E1 = W0 + (W1 - W0) * 0.41, W0 + (W1 - W0) * 0.59
AS, SS = E0 + (E1 - E0) * 0.22, E0 + (E1 - E0) * 0.78
FILL, EDGE, BH, HATCH = "#BDBDBD", "#404040", 0.22, "//////"
for i, (reg, lab) in enumerate(REG):
    yy = len(REG) - 1 - i
    shuf = {"shuf_none": (False, False, False), "shuf_window": (True, True, True),
            "shuf_element": (False, True, False), "shuf_flanks": (True, False, True)}[reg]
    for (lo_, hi_), sh_ in zip(((W0, E0), (E0, E1), (E1, W1)), shuf):
        p = box(axC, lo_, yy - BH, hi_ - lo_, 2 * BH, fc="white" if sh_ else FILL, ec=EDGE, lw=0.6, z=1)
        if sh_:
            p.set_hatch(HATCH)
            p.set_hatch_linewidth(0.5)
    for xa, col in ((AS, ANTI), (SS, SENSE)):
        axC.plot([xa, xa], [yy - 0.32, yy + 0.32], color=col, lw=1.8, zorder=3)
    axC.text(W0 - 1.5, yy, lab, ha="right", va="center", fontsize=6.6, linespacing=1.1)
yt = len(REG) - 1 + 0.45
arrow(axC, SS, yt, SS + 4.5, yt, SENSE, lw=1.0, head=4)
arrow(axC, AS, yt, AS - 4.5, yt, ANTI, lw=1.0, head=4)
axC.plot([E0, E1], [yt + 0.2, yt + 0.2], color=GENE, lw=0.6)
axC.text((E0 + E1) / 2, yt + 0.26, "element", ha="center", va="bottom", fontsize=6.3, color=GENE)
axC.text(W0, yt + 0.26, "2,114 nt seen by ProCapNet", ha="left", va="bottom", fontsize=6.3, color=GENE)
YL = -0.62
for x_, sh_, lab in ((W0, False, "original order"), (49.0, True, "shuffled (counts of bases and dinucleotides kept)")):
    p = box(axC, x_, YL - BH, 4.0, 2 * BH, fc="white" if sh_ else FILL, ec=EDGE, lw=0.6, z=1)
    if sh_:
        p.set_hatch(HATCH)
        p.set_hatch_linewidth(0.5)
    axC.text(x_ + 4.0 + 1.0, YL, lab, ha="left", va="center", fontsize=6.3, color=GENE)

axD = fig.add_axes([0.66, 0.11, 0.31, 0.34])
for i, (reg, _) in enumerate(REG):
    yy = len(REG) - 1 - i
    for line, col, dy in LINES:
        r = sh[(sh["label"] == reg) & (sh["line"] == line)].iloc[0]
        ci_point(axD, 100 * r["frac_linear"], yy - dy, 100 * r["frac_linear_lo95"], 100 * r["frac_linear_hi95"],
                 col, marker="o", ms=3.4, horizontal=True)
k0 = 100 * sh[(sh["label"] == "shuf_none") & (sh["line"] == "K562")]["frac_linear"].iloc[0]
axD.axvline(0, color=GREY, lw=0.6, zorder=0)
axD.axvline(k0, color=K562, lw=0.6, ls=":", zorder=0)
axD.set_yticks([])
axD.spines["left"].set_visible(False)
axD.set_ylim(-0.75, len(REG) - 0.05)
axD.set_xlim(-25, 75)
axD.set_xlabel("share predicted by the K562 ProCapNet models\nthat never saw the chromosome (%, 95% CI)")
for line, col, yy in (("K562", K562, 0.62), ("HCT116", HCT116, 0.52)):
    axD.text(0.98, yy, line, transform=axD.transAxes, color=col, fontsize=6.8, fontweight="bold", ha="right")
save(fig, OUT)
