"""Shared style for the PURGE figures: sans serif 6.5-8 pt, bold 10-pt panel letters,
open axes, Okabe-Ito colours with fixed meanings. Full page width 180 mm."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

MM = 1 / 25.4
WIDTH = 180 * MM

# Okabe-Ito, one meaning each
K562 = "#0072B2"        # blue
HCT116 = "#D55E00"      # vermilion
PUFFIN = "#CC79A7"      # reddish purple
PROCAP = "#009E73"      # bluish green
SENSE = "#E69F00"       # orange
ANTI = "#56B4E9"        # sky blue
OBS = "#000000"         # measured data
# neutrals
GREY = "#7F7F7F"
MID = "#A6A6A6"
LIGHT = "#D9D9D9"
PALE = "#F0F0F0"
GENE = "#404040"
MARK = {"Puffin-D": "s", "ProCapNet": "o"}   # squares for Puffin-D, circles for ProCapNet, everywhere

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
    "mathtext.default": "regular",
    "lines.linewidth": 1.0,
})


def fig_letter(fig, s, x, y):
    """Bold 10-pt panel letter at figure coordinates (top-left anchored)."""
    fig.text(x, y, s, fontsize=10, fontweight="bold", va="top", ha="left")


def arrow(ax, x0, y0, x1, y1, color, lw=1.2, head=6, style="-|>", z=3, ls="-"):
    st = style if style == "-" else f"{style},head_length={head * 0.12},head_width={head * 0.07}"
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=st, color=color, lw=lw, zorder=z, shrinkA=0, shrinkB=0,
                        linestyle=ls, mutation_scale=10)
    ax.add_patch(a)
    return a


def box(ax, x, y, w, h, fc, ec="none", lw=0.6, r=0.0, z=2, alpha=1.0, ls="-"):
    if r:
        p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z,
                           alpha=alpha, ls=ls)
    else:
        p = Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z, alpha=alpha, ls=ls)
    ax.add_patch(p)
    return p


def ci_point(ax, x, y, lo, hi, color, marker="o", ms=3.6, horizontal=False, filled=True, lw=0.9, z=3, label=None):
    """Point with a 95% CI bar; open symbols get a white face. Vertical: bar from lo to hi at x.
    Horizontal: x is the value, y the row, and the bar runs from lo to hi along x."""
    mfc = color if filled else "white"
    if horizontal:
        ax.plot([lo, hi], [y, y], color=color, lw=lw, zorder=z, solid_capstyle="butt")
    else:
        ax.plot([x, x], [lo, hi], color=color, lw=lw, zorder=z, solid_capstyle="butt")
    ax.plot([x], [y], marker=marker, ms=ms, mfc=mfc, mec=color, mew=0.8, ls="none", zorder=z + 1, label=label)


def log2_ratio_axis(ax, ticks=(0.5, 0.75, 1, 1.5, 2), axis="y"):
    """Strand ratios are drawn on a log2 scale with ratio labels (values below 1: fewer on the gene strand)."""
    t = np.log2(np.asarray(ticks, float))
    lab = [f"{v:g}×" for v in ticks]
    if axis == "y":
        ax.set_yticks(t)
        ax.set_yticklabels(lab)
    else:
        ax.set_xticks(t)
        ax.set_xticklabels(lab)


def save(fig, path):
    fig.savefig(path + ".png", dpi=600)
    fig.savefig(path + ".pdf")
    plt.close(fig)
