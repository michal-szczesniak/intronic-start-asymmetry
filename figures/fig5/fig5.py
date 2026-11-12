#!/usr/bin/env python3
"""Figure 5: the conclusions at three scales. Schematic; no data are plotted.

The numbers quoted in the text blocks are taken from the following results (Results 1, 2, 4 and 5):
  33% (K562) and 21% (HCT116), depth-matched bias (Fig. 1B);
  30-42%, share of the bias predicted by held-out sequence models (Fig. 2B);
  15-17%, deficit of strong measured sites in active genes of HCT116 and K562 (step 28, Results 4);
  46-56%, the same deficit within 5 kb of the 5' end of the gene (step 28c, Supplemental Fig. S5);
  one hundredth of the excess of singletons at canonical splice sites (Fig. 4C);
  76,215 genomes, gnomAD v4.1 (Methods).

Usage: python3 fig5.py [output_prefix]   (default: figures/Szczesniak_Fig5)
Writes <prefix>.png (600 dpi), <prefix>.pdf and <prefix>.svg (live text, for editing in a vector program).
Needs figstyle.py (shared style of the PURGE figures) and the fonts/ folder (Liberation Sans, SIL Open Font
License) in the same folder; the script stops if Liberation Sans cannot be found instead of silently switching
to a wider font. It prints the matplotlib and FreeType versions and the font files it uses.
The script stops with an error if any text runs outside its band, so a later wording change cannot
silently produce an overlapping figure.

Output is byte-for-byte reproducible with the same matplotlib and FreeType versions (no creation dates in the
PDF and SVG, fixed SVG ids); another FreeType version can change the bytes of the PNG without changing the figure.

Rebuilt on 10 Nov 2026 from the image of an earlier version of Figure 5 (30 Oct 2026), whose script was not
kept; changes: scale labels
spelled out ("Scale: ..."), positive wording of the population-test limit, "Supplementary" -> "Supplemental",
text column widened and band height derived from the text, label of the element moved clear of the ticks.
Version b (10 Nov 2026): fonts shipped with the script, because Liberation Sans was missing on the cluster.
Version c (10 Nov 2026): panel 2 redrawn after steps 28-28f: strong sites concentrated within 5 kb of the start
of the gene, a wedge for H3K36me3 rising toward the 3' end, and new text.
Version d (10 Nov 2026): supplemental figures renumbered in the order of citation; the position figure is now
Supplemental Fig. S5 (was S8).
"""
import hashlib
import os
import sys

import matplotlib
from matplotlib import font_manager as fm
from matplotlib import ft2font
from matplotlib.font_manager import FontProperties

from figstyle import ANTI, GENE, GREY, MM, PROCAP, SENSE, arrow, box, plt

HERE = os.path.dirname(os.path.abspath(__file__))
FONT_DIR = os.path.join(HERE, "fonts")


def md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


# Fonts. Liberation Sans (the font of all PURGE figures) ships with the script in fonts/, because it is missing on
# many Linux servers and matplotlib would otherwise fall back to the wider DejaVu Sans. The DNA letters use DejaVu
# Sans Mono, which comes with matplotlib.
if os.path.isdir(FONT_DIR):
    for f_ in sorted(os.listdir(FONT_DIR)):
        if f_.lower().endswith(".ttf"):
            fm.fontManager.addfont(os.path.join(FONT_DIR, f_))
SANS = {}
for w_ in ("normal", "bold"):
    try:
        SANS[w_] = fm.findfont(FontProperties(family="Liberation Sans", weight=w_), fallback_to_default=False)
    except ValueError:
        sys.exit("STOP: Liberation Sans not found; keep the fonts/ folder next to fig5.py")
MONO = {w_: fm.findfont(FontProperties(family="monospace", weight=w_)) for w_ in ("normal", "bold")}
print(f"matplotlib {matplotlib.__version__}, FreeType {ft2font.__freetype_version__}")
for p_ in list(SANS.values()) + list(MONO.values()):
    print(f"  font {os.path.basename(p_)}  md5 {md5(p_)}  ({p_})")

OUT = sys.argv[1] if len(sys.argv) > 1 else "figures/Szczesniak_Fig5"

# ------------------------------------------------------------------------------------------------ layout (mm)
W = 180.0
BG = "#F2F2F2"
INK = "#222222"
SCALE_GREY = "#737373"
FS_BODY, FS_HEAD, FS_SCALE, FS_LAB = 8.0, 8.5, 8.0, 8.0
LS_BODY = 1.35                           # line spacing of the body text
TX0, TX1 = 96.0, 177.0                   # text column
PAD_TOP, PAD_BOT = 3.0, 2.5              # inside a band, above the badge and below the last line of text
BAND_H, GAP, MARGIN = 37.0, 2.5, 1.0     # band height, gap between bands, outer margin
SCHEM_H = 33.0                           # height of the schematic drawings; centred in each band
H = 2 * MARGIN + 3 * BAND_H + 2 * GAP

fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")
renderer = fig.canvas.get_renderer()
PT = 25.4 / 72                           # mm per point

BANDS = [H - MARGIN - BAND_H - i * (BAND_H + GAP) for i in range(3)]   # bottom y of each band
OFF = (BAND_H - SCHEM_H) / 2             # vertical offset that centres a schematic in its band
TEXTS = {0: [], 1: [], 2: []}            # text artists of each band, for the fit check at the end


def T(band, *args, **kw):
    """ax.text that also records the artist for the fit check."""
    t = ax.text(*args, **kw)
    TEXTS[band].append(t)
    return t


def text_width_mm(s, size, weight="normal"):
    t = ax.text(0, 0, s, fontsize=size, fontweight=weight)
    bb = t.get_window_extent(renderer=renderer)
    t.remove()
    return bb.width / fig.dpi * 25.4


def wrap(s, size, width_mm, weight="normal"):
    """Greedy word wrap measured with the real font metrics; returns the lines."""
    lines, cur = [], ""
    for word in s.split():
        trial = (cur + " " + word).strip()
        if text_width_mm(trial, size, weight) <= width_mm or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = word
    lines.append(cur)
    return lines


def text_block(band, number, scale, heading, body):
    """Number badge, scale label (grey), bold heading and wrapped body text in the right-hand column."""
    y_top = BANDS[band] + BAND_H
    yb = y_top - PAD_TOP - 2.0
    ax.add_patch(plt.Circle((TX0 + 2.0, yb), 2.0, color="#333333", zorder=3))
    T(band, TX0 + 2.0, yb - 0.05, str(number), ha="center", va="center", fontsize=8.5, fontweight="bold",
      color="white", zorder=4)
    T(band, TX0 + 5.4, yb, scale, ha="left", va="center", fontsize=FS_SCALE, color=SCALE_GREY)
    head = wrap(heading, FS_HEAD, TX1 - TX0, "bold")
    y = yb - 4.5
    T(band, TX0, y, "\n".join(head), ha="left", va="top", fontsize=FS_HEAD, fontweight="bold", color=INK,
      linespacing=1.2)
    y -= len(head) * FS_HEAD * 1.2 * PT + 1.6
    T(band, TX0, y, "\n".join(wrap(body, FS_BODY, TX1 - TX0)), ha="left", va="top", fontsize=FS_BODY,
      color=INK, linespacing=LS_BODY)


for yb_ in BANDS:
    box(ax, MARGIN, yb_, W - 2 * MARGIN, BAND_H, fc=BG, r=1.5, z=0)

# ------------------------------------------------------------------------------------------------ 1: one element
y0 = BANDS[0] + OFF
yt, ybot = y0 + 19.5, y0 + 14.5                      # the two DNA strands; the upper one is the gene strand
ax.plot([20, 90], [yt, yt], color=GENE, lw=1.0)
ax.plot([20, 90], [ybot, ybot], color=GENE, lw=1.0)
T(0, 18.0, (yt + ybot) / 2, "gene\nstrand", ha="right", va="center", fontsize=FS_LAB, color=GREY, linespacing=1.1)
ex0, ex1 = 43.0, 71.0                                # the element
box(ax, ex0, ybot - 1.2, ex1 - ex0, yt - ybot + 2.4, fc="white", ec=GREY, lw=0.7, r=0.8, z=1)
for xf in (45.5, 49.5, 60.5, 64.5):                  # sequence features that the models read
    box(ax, xf - 1.1, (yt + ybot) / 2 - 1.1, 2.2, 2.2, fc=PROCAP, r=0.25, z=2, alpha=0.8)
# sense start: thin, on the gene strand, pointing in the direction of the gene
xs = 58.0
arrow(ax, xs, yt, xs, yt + 4.2, SENSE, lw=1.1, style="-", z=3)
arrow(ax, xs, yt + 4.2, xs + 6.5, yt + 4.2, SENSE, lw=1.1, head=6, z=3)
T(0, xs + 7.4, yt + 4.2, "sense start\n(weaker)", ha="left", va="center", fontsize=FS_LAB, fontweight="bold",
  color=SENSE, linespacing=1.0)
# antisense start: thick, on the other strand, pointing the other way
xa = 53.0
arrow(ax, xa, ybot, xa, ybot - 6.0, ANTI, lw=3.0, style="-", z=3)
arrow(ax, xa, ybot - 6.0, xa - 10.0, ybot - 6.0, ANTI, lw=3.0, head=9, z=3)
T(0, xa - 11.5, ybot - 6.0, "antisense\nstart", ha="right", va="center", fontsize=FS_LAB, fontweight="bold",
  color=ANTI, linespacing=1.0)
T(0, 55.5, ybot - 4.2, "sequence features\npredict 30–42%", ha="left", va="top", fontsize=FS_LAB,
  fontweight="bold", color=PROCAP, linespacing=1.0)
text_block(0, 1, "Scale: one element, a few hundred nucleotides long",
           "Weaker sense start, partly predicted from its sequence",
           "The sense start is weaker by a fifth to a third than at elements outside genes (33% in K562, "
           "21% in HCT116; Fig. 1B). Models that never saw the chromosome predict 30–42% of this from the "
           "element's own sequence (Fig. 2B), from the order of its bases (Fig. 2C).")

# ------------------------------------------------------------------------------------------------ 2: the host gene
y0 = BANDS[1] + OFF
yg = y0 + 15.5
ax.plot([5, 91], [yg, yg], color=GENE, lw=1.0)
for a_, b_ in ((7.0, 11.0), (31.0, 34.0), (68.0, 71.0), (84.0, 90.0)):    # exons
    box(ax, a_, yg - 1.4, b_ - a_, 2.8, fc=GENE, r=0.3, z=2)
arrow(ax, 7.0, yg + 1.4, 7.0, yg + 4.2, GENE, lw=1.1, style="-")
arrow(ax, 7.0, yg + 4.2, 11.5, yg + 4.2, GENE, lw=1.1, head=6)
T(1, 12.5, yg + 9.8, "active host gene", ha="left", va="center", fontsize=FS_LAB, color=GENE)
xe = 50.0                                                                # the element of panel 1
box(ax, xe - 1.2, yg - 1.1, 2.4, 2.2, fc="white", ec=GREY, lw=0.7, z=3)
T(1, xe, yg - 2.6, "element", ha="center", va="top", fontsize=FS_LAB, color=GREY)
NEAR0, NEAR1 = 12.0, 25.0                                                # the first 5 kb after the start
box(ax, NEAR0, yg - 7.4, NEAR1 - NEAR0, 14.8, fc="#E0E0E0", r=0.4, z=0)
T(1, (NEAR0 + NEAR1) / 2, yg - 8.2, "first 5 kb", ha="center", va="top", fontsize=FS_LAB, color=GREY)
for xs_ in (19.0, 41.0, 61.0, 79.0):                                     # strong sites on the gene strand
    ax.plot([xs_, xs_], [yg + 1.6, yg + 6.2], color=SENSE, lw=1.6, solid_capstyle="butt")
for xa_ in (14.0, 16.5, 20.5, 23.0, 44.0, 64.0, 81.5):                   # on the opposite strand
    ax.plot([xa_, xa_], [yg - 1.6, yg - 6.2], color=ANTI, lw=1.6, solid_capstyle="butt")
wy = yg + 8.6                                                            # H3K36me3: rises toward the 3' end
ax.fill([38.0, 90.0, 90.0], [wy, wy, wy + 3.2], color="#BDBDBD", lw=0, zorder=1)
T(1, 90.0, wy + 5.2, "H3K36me3", ha="right", va="center", fontsize=FS_LAB, color=GREY)
T(1, 48.0, y0 + 2.4, "strong start sites: fewer on the gene strand near the start", ha="center", va="center",
  fontsize=FS_LAB, color=GREY)
text_block(1, 2, "Scale: the host gene, tens of kilobases long",
           "Fewer strong sense sites near the start of active genes",
           "Strong start sites are less frequent on the gene strand in both cell lines (15–17% in active "
           "genes), mainly within 5 kb of the start of the gene (46–56%), not toward the 3′ end, where "
           "H3K36me3 is highest (Supplemental Fig. S5). The weaker sense start of the pairs is also largest "
           "near the start.")

# ------------------------------------------------------------------------------------------------ 3: all introns
y0 = BANDS[2] + OFF
# left: potential starts, as many on each strand
yt, ybot = y0 + 19.0, y0 + 16.5
ax.plot([8, 44], [yt, yt], color=GENE, lw=1.0)
ax.plot([8, 44], [ybot, ybot], color=GENE, lw=1.0)
arrow(ax, 20.0, yt, 20.0, yt + 4.0, SENSE, lw=1.6, style="-", z=3)
arrow(ax, 20.0, yt + 4.0, 26.0, yt + 4.0, SENSE, lw=1.6, head=7, z=3)
arrow(ax, 31.0, ybot, 31.0, ybot - 4.0, ANTI, lw=1.6, style="-", z=3)
arrow(ax, 31.0, ybot - 4.0, 25.0, ybot - 4.0, ANTI, lw=1.6, head=7, z=3)
T(2, 26.0, y0 + 4.0, "potential starts in introns:\nas many on each strand", ha="center", va="center",
  fontsize=FS_LAB, color=GREY, linespacing=1.15)
# right: a stretch of DNA with two variants
seq = "GACTTAGCA"
xs0, dx, ys = 55.0, 4.0, y0 + 17.5
for i, base in enumerate(seq):
    x = xs0 + i * dx
    var = {2: "T", 6: "A"}.get(i)
    if var:
        box(ax, x - 1.5, ys - 2.0, 3.0, 4.0, fc="white", ec=INK, lw=0.7, r=0.3, z=2)
        ax.plot([x, x], [ys + 2.0, ys + 3.4], color=INK, lw=0.7)
        box(ax, x - 1.5, ys + 3.4, 3.0, 4.0, fc="#333333", r=0.3, z=2)
        T(2, x, ys + 5.4, var, ha="center", va="center", fontsize=FS_LAB, fontweight="bold",
          family="monospace", color="white", zorder=3)
    T(2, x, ys, base, ha="center", va="center", fontsize=FS_LAB + 1, family="monospace",
      color=INK if var else GREY, fontweight="bold" if var else "normal", zorder=3)
T(2, xs0 + 4 * dx, y0 + 4.0, "variants that create sense\nstarts: not rarer", ha="center", va="center",
  fontsize=FS_LAB, color=GREY, linespacing=1.15)
text_block(2, 3, "Scale: all introns; variants in 76,215 genomes",
           "No sign of selection against sense starts",
           "Introns hold about as many potential sense as antisense starts (ProCapNet; Fig. 3B). Variants "
           "that create a sense start are not kept rarer than those that create an antisense start; the test "
           "would have detected a difference of one hundredth of the splice-site signal (Fig. 4C).")

# ------------------------------------------------------------------------------------------------ fit check
fig.canvas.draw()
to_mm = ax.transData.inverted()
problems = []
for band, texts in TEXTS.items():
    lo, hi = BANDS[band] + PAD_BOT, BANDS[band] + BAND_H - 1.0
    for t in texts:
        (x0, y0_), (x1, y1_) = to_mm.transform(t.get_window_extent(renderer=fig.canvas.get_renderer()))
        right_limit = TX1 + 0.3 if x0 >= TX0 - 0.1 else TX0 - 3.0     # schematic labels must stay clear of the text
        if x0 < MARGIN + 1.0 or x1 > right_limit or y0_ < lo or y1_ > hi:
            problems.append(f"band {band + 1}: {t.get_text()[:40]!r} at x {x0:.1f}-{x1:.1f}, y {y0_:.1f}-{y1_:.1f}")
if problems:
    sys.exit("text outside its band:\n" + "\n".join(problems))

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
plt.rcParams["svg.hashsalt"] = "PURGE-Fig5"          # fixed ids in the SVG
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf", metadata={"CreationDate": None, "ModDate": None})
fig.savefig(OUT + ".svg", metadata={"Date": None})
plt.close(fig)
print(f"wrote {OUT}.png/.pdf/.svg ({W:.0f} x {H:.0f} mm)")
