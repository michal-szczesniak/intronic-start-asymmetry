#!/usr/bin/env python3
"""Compares each 600-dpi figure in figures/ with the image of the same figure in the manuscript
(manuscript_images/, taken from the manuscript on 11 Nov 2026; the manuscript holds each figure at half
the pixel size of the 600-dpi file). Each figure is reduced to the size of the manuscript image (Lanczos filter,
Pillow) and compared pixel by pixel.

Usage: python3 check_figures.py [figures_folder]   (default: figures)
Prints, per figure, the number of pixels that differ and the largest difference (0-255); stops with an error if
any figure differs.
"""
import os
import sys

import numpy as np
from PIL import Image

FIG = sys.argv[1] if len(sys.argv) > 1 else "figures"
HERE = os.path.dirname(os.path.abspath(__file__))
bad = 0
for name in sorted(os.listdir(os.path.join(HERE, "manuscript_images"))):
    t = Image.open(os.path.join(HERE, "manuscript_images", name)).convert("RGB")
    p = os.path.join(FIG, name)
    if not os.path.exists(p):
        print(f"{name:24s} MISSING in {FIG}")
        bad += 1
        continue
    r = Image.open(p).convert("RGB").resize(t.size, Image.LANCZOS)
    d = np.abs(np.asarray(t).astype(int) - np.asarray(r).astype(int)).max(axis=2)
    print(f"{name:24s} {t.size[0]}x{t.size[1]}  differing pixels {int((d > 0).sum()):,}  largest difference {int(d.max())}")
    bad += int(d.max() > 0)
if bad:
    sys.exit(f"STOP: {bad} figure(s) differ from the manuscript images")
print("all figures identical to the manuscript images")
