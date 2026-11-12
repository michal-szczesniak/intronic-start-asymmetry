#!/usr/bin/env python3
"""Make 19b_selection_classes.py from 12c_selection_genomewide.py: the same test with the class
definitions (site thresholds, change cut-offs, chromosome groups) as parameters. Every edit is an
exact text replacement that must match the expected number of times; nothing else changes.
Usage: python3 patch_12c_to_19b.py <12c script> <19b script to write>"""
import sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
EDITS = [
    ('"""Step 12c: the selection test of steps 12/12b on many chromosomes, read one at a time.\n',
     '"""Step 19b: step 12c with the class definitions as parameters, so that the variant classes can come\n'
     'from a model on another scale (ProCapNet, step 19a). With the default parameters this is step 12c\n'
     'exactly (byte-identical output tables; checked on the 12c test data).\n\n'
     'Parameters that define the classes (defaults = Puffin-D, step 12c):\n'
     '  --thresholds 0.0439,0.1432  site thresholds (a variant creates / destroys a site)\n'
     '  --neutral 0.001             no-effect class: change below this on both strands\n'
     '  --mid 0.003, --strong 0.01  moderate / strong gain or loss on one strand ...\n'
     '  --other 0.003               ... while the other strand changes less than this\n'
     '  --cutoffs-file              step-19c table; takes neutral, mid, strong and other from it\n'
     '  --groups                    named groups of chromosomes for groups.tsv (default: heldout_test\n'
     '                              chr8,chr9; heldout_valid chr10; the rest "training"); for ProCapNet\n'
     '                              its seven folds\n'
     'Also written: class_definitions.tsv. The rest of this text is step 12c\'s.\n\n'
     'Step 12c: the selection test of steps 12/12b on many chromosomes, read one at a time.\n', 1),
    ('GROUPS = (("heldout_test", ("chr8", "chr9")), ("heldout_valid", ("chr10",)))\n',
     'GROUPS = (("heldout_test", ("chr8", "chr9")), ("heldout_valid", ("chr10",)))\n\n\n'
     'def set_classes(thresholds, strong, mid, other):\n'
     '    """Rebuild the class list for other thresholds (same order and names as above)."""\n'
     '    global THRESHOLDS, STRONG, MID, OTHER, CLASSES, BIT, PAIRS\n'
     '    THRESHOLDS, STRONG, MID, OTHER = tuple(thresholds), strong, mid, other\n'
     '    CLASSES = ["neutral"]\n'
     '    for T in THRESHOLDS:\n'
     '        CLASSES += [f"creates_sense_{T}", f"creates_anti_{T}", f"destroys_sense_{T}", f"destroys_anti_{T}"]\n'
     '    CLASSES += ["gain_sense_strong", "gain_anti_strong", "gain_sense_mid", "gain_anti_mid",\n'
     '                "loss_sense_strong", "loss_anti_strong"]\n'
     '    for T in THRESHOLDS:\n'
     '        CLASSES += [f"creates_sense_{T}_unique", f"creates_anti_{T}_unique"]\n'
     '    BIT = {k: np.uint32(1 << i) for i, k in enumerate(CLASSES)}\n'
     '    PAIRS = [(f"creates_sense_{T}", f"creates_anti_{T}") for T in THRESHOLDS] + \\\n'
     '            [(f"creates_sense_{T}_unique", f"creates_anti_{T}_unique") for T in THRESHOLDS] + \\\n'
     '            [(f"destroys_sense_{T}", f"destroys_anti_{T}") for T in THRESHOLDS] + \\\n'
     '            [("gain_sense_strong", "gain_anti_strong"), ("gain_sense_mid", "gain_anti_mid"),\n'
     '             ("loss_sense_strong", "loss_anti_strong")]\n', 1),
    ('    ap.add_argument("--neutral", type=float, default=0.001)\n',
     '    ap.add_argument("--neutral", type=float, default=0.001)\n'
     '    ap.add_argument("--thresholds", default="0.0439,0.1432", help="the two site thresholds")\n'
     '    ap.add_argument("--strong", type=float, default=0.01)\n'
     '    ap.add_argument("--mid", type=float, default=0.003)\n'
     '    ap.add_argument("--other", type=float, default=0.003)\n'
     '    ap.add_argument("--cutoffs-file", default="", help="19c table (overrides neutral, mid, strong, other)")\n'
     '    ap.add_argument("--groups", default="heldout_test:chr8,chr9;heldout_valid:chr10")\n'
     '    ap.add_argument("--rest-name", default="training", help="name of the group of all other chromosomes")\n', 1),
    ('    chroms = [c for c in a.chroms.split(",") if c]\n',
     '    thr = tuple(float(x) for x in a.thresholds.split(",") if x)\n'
     '    if len(thr) != 2:\n'
     '        sys.exit("STOP: --thresholds needs two values")\n'
     '    neutral, strong, mid, other = a.neutral, a.strong, a.mid, a.other\n'
     '    if a.cutoffs_file:\n'
     '        cf = pd.read_csv(a.cutoffs_file, sep="\\t").set_index("quantity")["procapnet_cutoff"]\n'
     '        neutral, mid, strong, other = (float(cf[k]) for k in ("neutral", "mid", "strong", "other"))\n'
     '    set_classes(thr, strong, mid, other)\n'
     '    groups = tuple((g.split(":")[0], tuple(x for x in g.split(":")[1].split(",") if x))\n'
     '                   for g in a.groups.split(";") if g)\n'
     '    print(f"classes: site thresholds {THRESHOLDS[0]}, {THRESHOLDS[1]}; no effect: change < {neutral:g} on both "\n'
     '          f"strands; moderate >= {MID:g}, strong >= {STRONG:g} with the other strand < {OTHER:g}; groups: "\n'
     '          + "; ".join(f"{g}={\',\'.join(cs)}" for g, cs in groups), flush=True)\n'
     '    chroms = [c for c in a.chroms.split(",") if c]\n', 1),
    ('    os.makedirs(a.out_dir, exist_ok=True)\n',
     '    os.makedirs(a.out_dir, exist_ok=True)\n'
     '    pd.DataFrame([{"quantity": "site_threshold_1", "value": THRESHOLDS[0]},\n'
     '                  {"quantity": "site_threshold_2", "value": THRESHOLDS[1]},\n'
     '                  {"quantity": "neutral", "value": neutral}, {"quantity": "mid", "value": MID},\n'
     '                  {"quantity": "strong", "value": STRONG}, {"quantity": "other", "value": OTHER}]\n'
     '                 ).to_csv(os.path.join(a.out_dir, "class_definitions.tsv"), sep="\\t", index=False)\n', 1),
    ('                                           a.chunk, a.neutral)\n',
     '                                           a.chunk, neutral)\n', 1),
    ('    group_defs = [(c, (c,)) for c in chroms] + [g for g in GROUPS if any(x in chroms for x in g[1])] + \\\n'
     '                 [("training", tuple(c for c in chroms if c not in {x for g in GROUPS for x in g[1]}))]\n',
     '    group_defs = [(c, (c,)) for c in chroms] + [g for g in groups if any(x in chroms for x in g[1])] + \\\n'
     '                 [(a.rest_name, tuple(c for c in chroms if c not in {x for g in groups for x in g[1]}))]\n', 1),
    ('            for k in ("splice_canonical", "splice_donor_3to6", "creates_sense_0.0439", "creates_anti_0.0439",\n'
     '                      "gain_sense_strong", "gain_anti_strong", "destroys_sense_0.0439", "destroys_anti_0.0439",\n',
     '            T0 = THRESHOLDS[0]\n'
     '            for k in ("splice_canonical", "splice_donor_3to6", f"creates_sense_{T0}", f"creates_anti_{T0}",\n'
     '                      "gain_sense_strong", "gain_anti_strong", f"destroys_sense_{T0}", f"destroys_anti_{T0}",\n', 1),
    ('            for k1 in ("creates_sense_0.0439", "creates_sense_0.1432", "gain_sense_strong", "destroys_sense_0.0439"):\n',
     '            for k1 in (f"creates_sense_{THRESHOLDS[0]}", f"creates_sense_{THRESHOLDS[1]}", "gain_sense_strong",\n'
     '                       f"destroys_sense_{THRESHOLDS[0]}"):\n', 1),
    ('                ("intron_protein_coding", "diff creates_sense_0.0439 - anti"),\n',
     '                ("intron_protein_coding", f"diff creates_sense_{THRESHOLDS[0]} - anti"),\n', 1),
    ('                ("intergenic", "diff creates_sense_0.0439 - anti")]\n',
     '                ("intergenic", f"diff creates_sense_{THRESHOLDS[0]} - anti")]\n', 1),
]
for old, new, times in EDITS:
    k = s.count(old)
    if k != times:
        sys.exit(f"STOP: expected {times} match(es), found {k}: {old[:70]!r}")
    s = s.replace(old, new)
open(dst, "w").write(s)
print(f"written {dst} ({len(EDITS)} edits)")
