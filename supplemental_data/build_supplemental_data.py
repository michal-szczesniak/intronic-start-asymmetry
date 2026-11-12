#!/usr/bin/env python3
"""Supplemental Data of the PURGE manuscript: every iSTART start-site pair of K562 and HCT116
with the per-pair measures used in the study, joined from the result tables of the analysis steps.

Every input is a result table of the project, read unchanged (md5 printed and written to the README sheet). Pairs are
identified by (line, chrom, mid); every per-pair table must hold exactly the pairs of the step-04 table, in the same
order, with the same read counts and DI, or the script stops. Nothing is recomputed except depth = sense + antisense.

Usage: python3 build_supplemental_data.py --root <PURGE project folder> --out Supplemental_Data.xlsx [--tsv pairs.tsv.gz]
"""
import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd

LINES = ["K562", "HCT116"]
KEY = ["line", "chrom", "mid"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Src:
    def __init__(self, root):
        self.root, self.used = root, {}

    def read(self, rel):
        p = os.path.join(self.root, rel)
        if not os.path.exists(p):
            sys.exit(f"STOP: missing {rel}")
        self.used[rel] = md5(p)
        return pd.read_csv(p, sep="\t", low_memory=False)


def same_pairs(base, d, rel, cols=("sens", "anty")):
    if len(d) != len(base) or not (d[KEY].to_numpy() == base[KEY].to_numpy()).all():
        sys.exit(f"STOP: {rel} does not hold the pairs of the step-04 table in the same order")
    for c in cols:
        if c in d.columns and not np.allclose(d[c].to_numpy(float), base[c].to_numpy(float), equal_nan=True):
            sys.exit(f"STOP: {rel}: column {c} differs from the step-04 table")


# column name, description, step, source (file pattern) ; filled in main()
DICT = []


def add(df, name, values, desc, step, src):
    df[name] = values
    DICT.append({"column": name, "description": desc, "step": step, "source": src})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tsv")
    a = ap.parse_args()
    S = Src(a.root)
    te_all = S.read("23_robustness/results/23c/23c_pair_te.tsv")
    sets = {k: S.read(f"23_robustness/results/23b/23b_{k}_matched_pairs.tsv") for k in ("all", "noTE", "noRepeat")}
    cg_all = S.read("28d_cpg/results/28d_pairs.tsv")
    parts = []
    for li, line in enumerate(LINES):
        rel = f"02_istart_baseline/results/arms_{line}.qc.tsv"
        b = S.read(rel)
        if (b["status"] != "ok").any() or b.duplicated(KEY).any() or set(b["line"]) != {line}:
            sys.exit(f"STOP: {rel}: unexpected status, duplicate pairs or line")
        o = pd.DataFrame()
        rec = li == 0  # write the dictionary once
        A = (lambda *x: add(o, *x)) if rec else (lambda n, v, *rest: o.__setitem__(n, v))
        A("pair_id", b["line"] + ":" + b["chrom"] + ":" + b["mid"].astype(str),
          "identifier of the pair: line:chromosome:midpoint", "this table", "")
        A("line", b["line"], "cell line of the PRO-cap data (K562 or HCT116)", "02-04", "arms_<line>.qc.tsv")
        A("chrom", b["chrom"], "chromosome (GRCh38)", "02-04", "arms_<line>.qc.tsv")
        A("mid", b["mid"], "midpoint of the pair, (plus-strand start + minus-strand start) // 2, as in the iSTART directionality table", "02-04",
          "arms_<line>.qc.tsv")
        A("pair_class", b["klasa"].map({"intron genu": "intronic", "międzygenowa": "intergenic"}),
          "intronic (in an intron of a host gene) or intergenic, as classified in the iSTART directionality table", "02-04", "arms_<line>.qc.tsv")
        if o["pair_class"].isna().any():
            sys.exit("STOP: unknown pair class")
        A("call_confidence", b["confidence"], "confidence of the PINTS divergent-peak call (iSTART)", "02-04",
          "arms_<line>.qc.tsv")
        A("pos_sense", b["pos_sense"], "position of the sense start (intronic: on the host strand; intergenic: "
          "plus-strand start, by convention)", "03", "arms_<line>.qc.tsv")
        A("pos_antisense", b["pos_antisense"], "position of the antisense start", "03", "arms_<line>.qc.tsv")
        A("strand_sense", b["strand_sense"], "strand of the sense start", "03", "arms_<line>.qc.tsv")
        A("arm_distance_nt", b["arm_distance"], "distance between the two starts (nt)", "03", "arms_<line>.qc.tsv")
        A("reads_sense", b["sens"], "PRO-cap reads at the sense start", "02-04", "arms_<line>.qc.tsv")
        A("reads_antisense", b["anty"], "PRO-cap reads at the antisense start", "02-04", "arms_<line>.qc.tsv")
        A("depth", b["sens"] + b["anty"], "reads_sense + reads_antisense", "this table", "")
        A("DI", b["DI"], "directionality index, log2((reads_sense + 1) / (reads_antisense + 1))", "02-04",
          "arms_<line>.qc.tsv")
        A("host_gene_id", b["host"], "host gene (GENCODE v47; intronic pairs)", "02-04", "arms_<line>.qc.tsv")
        A("host_gene_name", b["host_name"], "host gene name", "02-04", "arms_<line>.qc.tsv")
        A("host_strand", b["host_strand"], "strand of the host gene", "02-04", "arms_<line>.qc.tsv")
        A("host_activity_reads", b["host_signal"], "activity of the host gene: PRO-cap reads at its annotated starts "
          "(0 = silent, 1-4 = weakly active, >= 5 = active)", "02-04", "arms_<line>.qc.tsv")
        A("host_tx", b["host_tx"], "host transcription: mean Bru-seq signal on the host strand 5–1 kb upstream of the sense start, carried over from the iSTART directionality table "
          "(used only for the host-transcription thirds of step 02)", "02-04", "arms_<line>.qc.tsv")
        A("host_H3K36me3", b["k36"], "H3K36me3 peak within 1 kb of the pair (ENCODE replicated peaks)", "02-04",
          "arms_<line>.qc.tsv")
        for arm in ("sense", "antisense"):
            A(f"dist_tss_any_{arm}_nt", b[f"dist_tss_any_{arm}"], f"distance of the {arm} start to the nearest "
              "annotated start of any GENCODE v47 transcript (nt)", "04, 22c", "arms_<line>.qc.tsv")
            A(f"dist_tss_pc_{arm}_nt", b[f"dist_tss_pc_{arm}"], f"distance of the {arm} start to the nearest "
              "annotated start of a protein-coding transcript (nt)", "04, 22c", "arms_<line>.qc.tsv")
        A("dist_tss_any_min_nt", b["dist_tss_any_min"], "smaller of the two dist_tss_any distances (nt)", "04, 22c",
          "arms_<line>.qc.tsv")

        # Puffin-D (step 07)
        rel = f"07_puffin_pairs/results/puffin_pairs_{line}.tsv"
        d = S.read(rel)
        same_pairs(b, d, rel, ("sens", "anty", "DI"))
        A("gc_local", d["gc_local"], "GC fraction of the 1,000 nt centred on the sense start (used for the GC strata of the "
          "calibration)", "07", "puffin_pairs_<line>.tsv")
        A("puffin_status", d["pred_status"], "Puffin-D scoring status (ok = scored)", "07", "puffin_pairs_<line>.tsv")
        A("puffin_pred_sense", d["pred_sense"], "Puffin-D predicted initiation at the sense start", "07",
          "puffin_pairs_<line>.tsv")
        A("puffin_pred_antisense", d["pred_antisense"], "Puffin-D predicted initiation at the antisense start", "07",
          "puffin_pairs_<line>.tsv")
        A("puffin_pred_DI", d["pred_DI"], "Puffin-D predicted directionality (the predicted DI used for the shares, "
          "Fig. 2B)", "07", "puffin_pairs_<line>.tsv")

        # ProCapNet, K562 and Caco-2 models (steps 16, 21)
        for step, cell in (("16_procapnet", "K562"), ("21_caco2", "Caco2")):
            for kind in ("heldout", "trainfold"):
                rel = f"{step}/results/procapnet_{kind}_pairs_{line}.tsv"
                d = S.read(rel)
                same_pairs(b, d, rel, ("sens", "anty", "DI"))
                p = f"pc_{cell}_{kind}"
                what = ("the model that never saw the chromosome of the pair" if kind == "heldout"
                        else "a model trained on the chromosome of the pair")
                st = step.split("_")[0]
                A(f"{p}_status", d["pred_status"], f"ProCapNet ({cell} models), {what}: scoring status", st,
                  f"procapnet_{kind}_pairs_<line>.tsv")
                if kind == "heldout":
                    A(f"{p}_fold", d["pc_fold_heldout"], f"ProCapNet ({cell}): model (fold) that never saw the "
                      "chromosome", st, f"procapnet_{kind}_pairs_<line>.tsv")
                    A(f"{p}_pred_sense", d["pred_sense"], f"ProCapNet ({cell}), {what}: predicted counts at the "
                      "sense start", st, f"procapnet_{kind}_pairs_<line>.tsv")
                    A(f"{p}_pred_antisense", d["pred_antisense"], f"ProCapNet ({cell}), {what}: predicted counts "
                      "at the antisense start", st, f"procapnet_{kind}_pairs_<line>.tsv")
                A(f"{p}_pred_DI", d["pred_DI"], f"ProCapNet ({cell}), {what}: predicted DI", st,
                  f"procapnet_{kind}_pairs_<line>.tsv")

        # shuffles (step 20, K562 held-out models)
        for reg, txt in (("none", "original sequence"), ("window", "whole 2,114-nt window shuffled"),
                         ("element", "element shuffled"), ("flanks", "flanks shuffled, element intact")):
            rel = f"20_pair_shuffle/results/procapnet_shuf_{reg}_pairs_{line}.tsv"
            d = S.read(rel)
            same_pairs(b, d, rel, ("sens", "anty", "DI"))
            A(f"pc_K562_heldout_shuf_{reg}_pred_DI", d["pred_DI"], f"ProCapNet (K562 held-out models), {txt} "
              "(dinucleotide counts kept; mean over 5 shuffles): predicted DI (Fig. 2C-D)", "20",
              f"procapnet_shuf_{reg}_pairs_<line>.tsv")

        # transposable elements and repeats (step 23c; one table for both lines)
        te = te_all[te_all["line"] == line].reset_index(drop=True)
        same_pairs(b, te, "23c_pair_te.tsv", ("sens", "anty", "DI"))
        for c, desc in (("in_te", "pair overlaps a transposable element (RepeatMasker classes SINE, LINE, LTR, DNA, "
                         "Retroposon, RC)"),
                        ("in_any_repeat", "pair overlaps any RepeatMasker element"),
                        ("sense_start_in_te", "sense start lies in a transposable element"),
                        ("antisense_start_in_te", "antisense start lies in a transposable element"),
                        ("te_class", "class of the overlapping element"), ("te_family", "family of the element"),
                        ("te_name", "name of the element"), ("te_strand", "strand of the element"),
                        ("te_overlap_nt", "overlap with the element (nt)"),
                        ("te_orientation", "orientation of the element: same or opposite to the host strand (intronic pairs) "
                         "or to the sense start (intergenic pairs)"),
                        ("DI_te", "log2((reads at the start on the strand of the element + 1) / (reads at the other "
                         "start + 1))")):
            A(c, te[c], desc, "23c", "23c_pair_te.tsv")

        # depth-matched comparison (step 23b) and robustness sets
        m = sets["all"][sets["all"]["line"] == line]
        mk = m["chrom"] + ":" + m["mid"].astype(str)
        ck = m["ctrl_chrom"] + ":" + m["ctrl_mid"].astype(str)
        k_b = b["chrom"] + ":" + b["mid"].astype(str)
        if not mk.isin(k_b).all() or not ck.isin(k_b).all():
            sys.exit("STOP: matched pairs or their controls not found among the pairs")
        ctrl_of = dict(zip(mk, line + ":" + ck))
        used_ctrl = set(ck)
        A("matched_intronic", k_b.isin(set(mk)), "intronic pair of an active host (>= 5 reads) that entered the "
          "depth-matched comparison (Fig. 1B; 1,610 pairs in K562, 2,489 in HCT116)", "23b",
          "23b_all_matched_pairs.tsv")
        A("matched_control_pair_id", k_b.map(ctrl_of), "intergenic pair matched to it on read depth", "23b",
          "23b_all_matched_pairs.tsv")
        A("matched_control", k_b.isin(used_ctrl), "intergenic pair used as a matched control", "23b",
          "23b_all_matched_pairs.tsv")
        for k, txt in (("noTE", "without pairs overlapping a transposable element"),
                       ("noRepeat", "without pairs overlapping any repeat")):
            mm = sets[k][sets[k]["line"] == line]
            A(f"matched_intronic_{k}", k_b.isin(set(mm["chrom"] + ":" + mm["mid"].astype(str))),
              f"intronic pair in the depth-matched comparison {txt}", "23b", f"23b_{k}_matched_pairs.tsv")

        # position along the host gene and distance to CpG islands (steps 28b, 28d; matched intronic pairs)
        rel = f"28b_position/results/28b_{line}_pairs.tsv"
        p = S.read(rel)
        pk = p["chrom"] + ":" + p["mid"].astype(str)
        if set(pk) != set(mk):
            sys.exit(f"STOP: {rel} does not hold the matched pairs")
        pi = p.set_index(pk)
        for c, name, desc in (("gene_len_kb", "host_gene_length_kb", "length of the annotated host gene (kb)"),
                              ("rel", "position_in_gene", "position of the midpoint along the host gene, from 0 "
                               "(5' end) to 1 (3' end)"),
                              ("fifth", "gene_fifth", "fifth of the host gene (P1 = 5' fifth)"),
                              ("dist5_kb", "dist_from_5prime_end_kb", "distance from the annotated 5' end of the "
                               "host gene (kb)"),
                              ("dist_bin", "dist_from_5prime_bin", "distance bin (kb) used in Supplemental Fig. S5")):
            A(name, k_b.map(pi[c]), desc + "; matched intronic pairs", "28b", "28b_<line>_pairs.tsv")
        cg = cg_all[cg_all["line"] == line]
        ci = cg.set_index(cg["chrom"] + ":" + cg["mid"].astype(str))["cgi_dist_kb"]
        if set(ci.index) != set(mk):
            sys.exit("STOP: 28d_pairs.tsv does not hold the matched pairs")
        A("cgi_dist_kb", k_b.map(ci), "distance of the midpoint to the nearest CpG island (UCSC cpgIslandExt, kb); "
          "matched intronic pairs", "28d", "28d_pairs.tsv")
        parts.append(o)

    T = pd.concat(parts, ignore_index=True)
    # summary checks against the numbers in the text
    n = T.groupby(["line", "pair_class"]).size()
    print(n.to_string())
    print("matched intronic:", T.groupby("line")["matched_intronic"].sum().to_dict())
    print("intronic host genes:", T[T.pair_class == "intronic"].groupby("line")["host_gene_id"].nunique().to_dict())
    if T.groupby("line")["matched_intronic"].sum().to_dict() != {"HCT116": 2489, "K562": 1610}:
        sys.exit("STOP: matched counts differ from the text")
    D = pd.DataFrame(DICT)
    R = pd.DataFrame(
        [{"item": "content", "value": "Every iSTART start-site pair of K562 and HCT116 (sheet 'pairs', one row per "
          "pair) with the per-pair measures used in the study; sheet 'columns' describes each column and names the "
          "analysis step and result table it comes from."},
         {"item": "pairs", "value": "; ".join(f"{l} {c}: {v:,}" for (l, c), v in n.items())},
         {"item": "values", "value": "copied from the result tables without change; empty = not applicable "
          "(e.g. host columns of intergenic pairs) or not scored"},
         {"item": "built", "value": "build_supplemental_data.py (Supplemental Code), 12 Nov 2026"}]
        + [{"item": f"input {k}", "value": f"md5 {v}"} for k, v in sorted(S.used.items())])
    with pd.ExcelWriter(a.out, engine="openpyxl") as xw:
        R.to_excel(xw, sheet_name="README", index=False)
        D.to_excel(xw, sheet_name="columns", index=False)
        T.to_excel(xw, sheet_name="pairs", index=False)
        for ws, widths in ((xw.sheets["README"], {"A": 40, "B": 140}),
                           (xw.sheets["columns"], {"A": 34, "B": 110, "C": 10, "D": 38})):
            for col, w in widths.items():
                ws.column_dimensions[col].width = w
        xw.sheets["pairs"].freeze_panes = "B2"
    if a.tsv:
        T.to_csv(a.tsv, sep="\t", index=False)
    print(f"written {a.out}: {len(T):,} pairs x {T.shape[1]} columns")


if __name__ == "__main__":
    main()
