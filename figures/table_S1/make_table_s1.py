#!/usr/bin/env python3
"""Supplemental Table S1 (11 Nov 2026): robustness checks, the DNA methylation
analysis and the population test, read directly from the result tables of the analysis steps. Every row
carries the result file it comes from (path relative to the PURGE project folder) and that file's md5
checksum. Values are copied, not recomputed; display formats round them, full precision is kept.

Inputs
  --data  the archive 23_robustness/results/23d_figure_tables.tar.gz (md5 fb044b571501ff1deb44a5f00c7097dc),
          unpacked; sheets S1A-S1G and S1I
  --meth  a folder holding 26_methylation/results/ with the step-26b strata tables (26b_meth_primary_strata.tsv,
          26b_meth_alternative_strata.tsv) and the step-26c tables written by 26c_posthoc_rebuild.py
          (26c_posthoc_<set>_contrasts.tsv, _cells.tsv, _rho.tsv for set = primary, alternative); sheet S1H

History. An earlier version of this script (md5 e4f2b94bef4dd621e195a64ade82e5f3, 23 Oct 2026) wrote sheets
S1A-S1H, S1H being the population test; with the identical 23d archive it reproduces that workbook cell for cell.
This version adds the methylation sheet as S1H (the manuscript cites "Supplemental Table S1 (sheet S1H)" for it),
moves the population test to S1I and uses the figure numbering of the manuscript in titles and notes. Sheets
S1A-S1G and S1I hold the same values as the earlier version.

Usage: python3 make_table_s1.py --data <unpacked 23d archive> --meth <folder with 26_methylation/results>
           --out Supplemental_Table_S1.xlsx
"""
import argparse
import hashlib
import os

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Reader:
    """Reads a result table by its path relative to the project folder, from the first root that holds it."""

    def __init__(self, *roots):
        self.roots = roots

    def __call__(self, rel):
        hits = [os.path.join(r, rel) for r in self.roots if os.path.isfile(os.path.join(r, rel))]
        if len(hits) != 1:
            raise SystemExit(f"STOP: {rel} found in {len(hits)} of the input folders (expected exactly 1)")
        p = hits[0]
        d = pd.read_csv(p, sep="\t")
        d["source_file"] = rel
        d["source_md5"] = md5(p)
        return d


CI_HOST = "95% CI: 10,000 resamples of host genes (intergenic pairs within depth bands: 1-Mb blocks)"


def sheets(R):
    out = {}

    # ---- S1A: the asymmetry, by depth band and after matching; two ways of resampling
    st = R("23_robustness/results/23b/23b_all_strata.tsv")
    a = st[["line", "stratum", "depth_lo", "depth_hi", "n_intron_active", "n_host_genes", "n_intergenic", "n_blocks",
            "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95", "source_file", "source_md5"]].copy()
    a.insert(1, "comparison", "depth band (fifths of S + A)")
    a["stratum"] = a["stratum"] + 1
    a["resampling"] = CI_HOST
    mt = R("23_robustness/results/23b/23b_all_matched.tsv")
    b = mt[["line", "n_pairs", "n_host_genes", "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95",
            "wilcoxon_signed_rank_p", "source_file", "source_md5"]].copy()
    b.insert(1, "comparison", "1:1 depth-matched")
    b["resampling"] = "95% CI: 10,000 resamples of host genes; P: two-sided Wilcoxon signed-rank test on matched pairs"
    c = R("02_istart_baseline/results/istart_depth_matched.tsv")
    c = c.rename(columns={"n_matched": "n_pairs", "sense_weaker_pct_lo95": "pct_lo95",
                          "sense_weaker_pct_hi95": "pct_hi95"})
    c = c[["line", "n_pairs", "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95", "mannwhitney_p",
           "source_file", "source_md5"]]
    c.insert(1, "comparison", "1:1 depth-matched")
    c["resampling"] = "95% CI: 10,000 resamples of single pairs; P: two-sided Mann-Whitney U test"
    out["S1A Asymmetry"] = ("Sense start of intronic pairs in active hosts weaker than in intergenic pairs: by depth band "
                            "and after 1:1 depth matching, with confidence intervals from resampling host genes or "
                            "single pairs. sense_weaker_pct = 100 x (1 - 2^median_diff).",
                            [("Depth bands, host-gene resampling (Figure 1B)", a.drop(columns=["comparison"])),
                             ("1:1 depth-matched, host-gene resampling (Figure 1B)", b.drop(columns=["comparison"])),
                             ("1:1 depth-matched, resampling single pairs (for comparison)", c.drop(columns=["comparison"]))])

    # ---- S1B: pairs near annotated transcript starts removed
    rows = []
    for m in (0, 100, 1000):
        d = R(f"22_open_items/results/22f_tables/depth_matched.min{m}.tsv")
        d.insert(1, "removed_within_nt_of_any_TSS", m if m else "none")
        rows.append(d)
    f = pd.concat(rows, ignore_index=True)
    cnt = R("22_open_items/results/22f_tables/22f_filter_counts.tsv")
    ig = cnt[cnt["group"] == "intergenic"][["line", "min_dist_nt", "n_kept", "pct_removed"]].rename(
        columns={"n_kept": "n_intergenic_kept", "pct_removed": "pct_intergenic_removed"})
    ia = cnt[cnt["group"] == "intron_active"][["line", "min_dist_nt", "pct_removed"]].rename(
        columns={"pct_removed": "pct_intronic_removed"})
    f["min_dist_nt"] = f["removed_within_nt_of_any_TSS"].replace("none", 0).astype(int)
    f = f.merge(ia, on=["line", "min_dist_nt"], how="left").merge(ig, on=["line", "min_dist_nt"], how="left")
    f = f.drop(columns=["min_dist_nt"])
    f["resampling"] = "95% CI: 10,000 resamples of single pairs (step 02c, unchanged); P: two-sided Mann-Whitney U test"
    td = R("22_open_items/results/22c_tss_distance.tsv")
    out["S1B Promoter filter"] = ("Distance of the pairs from annotated transcript starts, and the depth-matched comparison "
                                  "after removing, from both groups, every pair with either start within 100 nt or 1 kb of the "
                                  "start of any GENCODE v47 transcript. With nothing removed the run reproduced the original "
                                  "tables byte for byte. Distances: dist_tss_any_min, the start of the pair closer to any "
                                  "annotated transcript start; dist_tss_pc_min, the same for protein-coding transcripts; "
                                  "_sense / _antisense, each start separately. Supplemental Results.",
                                  [("Distance from annotated transcript starts (nt)", td),
                                   ("Depth-matched comparison after removing pairs near annotated transcript starts", f)])

    # ---- S1C: host activity, trend and H3K36me3
    ha = R("23_robustness/results/23b/23b_all_host_activity.tsv").drop(columns=["tag", "ci"])
    ha.insert(1, "analysis", "host activity group")
    tr = R("23_robustness/results/23b/23b_all_trend.tsv").drop(columns=["tag", "ci"])
    tr.insert(1, "analysis", "trend with host activity")
    k = R("23_robustness/results/23b/23b_all_k36.tsv").drop(columns=["tag", "ci"])
    k.insert(1, "analysis", "H3K36me3 within 1 kb")
    out["S1C Host activity, H3K36me3"] = (
        "Intronic pairs matched to intergenic pairs of the same depth, by activity of the host gene (PRO-cap reads at its "
        "annotated starts: silent 0, weak 1-4, thirds of active genes), the Spearman correlation between host activity "
        "and the within-pair difference, and pairs with or without an H3K36me3 peak within 1 kb (active hosts). "
        + CI_HOST + ".", [("By activity of the host gene (Figure 2A)", ha.drop(columns=["analysis"])),
                          ("Trend with host activity (all intronic pairs)", tr.drop(columns=["analysis"])),
                          ("H3K36me3 peak within 1 kb of the pair (active hosts)", k.drop(columns=["analysis"]))])

    # ---- S1D: share of the asymmetry predicted from sequence; composition
    sh = R("23_robustness/results/23b/23b_shares.tsv")
    lab = {"puffin_D": "Puffin-D", "pc_K562_heldout": "ProCapNet K562, held-out",
           "pc_K562_trainfold": "ProCapNet K562, train-fold", "pc_Caco2_heldout": "ProCapNet Caco-2, held-out",
           "pc_Caco2_trainfold": "ProCapNet Caco-2, train-fold", "shuf_none": "K562 held-out, original sequence (shuffle run)",
           "shuf_window": "K562 held-out, whole 2,114-nt window shuffled", "shuf_element": "K562 held-out, element shuffled",
           "shuf_flanks": "K562 held-out, flanks shuffled"}
    sh = sh[sh["label"].isin(lab)].copy()
    sh.insert(0, "predictor", sh["label"].map(lab))
    keep = ["predictor", "line", "group", "n_pairs", "n_host_genes", "n_intergenic", "n_blocks", "frac_linear",
            "frac_linear_lo95", "frac_linear_hi95", "frac_rank", "frac_rank_lo95", "frac_rank_hi95", "frac_depth",
            "frac_depth_lo95", "frac_depth_hi95", "source_file", "source_md5"]
    sh = sh[keep]
    cm = R("15_composition/results/15_comp_models.tsv")
    cm = cm[(cm["subset"] == "all") & cm["group"].isin(["intron_active", "intron_all"])].copy()
    names = {"pred": "Puffin-D alone", "mono": "composition: T/A and G/C skews", "mono+di": "composition: + dinucleotide skews",
             "pred+mono": "Puffin-D + T/A and G/C skews", "pred+mono+di": "Puffin-D + all skews"}
    cm.insert(0, "predictor", cm["model"].map(names) + cm["window"].map(lambda w: "" if w == "-" else f" (+-{w} nt)"))
    cm = cm.rename(columns={"n_case": "n_pairs", "n_ctrl": "n_intergenic", "share": "frac_linear",
                            "share_lo95": "frac_linear_lo95", "share_hi95": "frac_linear_hi95"})
    cm = cm[["predictor", "line", "group", "n_pairs", "n_intergenic", "cv_r2_ctrl", "frac_linear", "frac_linear_lo95",
             "frac_linear_hi95", "source_file", "source_md5"]]
    out["S1D Share from sequence"] = (
        "Share of the measured asymmetry of intronic pairs predicted from sequence. frac_linear: calibration by least "
        "squares through the origin on intergenic pairs taken in both orientations (Figure 2B); frac_rank: (AUC - 0.5) "
        "predicted / (AUC - 0.5) measured; frac_depth: calibration within thirds of read depth. Held-out: the ProCapNet "
        "model trained without the pair's chromosome; train-fold: mean of the six models that saw it. 95% CI: 2,000 "
        "resamples of host genes, intergenic pairs in 1-Mb blocks. Composition rows (step 15): linear models fitted on "
        "intergenic pairs, 1,000 resamples of host genes and of intergenic pairs; cv_r2_ctrl = five-fold cross-validated "
        "R2 among intergenic pairs.", [("Sequence models and shuffles (Figure 2B, 2C; Supplemental Fig. S1B)", sh),
                                         ("Base composition, alone and added to Puffin-D (step 15)", cm)])

    # ---- S1E: transposable elements
    ov = R("23_robustness/results/23c/23c_overlap.tsv")
    ov.insert(1, "table", "overlap")
    orr = R("23_robustness/results/23c/23c_orientation.tsv")
    orr.insert(1, "table", "orientation relative to the host gene (intronic pairs in active hosts)")
    td = R("23_robustness/results/23c/23c_te_direction.tsv")
    td.insert(1, "table", "direction imposed by the TE: DI_TE = log2((reads at the start on the TE strand + 1)/(reads at the other start + 1))")
    dbt = R("23_robustness/results/23c/23c_di_by_te.tsv")
    dbt.insert(1, "table", "measured DI of intronic pairs by TE orientation")
    rows = []
    for tag in ("noTE", "noRepeat"):
        d = R(f"23_robustness/results/23b/23b_{tag}_matched.tsv").drop(columns=["ci"])
        d.insert(1, "table", "depth-matched comparison without pairs overlapping a TE" if tag == "noTE"
                 else "depth-matched comparison without pairs overlapping any RepeatMasker element")
        rows.append(d)
    fc = R("23_robustness/results/23c/23c_filter_counts.tsv")
    fc.insert(1, "table", "pairs removed by each filter")
    out["S1E Transposable elements"] = (
        "Transposable elements (RepeatMasker classes SINE, LINE, LTR, DNA, Retroposon, RC). A pair overlaps a TE when its "
        "element, extended by 50 nt beyond each start, overlaps one; the TE with the largest overlap gives class and "
        "orientation. Background orientation: all TEs with their midpoint in the bodies of the same host genes. "
        "Depth-matched rows: " + CI_HOST + ".",
        [("Overlap of the elements with TEs (% of pairs)", ov.drop(columns=["table"])),
         ("Orientation of the TE relative to the host gene: intronic pairs vs TEs in the same gene bodies (background)",
          orr.drop(columns=["table"])),
         ("Direction imposed by the TE: DI_TE = log2((reads at the start on the TE strand + 1)/(reads at the other start + 1))",
          td.drop(columns=["table"])),
         ("Measured DI of intronic pairs by TE orientation", dbt.drop(columns=["table"])),
         ("Depth-matched comparison without pairs overlapping a TE", rows[0].drop(columns=["table"])),
         ("Depth-matched comparison without pairs overlapping any RepeatMasker element", rows[1].drop(columns=["table"])),
         ("Pairs removed by each filter", fc.drop(columns=["table"]))])

    # ---- S1F: genome-wide census of potential start sites
    rows = []
    for model, rel in (("Puffin-D", "13_genomewide/results/13c_shuffle_gw_summary.tsv"),
                       ("ProCapNet", "16_procapnet/results/13c_pc_summary.tsv")):
        d = R(rel)
        d = d[(d["group"].isin(["genome", "heldout", "training"])) & (d["seed"] == 1)].copy()
        d.insert(0, "model", model)
        rows.append(d.drop(columns=["seed"]))
    out["S1F Genome-wide census"] = (
        "Potential start sites (highest 0.1%, 0.03% and 0.01% of predictions: top0.001, top0.0003, top0.0001) and summed "
        "signal (sum) on the gene strand and the opposite strand of introns (intergenic DNA: plus and minus strand), in "
        "real sequence and after shuffling with base and dinucleotide counts kept. Deficit = 100 x (1 - gene strand / "
        "opposite strand); drop_pct_points = real minus shuffled deficit. seq all: all eligible positions; seq unique: "
        "positions not soft-masked as repeats. group genome: all autosomes; heldout / training: chromosomes the model was "
        "or was not trained on. 95% CI: 2,000 resamples of genes (intergenic: 1-Mb blocks).",
        [("", pd.concat(rows, ignore_index=True))])

    # ---- S1G: silent genes and their neighbours
    nb = R("17_observed/results/17b_neighbors_strata.tsv")
    names = {"rt20k_pr5k": "own strand to 20 kb past the end, opposite strand 5 kb upstream (Supplemental Fig. S3)",
             "rt50k_pr10k": "own strand to 50 kb past the end, opposite strand 10 kb upstream",
             "strict_rt50k_pr10k": "as rt50k_pr10k, neighbours counted as active from one read at their starts",
             "iso_keto": "silent genes without active neighbours, by keto skew"}
    nb.insert(1, "coverage_rule", nb["analysis"].map(names))
    out["S1G Neighbours of silent genes"] = (
        "Silent protein-coding genes classed by the strand on which transcription from active neighbouring genes reaches "
        "them, under three coverage rules. Deficit = 100 x (1 - gene strand / opposite strand) for measured sites with at "
        "least 1, 3 or 10 reads and for summed predicted signal. 95% CI: 2,000 resamples of genes. Supplemental "
        "Fig. S3; Supplemental Results.", [("", nb)])

    # ---- S1H: DNA methylation (step 26b; post hoc weighting, step 26c)
    mb = []
    for st, lab in (("primary", "ENCODE default processing (ENCFF660IHA, ENCFF328NMN)"),
                    ("alternative", "earlier ENCODE processing of the same reads (ENCFF867JRG, ENCFF721JMB)")):
        d = R(f"26_methylation/results/26b_meth_{st}_strata.tsv")
        d.insert(0, "processing", lab)
        mb.append(d)
    pc, cl, rh = [], [], []
    for st in ("primary", "alternative"):
        pc.append(R(f"26_methylation/results/26c_posthoc_{st}_contrasts.tsv"))
        cl.append(R(f"26_methylation/results/26c_posthoc_{st}_cells.tsv"))
        rh.append(R(f"26_methylation/results/26c_posthoc_{st}_rho.tsv"))
    out["S1H DNA methylation"] = (
        "Gene-body DNA methylation of K562 (ENCODE whole-genome bisulfite sequencing, ENCSR765JPC) and strong measured "
        "start sites. Methylation of a gene: mean over the CpGs with at least 10 reads at its eligible intronic positions; "
        "genes with at least 20 such CpGs. Methylation thirds M1 (low) to M3 (high) were formed within each group of host "
        "activity (silent, low = weakly active, A1-A3 = thirds of active genes). Deficit = 100 x (1 - gene strand / "
        "opposite strand), summed over the genes of a group; metrics: measured sites with at least 1, 3 or 10 reads "
        "(obs r1, r3, r10) and summed reads (obs sum), and weak, strong and summed predicted signal of Puffin-D and "
        "ProCapNet. Planned comparison (specified before the data were seen): subset pc_active, metric obs r10, "
        "contrast M3-M1, i.e. the methylation thirds pooled over A1-A3. Equal-weight comparison (chosen after seeing "
        "the data): mean of the within-third contrasts. Negative control: intergenic 1-Mb blocks in thirds of their "
        "methylation, plus strand divided by minus strand. 95% CI: 2,000 resamples of genes within each combination of "
        "activity and methylation third (intergenic: blocks); p: two-sided bootstrap. Supplemental Fig. S4; "
        "Supplemental Results.",
        [("Step 26b: deficits by methylation third, every activity group and metric", pd.concat(mb, ignore_index=True)),
         ("Step 26c: pooled (planned) and equal-weight (post hoc) contrasts", pd.concat(pc, ignore_index=True)),
         ("Step 26c: cells of activity third x methylation third, strong sites (>= 10 reads); weight_in_pooled = share "
          "of the opposite-strand strong sites of the methylation third contributed by the cell (Supplemental Fig. S4C)",
          pd.concat(cl, ignore_index=True)),
         ("Step 26c: Spearman correlation of methylation with host activity and with eligible length, inside each "
          "activity third", pd.concat(rh, ignore_index=True))])

    # ---- S1I: population test, every class, set and setting
    cls, dif = [], []
    for model, folder in (("Puffin-D", "12_selection/results/gw"), ("ProCapNet", "19_procapnet_gnomad/results/gw_procapnet")):
        for setting in ("base", "noCpG", "gc5", "noCpG_gc5"):
            d = R(f"{folder}/{setting}_pairs.tsv")
            d.insert(0, "model", model)
            d.insert(1, "expectation", setting)
            dif.append(d)
            c = R(f"{folder}/{setting}_classes.tsv")
            c.insert(0, "model", model)
            c.insert(1, "expectation", setting)
            cls.append(c)
    out["S1I Population test"] = (
        "Excess of singletons for every class of variants (sets: intronic variants of protein-coding genes, all; lncRNA "
        "intronic, 30% subset; intergenic, 20% subset) and gene-strand minus opposite-strand differences (intergenic: "
        "plus minus minus strand), with classes from Puffin-D and from ProCapNet. Expectation: base, same mutation type "
        "and three-base context; noCpG, CpG sites excluded; gc5, contexts split into five bins of local GC content; "
        "noCpG_gc5, both. Class names: creates/destroys_<strand>_<threshold> (site created or removed at the weak or "
        "medium threshold); gain_<strand>_strong / _mid (strong or moderate increase); loss_<strand>_strong (strong "
        "decrease); _unique, variants outside repeats. Weak and medium thresholds: Puffin-D 0.0439 and 0.1432, "
        "ProCapNet 1.131 and 2.658. 95% CI: 1,000 resamples of genes (intergenic: 1-Mb blocks). Figure 4B, C; "
        "Supplemental Fig. S7.",
        [("Differences between strands", pd.concat(dif, ignore_index=True)),
         ("Excess of singletons per class", pd.concat(cls, ignore_index=True))])
    return out


def write(out, path):
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    F = "Arial"
    ws["A1"] = "Supplemental Table S1. Robustness checks, DNA methylation and the population test"
    ws["A1"].font = Font(name=F, bold=True, size=12)
    ws["A3"] = ("Each sheet is copied from the result tables of the analysis (no values recomputed). Columns source_file "
                "and source_md5 name the result file (relative to the project folder) and its md5 checksum. Built by "
                "make_table_s1.py.")
    ws["A3"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 32
    ws.column_dimensions["B"].width = 110
    ws.merge_cells("A3:B3")
    ws.row_dimensions[3].height = 45
    r = 5
    for name, (desc, _) in out.items():
        ws.cell(r, 1, name).font = Font(name=F, bold=True)
        c = ws.cell(r, 2, desc)
        c.font = Font(name=F)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = max(30, 15 * (len(desc) // 110 + 1))
        r += 1
    for row in ws.iter_rows():
        for c in row:
            if c.font.name != F:
                c.font = Font(name=F, bold=c.font.bold, size=c.font.size)
    head = PatternFill("solid", fgColor="F0F0F0")
    for name, (desc, blocks) in out.items():
        s = wb.create_sheet(name[:31])
        s["A1"] = name
        s["A1"].font = Font(name=F, bold=True, size=11)
        s["A2"] = desc
        s["A2"].font = Font(name=F, italic=True)
        s["A2"].alignment = Alignment(wrap_text=True, vertical="top")
        s.merge_cells(start_row=2, start_column=1, end_row=2, end_column=12)
        s.row_dimensions[2].height = max(30, 15 * (len(desc) // 150 + 1))
        widths = {}
        r = 4
        for title, df in blocks:
            if title:
                s.cell(r, 1, title).font = Font(name=F, bold=True)
                r += 1
            for j, col in enumerate(df.columns, start=1):
                c = s.cell(r, j, col)
                c.font = Font(name=F, bold=True)
                c.fill = head
                widths[j] = max(widths.get(j, 10), min(40, len(str(col)) + 2))
            r += 1
            for rec in df.itertuples(index=False):
                for j, v in enumerate(rec, start=1):
                    if isinstance(v, float) and pd.isna(v):
                        v = None
                    elif hasattr(v, "item"):
                        v = v.item()
                    c = s.cell(r, j, v)
                    c.font = Font(name=F)
                    if isinstance(v, float):
                        if v != 0 and abs(v) < 0.001:
                            c.number_format = "0.00E+00"
                        else:
                            c.number_format = "0.0000" if abs(v) < 10 else "0.00"
                r += 1
            r += 1
        for j, w in widths.items():
            s.column_dimensions[get_column_letter(j)].width = w
    wb.save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--meth", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = sheets(Reader(a.data, a.meth))
    write(out, a.out)
    for k, (_, blocks) in out.items():
        print(k, "|", "; ".join(f"{t[:40] or 'table'}: {d.shape[0]}x{d.shape[1]}" for t, d in blocks))


if __name__ == "__main__":
    main()
