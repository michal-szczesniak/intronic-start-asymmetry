#!/usr/bin/env python3
"""Source Data for the figures, read directly from the result tables of the analysis steps.
Every CSV carries a 'source' column: the result file (path relative to the PURGE project folder) and its md5 checksum.

Inputs
  --data    folder holding the result tables of steps 02-23 at their project paths (the PURGE project folder,
            or the unpacked archive 23_robustness/results/23d_figure_tables.tar.gz)
  --meth    folder holding 26_methylation/results (steps 26b and 26c; the PURGE project folder)
  --old-sd  18_manuscript/source_data (two inputs that exist only as log lines, see below)
  --position  figures/figS_position/data: the strata tables of steps 28 and 29 (checked copies, see
            README_figS_position.txt there), for Figure 3E (HCT116 PRO-cap and K562 GRO-cap)
  --metagene  28g_metagene/results: strong sites in 1-kb bins from the 5' end of active genes (step 28g), Figure 3F

Two inputs exist only as log lines (no result table was written): the ProCapNet chromosome folds
(step 16b log) and the overlap of the variant classes of the two models (step 19c log). They are
copied from the earlier Source Data files built from those logs, and marked as such.

File names follow the figure numbers of the manuscript.
Supplemental Figs. S2 and S5 and Figure 5 have their own folders in figures/.

Usage: python3 make_source_data.py --data <PURGE> --meth <PURGE> --old-sd <PURGE>/18_manuscript/source_data \
           --position <PURGE>/figures/figS_position/data --metagene <PURGE>/28g_metagene/results --out source_data
"""
import argparse
import hashlib
import os

import numpy as np
import pandas as pd


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Reader:
    def __init__(self, root, label=""):
        self.root, self.label = root, label

    def __call__(self, rel):
        p = os.path.join(self.root, rel)
        return pd.read_csv(p, sep="\t"), f"{self.label}{rel} (md5 {md5(p)})"


def ratio_cols(df, d, lo, hi, prefix="ratio"):
    """deficit (%) = 100 (1 - S/A) -> ratio S/A with its CI (the CI of the deficit flips)."""
    df[prefix] = 1 - df[d] / 100
    df[prefix + "_lo95"] = 1 - df[hi] / 100
    df[prefix + "_hi95"] = 1 - df[lo] / 100
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--meth", required=True)
    ap.add_argument("--old-sd", required=True)
    ap.add_argument("--position", required=True)
    ap.add_argument("--metagene", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = Reader(a.data)
    os.makedirs(a.out, exist_ok=True)
    out = lambda df, name: df.to_csv(os.path.join(a.out, name), index=False, float_format="%.6g")

    # ---------------- Figure 1B: depth bands and matched comparison, host-gene CIs (step 23b)
    st, src = R("23_robustness/results/23b/23b_all_strata.tsv")
    mt, src2 = R("23_robustness/results/23b/23b_all_matched.tsv")
    f1 = st[["line", "stratum", "depth_lo", "depth_hi", "n_intron_active", "n_host_genes", "n_intergenic", "n_blocks",
             "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95"]].copy()
    f1["comparison"] = "depth band"
    f1["source"] = src
    m1 = mt[["line", "n_pairs", "n_host_genes", "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95",
             "wilcoxon_signed_rank_p"]].copy()
    m1["comparison"] = "matched"
    m1["source"] = src2
    out(pd.concat([f1, m1], ignore_index=True), "fig1_asymmetry.csv")

    # ---------------- Figure 2A: host activity (step 23b)
    ha, src = R("23_robustness/results/23b/23b_all_host_activity.tsv")
    f2a = ha[["line", "group", "n_pairs", "n_host_genes", "host_signal_min", "host_signal_median", "host_signal_max",
              "median_diff", "sense_weaker_pct", "pct_lo95", "pct_hi95", "wilcoxon_signed_rank_p"]].copy()
    f2a["source"] = src
    tr, src2 = R("23_robustness/results/23b/23b_all_trend.tsv")
    tr = tr.drop(columns=["ci", "tag"])
    tr["source"] = src2
    out(f2a, "fig2_host_activity.csv")
    out(tr, "fig2_host_trend.csv")

    # ---------------- Figure 2B-C: shares (step 23b) and composition (step 15)
    sh, src = R("23_robustness/results/23b/23b_shares.tsv")
    keep = ["puffin_D", "pc_K562_heldout", "pc_K562_trainfold", "pc_Caco2_heldout", "pc_Caco2_trainfold",
            "shuf_none", "shuf_window", "shuf_element", "shuf_flanks"]
    f2b = sh[sh["label"].isin(keep)][["label", "line", "group", "n_pairs", "n_host_genes", "n_intergenic", "n_blocks",
                                      "frac_linear", "frac_linear_lo95", "frac_linear_hi95", "frac_rank",
                                      "frac_rank_lo95", "frac_rank_hi95", "frac_depth", "frac_depth_lo95",
                                      "frac_depth_hi95"]].copy()
    f2b["source"] = src
    cm, src2 = R("15_composition/results/15_comp_models.tsv")
    c = cm[(cm["subset"] == "all") & (cm["window"].astype(str) == "300") & (cm["model"] == "mono")
           & cm["group"].isin(["intron_active", "intron_all"])]
    c2 = pd.DataFrame({"label": "composition_300nt_mono", "line": c["line"], "group": c["group"], "n_pairs": c["n_case"],
                       "n_intergenic": c["n_ctrl"], "frac_linear": c["share"], "frac_linear_lo95": c["share_lo95"],
                       "frac_linear_hi95": c["share_hi95"], "cv_r2_ctrl": c["cv_r2_ctrl"], "source": src2})
    out(pd.concat([f2b, c2], ignore_index=True), "fig2_shares.csv")

    # ---------------- Figure 3B: potential start sites, real and shuffled (step 13c, both models)
    rows = []
    for model, rel in (("Puffin-D", "13_genomewide/results/13c_shuffle_gw_summary.tsv"),
                       ("ProCapNet", "16_procapnet/results/13c_pc_summary.tsv")):
        d, src = R(rel)
        s = d[(d["group"] == "genome") & (d["seed"] == 1) & d["metric"].isin(["top0.001", "top0.0003", "top0.0001", "sum"])
              & d["set"].isin(["intron_protein_coding", "intergenic_plus_vs_minus"])].copy()
        s = ratio_cols(s, "real_deficit_pct", "real_lo95", "real_hi95", "real_ratio")
        s = ratio_cols(s, "shuf_deficit_pct", "shuf_lo95", "shuf_hi95", "shuf_ratio")
        s["model"] = model
        s["source"] = src
        rows.append(s[["model", "set", "seq", "metric", "n_units", "real_sense", "real_anti", "real_deficit_pct",
                       "real_lo95", "real_hi95", "real_ratio", "real_ratio_lo95", "real_ratio_hi95", "shuf_deficit_pct",
                       "shuf_ratio", "shuf_ratio_lo95", "shuf_ratio_hi95", "drop_pct_points", "drop_lo95", "drop_hi95",
                       "source"]])
    out(pd.concat(rows, ignore_index=True), "fig3_predicted_sites.csv")

    # ---------------- Figure 3C-E: measured sites and summed signal (step 17)
    ob, src = R("17_observed/results/17_obs_strata.tsv")
    src_ob = src
    h = ob[(ob["analysis"] == "host") & (ob["subset"] == "pc") & ob["bin"].isin(["silent", "low", "A1", "A2", "A3"])].copy()
    h = ratio_cols(h, "deficit", "lo", "hi")
    h["source"] = src
    out(h[["bin", "metric", "n_units", "sense", "anti", "deficit", "lo", "hi", "ratio", "ratio_lo95", "ratio_hi95",
           "median", "source"]].rename(columns={"median": "host_reads_median"}), "fig3_by_activity.csv")
    o = ob[(ob["analysis"] == "overall") & (ob["subset"] == "all")
           & ob["bin"].isin(["pc_all", "pc_active", "pc_silent", "lncRNA", "intergenic"]) & (ob["metric"] == "obs r10")].copy()
    o = ratio_cols(o, "deficit", "lo", "hi")
    o["source"] = src
    out(o[["bin", "metric", "n_units", "sense", "anti", "deficit", "lo", "hi", "ratio", "ratio_lo95", "ratio_hi95",
           "source"]], "fig3_strong_sites.csv")

    # ---------------- Figure 3E (replication): steps 28 (HCT116 PRO-cap) and 29 (K562 GRO-cap)
    P = Reader(a.position, "figures/figS_position/data/")
    rep = []
    for data_set, rel in (("HCT116 PRO-cap", "28_HCT116_strata.tsv"), ("K562 GRO-cap", "29_K562_GROcap_strata.tsv")):
        d, src = P(rel)  # K562 PRO-cap in Figure 3E comes from step 17 (fig3_strong_sites.csv)
        r = d[(d["metric"] == "strong") & d["subset"].isin(["active", "intergenic"]) & (d["stratum"] == "all")].copy()
        r = ratio_cols(r, "deficit", "lo", "hi")
        r["data_set"], r["source"] = data_set, src
        rep.append(r)
    cols = ["data_set", "subset", "stratum", "metric", "min_reads", "n_units", "sense", "anti", "deficit", "lo", "hi",
            "p_boot", "ratio", "ratio_lo95", "ratio_hi95", "source"]
    out(pd.concat(rep, ignore_index=True)[cols], "fig3_replication.csv")
    MG = Reader(a.metagene, "28g_metagene/results/")
    mg = []
    for rel in ("28g_K562_bins.tsv", "28g_HCT116_bins.tsv", "28g_GROcap_bins.tsv"):
        d, src = MG(rel)
        d["source"] = src
        mg.append(d)
    out(pd.concat(mg, ignore_index=True), "fig3_metagene.csv")

    # ---------------- Figure 4: population test (12c Puffin-D, 19b ProCapNet)
    cls, dif = [], []
    for model, folder in (("Puffin-D", "12_selection/results/gw"), ("ProCapNet", "19_procapnet_gnomad/results/gw_procapnet")):
        for setting in ("base", "noCpG", "gc5", "noCpG_gc5"):
            d, src = R(f"{folder}/{setting}_classes.tsv")
            d["model"], d["setting"], d["source"] = model, setting, src
            cls.append(d)
            p, src2 = R(f"{folder}/{setting}_pairs.tsv")
            p["model"], p["setting"], p["source"] = model, setting, src2
            dif.append(p)
    out(pd.concat(cls, ignore_index=True), "fig4_classes.csv")
    out(pd.concat(dif, ignore_index=True), "fig4_differences.csv")

    # ---------------- Supplementary figures
    nb, src = R("17_observed/results/17b_neighbors_strata.tsv")
    s2 = nb[(nb["analysis"] == "rt20k_pr5k") & nb["subset"].isin(["obs", "models"])
            & nb["bin"].isin(["isolated", "antisense_only", "sense_only", "both"])
            & nb["metric"].isin(["obs r1", "obs r3", "obs r10", "Puffin-D sum", "ProCapNet sum"])].copy()
    s2 = ratio_cols(s2, "deficit", "lo", "hi")
    s2["source"] = src
    out(s2, "figS3_neighbours.csv")

    rows = []
    for model, rel in (("Puffin-D", "13_genomewide/results/13d_gc_strata.tsv"), ("ProCapNet", "16_procapnet/results/13d_pc_strata.tsv")):
        d, src = R(rel)
        s = d[(d["analysis"] == "pc_by_keto") & d["metric"].isin(["top0.001", "top0.0001"])].copy()
        s = ratio_cols(s, "real", "real_lo", "real_hi")
        s["model"], s["source"] = model, src
        rows.append(s[["model", "bin", "metric", "n_units", "sites", "keto_median", "real", "real_lo", "real_hi",
                       "ratio", "ratio_lo95", "ratio_hi95", "source"]])
    k = ob[(ob["analysis"] == "keto") & (ob["subset"] == "all") & (ob["metric"] == "obs r10")].copy()
    k = ratio_cols(k, "deficit", "lo", "hi")
    k["model"], k["source"] = "measured", src_ob
    k = k.rename(columns={"deficit": "real", "lo": "real_lo", "hi": "real_hi", "median": "keto_median"})
    k["metric"] = "obs r10"
    rows.append(k[["model", "bin", "metric", "n_units", "keto_median", "real", "real_lo", "real_hi", "ratio",
                   "ratio_lo95", "ratio_hi95", "source"]])
    s3 = pd.concat(rows, ignore_index=True)
    dif_rows = s3["bin"].astype(str).str.startswith("Q5-Q1")
    s3.loc[dif_rows, ["ratio", "ratio_lo95", "ratio_hi95"]] = np.nan  # a difference of deficits has no strand ratio
    out(s3, "figS6_keto.csv")

    fo = pd.read_csv(os.path.join(a.old_sd, "fig2_folds.csv"))
    fo["source"] = fo["source"] + " [log-derived; no result table]"
    out(fo, "figS1_folds.csv")
    ov = pd.read_csv(os.path.join(a.old_sd, "fig5_class_overlap.csv"))
    ov["source"] = ov["source"] + " [log-derived; no result table]"
    out(ov, "figS8_class_overlap.csv")

    # ---------------- Supplemental Fig. S4: gene-body DNA methylation (steps 26b and 26c, K562)
    M = Reader(a.meth)
    g, src = M("26_methylation/results/26b_meth_primary_genes.tsv")
    rows = []
    for grp, lab in (("silent", "silent"), ("low", "weakly active"), ("A1", "active, least active third"),
                     ("A2", "active, middle third"), ("A3", "active, most active third")):
        m = g[(g["gene_type"] == "protein_coding") & (g["activity"] == grp) & (g["meth_cpg"] >= 20)]["meth"].dropna()
        p10, q1, med, q3, p90 = np.percentile(m.to_numpy(), [10, 25, 50, 75, 90])
        rows.append({"group": grp, "label": lab, "n_genes": len(m), "p10": p10, "q1": q1, "median": med, "q3": q3,
                     "p90": p90, "source": src})
    out(pd.DataFrame(rows), "figS4_meth_by_activity.csv")
    stt, src = M("26_methylation/results/26b_meth_primary_strata.tsv")
    s4 = stt[(stt["metric"] == "obs r10") & stt["subset"].isin(["A1", "A2", "A3", "pc_active", "intergenic"])
             & stt["bin"].isin(["M1", "M2", "M3"])].copy()
    s4 = ratio_cols(s4, "deficit", "lo", "hi")
    s4["source"] = src
    out(s4[["subset", "bin", "metric", "n_units", "sense", "anti", "deficit", "lo", "hi", "ratio", "ratio_lo95",
            "ratio_hi95", "median_meth", "source"]], "figS4_strand_ratio.csv")
    rows = []
    for proc, lab in (("primary", "ENCODE default processing"), ("alternative", "earlier processing")):
        d, src = M(f"26_methylation/results/26b_meth_{proc}_strata.tsv")
        r = d[(d["subset"] == "pc_active") & (d["bin"] == "M3-M1") & (d["metric"] == "obs r10")].iloc[0]
        rows.append({"processing": proc, "label": lab, "comparison": "pooled over activity thirds (pre-specified)",
                     "contrast": "M3-M1", "metric": "obs r10", "value": r["deficit"], "lo95": r["lo"], "hi95": r["hi"],
                     "p": r["p_boot"], "source": src})
        d, src = M(f"26_methylation/results/26c_posthoc_{proc}_contrasts.tsv")
        r = d[(d["contrast"] == "M3-M1") & (d["metric"] == "obs r10")].iloc[0]
        rows.append({"processing": proc, "label": lab, "comparison": "activity thirds weighted equally (post hoc)",
                     "contrast": "M3-M1", "metric": "obs r10", "value": r["equal_weight"], "lo95": r["equal_lo95"],
                     "hi95": r["equal_hi95"], "p": r["p_equal"], "source": src})
    out(pd.DataFrame(rows), "figS4_contrasts.csv")
    print("written:", sorted(os.listdir(a.out)))


if __name__ == "__main__":
    main()
