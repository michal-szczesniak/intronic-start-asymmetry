#!/usr/bin/env python3
"""Checks every number that the manuscript (main text, Supplemental Results and Methods, legends) takes from
Supplemental Table S1 against the cells of the workbook itself, at the precision printed in the text.

Usage: python3 check_table_s1.py Supplemental_Table_S1.xlsx
"""
import sys

import openpyxl
import pandas as pd


def blocks(ws):
    """Tables of one sheet: a header row is a row whose first cell is bold and filled; rows follow until a blank row."""
    out, rows = [], list(ws.iter_rows())
    i = 3
    while i < len(rows):
        c = rows[i][0]
        if c.value is not None and c.font.bold and c.fill.fgColor.rgb not in (None, "00000000"):
            hdr = [x.value for x in rows[i] if x.value is not None]
            data = []
            j = i + 1
            while j < len(rows) and any(x.value is not None for x in rows[j]):
                data.append([x.value for x in rows[j]][: len(hdr)])
                j += 1
            out.append(pd.DataFrame(data, columns=hdr))
            i = j
        i += 1
    return out


wb = openpyxl.load_workbook(sys.argv[1])
B = {ws.title: blocks(ws) for ws in wb.worksheets[1:]}
ok, bad = 0, []


def check(label, value, expected, digits):
    global ok
    if value is None or round(float(value), digits) != round(expected, digits):
        bad.append(f"{label}: table {value!r}, text {expected}")
    else:
        ok += 1


def one(df, **kw):
    m = df
    for k, v in kw.items():
        m = m[m[k] == v]
    if len(m) != 1:
        raise SystemExit(f"STOP: {kw} matches {len(m)} rows")
    return m.iloc[0]


# S1A: matched comparison, host genes vs single pairs (main text: "almost the same intervals")
host, single = B["S1A Asymmetry"][1], B["S1A Asymmetry"][2]
for ln, pct, lo, hi in (("K562", 33, 13, 46), ("HCT116", 21, 12, 29)):
    r = one(host, line=ln)
    check(f"S1A host {ln} pct", r["sense_weaker_pct"], pct, 0)
    check(f"S1A host {ln} lo", r["pct_lo95"], lo, 0)
    check(f"S1A host {ln} hi", r["pct_hi95"], hi, 0)
    s = one(single, line=ln)
    check(f"S1A single {ln} pct", s["sense_weaker_pct"], pct, 0)
# S1B: distances and the promoter filter (Supplemental Results)
td, pf = B["S1B Promoter filter"]
for ln, pct_all, pct_act, med_kb, ig in (("K562", 7.9, 6.5, 4, 37), ("HCT116", 6.0, 5.3, 5, 28)):
    r = one(td, line=ln, group="intron_all", measure="dist_tss_any_min")
    check(f"S1B {ln} % within 100 nt, all intronic (6-8%)", r["pct_within_100"], pct_all, 1)
    check(f"S1B {ln} median distance kb (4-5 kb)", r["median_nt"] / 1000, med_kb, 0)
    r = one(td, line=ln, group="intron_active", measure="dist_tss_any_min")
    check(f"S1B {ln} % within 100 nt, active hosts (5.3-6.5%)", r["pct_within_100"], pct_act, 1)
    r = one(td, line=ln, group="intergenic", measure="dist_tss_any_min")
    check(f"S1B {ln} intergenic % within 100 nt (37%, 28%)", r["pct_within_100"], ig, 0)
for ln, d, pct in (("K562", 100, 35), ("HCT116", 100, 22), ("K562", 1000, 32), ("HCT116", 1000, 19)):
    r = one(pf, line=ln, removed_within_nt_of_any_TSS=d)
    check(f"S1B {ln} removed within {d} nt", r["sense_weaker_pct"], pct, 0)
# S1C: silent hosts, H3K36me3
ha, tr, k = B["S1C Host activity, H3K36me3"]
for ln, pct in (("K562", 47), ("HCT116", 24)):
    check(f"S1C silent {ln}", one(ha, line=ln, group="silent")["sense_weaker_pct"], pct, 0)
for ln, rho, lo, hi in (("K562", 0.00, -0.04, 0.04), ("HCT116", -0.02, -0.05, 0.02)):
    r = one(tr, line=ln)
    check(f"S1C rho {ln}", r["spearman_rho_hostsignal_vs_withinpair_diff"], rho, 2)
    check(f"S1C rho lo {ln}", r["rho_lo95"], lo, 2)
    check(f"S1C rho hi {ln}", r["rho_hi95"], hi, 2)
for ln in ("K562", "HCT116"):
    r = one(k, line=ln, set="k36 minus no_k36")
    if not (r["diff_lo95"] < 0 < r["diff_hi95"]):
        bad.append(f"S1C H3K36me3 {ln}: interval excludes zero")
    else:
        ok += 1
# S1D: composition
sh, cm = B["S1D Share from sequence"]
for ln, share, r2 in (("K562", -0.24, 0.002), ("HCT116", -0.11, 0.006)):
    r = one(cm, line=ln, group="intron_active", predictor="composition: T/A and G/C skews (+-300 nt)")
    check(f"S1D composition share {ln}", r["frac_linear"], share, 2)
    check(f"S1D cv R2 {ln}", r["cv_r2_ctrl"], r2, 3)
check("S1D Puffin-D alone K562 (35%)", one(cm, line="K562", group="intron_active", predictor="Puffin-D alone")["frac_linear"], 0.35, 2)
check("S1D Puffin-D + composition K562 (29%)",
      one(cm, line="K562", group="intron_active", predictor="Puffin-D + T/A and G/C skews (+-300 nt)")["frac_linear"], 0.29, 2)
# S1E: transposable elements
te = B["S1E Transposable elements"]
for ln, pct, lo, hi in (("K562", 31, 11, 47), ("HCT116", 25, 14, 35)):
    r = one(te[4], line=ln)
    check(f"S1E noTE {ln}", r["sense_weaker_pct"], pct, 0)
    check(f"S1E noTE lo {ln}", r["pct_lo95"], lo, 0)
    check(f"S1E noTE hi {ln}", r["pct_hi95"], hi, 0)
for ln, pct in (("K562", 25), ("HCT116", 22)):
    check(f"S1E noRepeat {ln}", one(te[5], line=ln)["sense_weaker_pct"], pct, 0)
check("S1E noRepeat K562 lo (-1)", one(te[5], line="K562")["pct_lo95"], -1, 0)
check("S1E noRepeat K562 hi (41)", one(te[5], line="K562")["pct_hi95"], 41, 0)
for ln, pct in (("K562", 63), ("HCT116", 60)):
    check(f"S1E TE against gene {ln}", one(te[1], line=ln, te_class="all")["pct_opposite_pairs"], pct, 0)
check("S1E background (57%)", one(te[1], line="K562", te_class="all")["pct_opposite_background"], 57, 0)
# S1F: census shuffle, ProCapNet, all positions: 0.8 and 1.5 percentage points
cs = B["S1F Genome-wide census"][0]
for met, v in (("top0.001", 0.8), ("top0.0003", 1.5)):
    r = one(cs, model="ProCapNet", group="genome", set="intron_protein_coding", seq="all", metric=met)
    check(f"S1F ProCapNet drop {met}", r["drop_pct_points"], v, 1)
for met in ("top0.001", "top0.0003", "top0.0001"):
    r = one(cs, model="ProCapNet", group="genome", set="intron_protein_coding", seq="unique", metric=met)
    if not (r["drop_lo95"] < 0 < r["drop_hi95"]):
        bad.append(f"S1F ProCapNet unique {met}: interval excludes zero")
    else:
        ok += 1
for met, v in (("top0.0003", 18), ("top0.0001", 23)):
    r = one(cs, model="Puffin-D", group="genome", set="intron_protein_coding", seq="all", metric=met)
    check(f"S1F Puffin-D deficit {met}", r["real_deficit_pct"], v, 0)
# S1G: silent genes, three rules present; 28% deficit of weak sites is in the main step-17 tables, not here
nb = B["S1G Neighbours of silent genes"][0]
for an in ("rt20k_pr5k", "rt50k_pr10k", "strict_rt50k_pr10k"):
    if (nb["analysis"] == an).sum() == 0:
        bad.append(f"S1G {an} missing")
    else:
        ok += 1
# S1H: methylation
st, pc, cl, rh = B["S1H DNA methylation"]
pri = "ENCODE default processing (ENCFF660IHA, ENCFF328NMN)"
alt = "earlier ENCODE processing of the same reads (ENCFF867JRG, ENCFF721JMB)"
r = one(st, processing=pri, analysis="meth_within_activity", subset="pc_active", bin="M3-M1", metric="obs r10")
check("S1H planned 4.5", r["deficit"], 4.5, 1)
check("S1H planned lo -12.9", r["lo"], -12.9, 1)
check("S1H planned hi 20.7", r["hi"], 20.7, 1)
check("S1H planned P 0.61", r["p_boot"], 0.61, 2)
r = one(st, processing=alt, analysis="meth_within_activity", subset="pc_active", bin="M3-M1", metric="obs r10")
check("S1H alt planned 1.7", r["deficit"], 1.7, 1)
check("S1H alt planned lo -15.5", r["lo"], -15.5, 1)
check("S1H alt planned hi 17.6", r["hi"], 17.6, 1)
for b_, d, lo, hi in (("M1", 16, 7, 24), ("M2", 19, 10, 28), ("M3", 20, 5, 34)):
    r = one(st, processing=pri, analysis="meth_within_activity", subset="pc_active", bin=b_, metric="obs r10")
    check(f"S1H pooled {b_}", r["deficit"], d, 0)
    check(f"S1H pooled {b_} lo", r["lo"], lo, 0)
    check(f"S1H pooled {b_} hi", r["hi"], hi, 0)
    if int(r["n_units"]) not in (3058, 3059):
        bad.append(f"S1H pooled {b_} genes {r['n_units']}")
    else:
        ok += 1
nb_ = st[(st["processing"] == pri) & (st["analysis"] == "meth_blocks") & (st["bin"].isin(["M1", "M2", "M3"]))]["n_units"]
if set(nb_.astype(int)) <= {848, 849}:
    ok += 1
else:
    bad.append(f"S1H blocks per third {sorted(set(nb_))}")
r = one(pc, set="primary", metric="obs r10", contrast="M3-M1")
check("S1H equal 24.5", r["equal_weight"], 24.5, 1)
check("S1H equal lo 1.4", r["equal_lo95"], 1.4, 1)
check("S1H equal hi 48.0", r["equal_hi95"], 48.0, 1)
check("S1H equal P 0.04", r["p_equal"], 0.04, 2)
r = one(pc, set="alternative", metric="obs r10", contrast="M3-M1")
check("S1H alt equal 19.6", r["equal_weight"], 19.6, 1)
check("S1H alt equal lo -4.6", r["equal_lo95"], -4.6, 1)
check("S1H alt equal hi 42.7", r["equal_hi95"], 42.7, 1)
check("S1H weight A1-M3 (63%)", one(cl, set="primary", cell="A1-M3")["weight_in_pooled"], 0.63, 2)
check("S1H weight A1-M1 (8%)", one(cl, set="primary", cell="A1-M1")["weight_in_pooled"], 0.08, 2)
for c_, m_ in (("A1-M1", 0.07), ("A1-M2", 0.41), ("A1-M3", 0.75)):
    check(f"S1H median meth {c_}", one(cl, set="primary", cell=c_)["meth_median"], m_, 2)
mm = cl[(cl["set"] == "primary") & (~cl["activity_third"].eq("A1"))]["meth_median"]
check("S1H median meth other thirds min (0.58)", mm.min(), 0.58, 2)
check("S1H median meth other thirds max (0.89)", mm.max(), 0.89, 2)
g = cl[cl["set"] == "primary"]["genes"].astype(int)
if g.min() == 1007 and g.max() == 1028 and g.sum() == 9176:
    ok += 1
else:
    bad.append(f"S1H genes per cell {g.min()}-{g.max()}, total {g.sum()}")
check("S1H rho A1 (0.54)", one(rh, set="primary", activity_third="A1")["rho_meth_host_activity"], 0.54, 2)
# S1I: population test
dif, cls = B["S1I Population test"]
r = one(dif, model="ProCapNet", expectation="base", set="intron_protein_coding", sense_class="creates_sense_1.131")
check("S1I ProCapNet create weak", r["excess_diff"], -0.0014, 4)
check("S1I ProCapNet create weak lo", r["diff_lo95"], -0.0025, 4)
check("S1I ProCapNet create weak hi", r["diff_hi95"], -0.0002, 4)
r = one(dif, model="ProCapNet", expectation="base", set="intergenic", sense_class="destroys_sense_2.658")
check("S1I intergenic remove medium", r["excess_diff"], -0.0099, 4)
check("S1I intergenic remove medium lo", r["diff_lo95"], -0.0173, 4)
check("S1I intergenic remove medium hi", r["diff_hi95"], -0.0022, 4)
r = one(dif, model="Puffin-D", expectation="base", set="intron_protein_coding", sense_class="gain_sense_mid")
check("S1I Puffin-D moderate raise", r["excess_diff"], 0.0034, 4)
check("S1I Puffin-D moderate raise lo", r["diff_lo95"], 0.0020, 4)
check("S1I Puffin-D moderate raise hi", r["diff_hi95"], 0.0049, 4)
check("S1I Puffin-D moderate raise n sense (1.07 M)", r["n_sense"] / 1e6, 1.07, 2)
check("S1I Puffin-D moderate raise n anti (0.73 M)", r["n_anti"] / 1e6, 0.73, 2)
r = one(dif, model="ProCapNet", expectation="base", set="intron_protein_coding", sense_class="gain_sense_mid")
check("S1I ProCapNet moderate raise", r["excess_diff"], -0.0023, 4)
check("S1I ProCapNet moderate raise lo", r["diff_lo95"], -0.0036, 4)
check("S1I ProCapNet moderate raise hi", r["diff_hi95"], -0.0010, 4)
for model, cl_, n in (("ProCapNet", "creates", 2610417), ("ProCapNet", "gain_strong", 911755), ("ProCapNet", "gain_mid", 2249974),
                      ("Puffin-D", "creates", 275674), ("Puffin-D", "gain_strong", 642746), ("Puffin-D", "gain_mid", 1800949)):
    sc = {"creates": "creates_sense_1.131" if model == "ProCapNet" else "creates_sense_0.0439",
          "gain_strong": "gain_sense_strong", "gain_mid": "gain_sense_mid"}[cl_]
    r = one(dif, model=model, expectation="base", set="intron_protein_coding", sense_class=sc)
    check(f"S1I {model} {cl_} n (Fig. 4 legend)", r["n_sense"] + r["n_anti"], n, 0)
if (dif["set"] == "intron_lncRNA").sum() > 0:
    ok += 1
else:
    bad.append("S1I lncRNA rows missing")
print(f"checked {ok + len(bad)} values: {ok} agree, {len(bad)} differ")
for b_ in bad:
    print("DIFFERS:", b_)
sys.exit(1 if bad else 0)
