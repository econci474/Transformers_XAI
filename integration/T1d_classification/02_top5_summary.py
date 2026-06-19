"""
02_top5_summary.py
==================
Summarise ALL T1d fusion strategies (cross-sectional + temporal, both MRI models) and render the
TOP-N by TEST balanced accuracy, with macro-AUC / macro-F1 / per-class F1 and explicit Clinical +
MRI model columns, plus clinical-only / MRI-only reference rows.

Clinical encoder for all T1d fusions = ModernBERT-large (full_ft); MRI model named per row.
Meta-learners fit per-seed on VAL, reported on TEST; mean ± std across seeds.

Reads every `fusion_metrics_per_seed.csv` under `integration/T1d_classification/outputs/**`.

Run:
  python integration/T1d_classification/02_top5_summary.py [--top 5]
Outputs:
  integration/T1d_classification/outputs/SUMMARY_top_fusion_by_test_bacc.{png,csv}
  integration/T1d_classification/outputs/SUMMARY_all_fusion_by_test_bacc.csv
"""
from __future__ import annotations

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
METRIC_COLS = ["bacc", "macro_auc", "macro_f1", "f1_sMCI", "f1_pMCI"]
REF_VARIANTS = {"clinical_only", "mri_only", "clinical_temporal_EN"}

METHOD = {
    ("convex_present", "equal"): "class-wise present-only (equal weights)",
    ("convex_present", "tuned_bacc"): "class-wise present-only (val-tuned w, bACC)",
    ("convex_present", "tuned_youden"): "class-wise present-only (val-tuned w, Youden)",
    ("convex_complete", "equal"): "class-wise complete-case (equal weights)",
    ("convex_complete", "tuned_bacc"): "class-wise complete-case (val-tuned w, bACC)",
    ("weighted_avg", "complete_case"): "weighted-avg 2-source (val-tuned w, bACC)",
    ("weighted_avg_J", "complete_case"): "weighted-avg 2-source (val-tuned w, Youden)",
    ("elastic_net", "complete_case"): "elastic-net stacker (complete-case)",
    ("elastic_net", "impute_indicator"): "elastic-net stacker (impute+indicator)",
    ("elastic_net", "present_only"): "elastic-net stacker (present-only)",
    ("elastic_net", "zero_fallback"): "elastic-net stacker (zero-fallback)",
    ("lr", "complete_case"): "logistic-regression stacker (complete-case)",
    ("avg", "equal"): "hierarchical avg (equal 0.5/0.5)",
    ("avg", "tuned_bacc"): "hierarchical avg (val-tuned w, bACC)",
    ("avg", "tuned_youden"): "hierarchical avg (val-tuned w, Youden)",
    ("avg", "present_only_equal"): "hierarchical avg present-only (0.5/0.5)",
    ("clinical_only", "-"): "clinical-only (reference)",
    ("mri_only", "-"): "MRI-only (reference)",
    ("clinical_temporal_EN", "impute_indicator"): "clinical-temporal EN (reference)",
}


def mri_model(group):
    g = group.lower()
    if "braindino" in g:
        return "BrainDINO frozen-none"
    if "brainmvp" in g:
        return "BrainMVP plus_original"
    return "?"


def clin_visits(group, variant):
    if variant == "clinical_only":
        return "bl" if group.startswith("CLbl") else "m12"
    if group.startswith("CLbl"):
        return "bl"
    if group.startswith("CLm12"):
        return "m12"
    if group.split("/")[0] in ("flat", "hierarchical", "class_wise_fusion"):
        return "bl+m06+m12"
    return "?"


def collect():
    rows = []
    for f in glob.glob(os.path.join(OUT_ROOT, "**", "fusion_metrics_per_seed.csv"), recursive=True):
        rel = os.path.relpath(os.path.dirname(f), OUT_ROOT).replace("\\", "/")
        d = pd.read_csv(f)
        for c in METRIC_COLS + ["n"]:
            if c not in d.columns:
                d[c] = np.nan
        d["group"] = rel
        rows.append(d[["group", "variant", "missing", "seed", "n"] + METRIC_COLS])
    if not rows:
        raise SystemExit(f"No fusion_metrics_per_seed.csv under {OUT_ROOT}")
    return pd.concat(rows, ignore_index=True)


# --------------------------------------------------------------------------- #
# Cleaned house-style table (Method / MRI / Aggregation), mirrors the T2 cleaned table.
# --------------------------------------------------------------------------- #
CLEAN_METRICS = [("bacc", "Test bACC"), ("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1"),
                 ("f1_sMCI", "F1 sMCI"), ("f1_pMCI", "F1 pMCI")]
# Per-class TEST composition (mean per seed), verified from the canonical splits.
N_TEST_SUBTITLE = (
    "N(test), mean per seed.   CL_baseline: n=33 (sMCI 21 / pMCI 12)   ·   "
    "CL_12 months: n=30 (sMCI 19 / pMCI 11)   ·   MRI_12 months: n=24 (sMCI 19 / pMCI 5)   ·   "
    "all-available: n=32   ·   complete-case: n=20")


def cleaned_name(variant, clinv, mri_used):
    """Compact thesis Method name; clinical = ModernBERT-large/full_ft -> 'MBERT-L-ft'."""
    if variant == "mri_only":
        return "MRI-only [T1d @m12]"
    if variant == "clinical_only":
        return f"MBERT-L-ft T1d @{clinv}"
    if variant == "clinical_temporal_EN":
        return "MBERT-L-ft T1d @bl+m06+m12"
    base = f"MBERT-L-ft T1d @{clinv}"
    return base + (" + MRI @m12" if mri_used else "")


def agg_method(variant, missing):
    """The fusion combiner / aggregation architecture (its own column)."""
    if variant in ("clinical_only", "mri_only"):
        return "—"
    if variant == "clinical_temporal_EN":
        return "elastic-net (temporal)"
    if variant == "convex_present":
        return {"equal": "class-wise, mean",
                "tuned_bacc": "class-wise, all-available (val bACC)",
                "tuned_youden": "class-wise, all-available (Youden J)"}.get(missing, "class-wise, all-available")
    if variant == "convex_complete":
        return {"equal": "class-wise, complete-case (equal)",
                "tuned_bacc": "class-wise, complete-case (val bACC)"}.get(missing, "class-wise, complete-case")
    if variant == "weighted_avg":
        return "weighted-avg (val bACC)"
    if variant == "weighted_avg_J":
        return "weighted-avg (Youden J)"
    if variant == "avg":
        return {"equal": "hierarchical, mean (complete-case)", "tuned_bacc": "hierarchical (val bACC)",
                "tuned_youden": "hierarchical (Youden J)",
                "present_only_equal": "hierarchical, mean"}.get(missing, "hierarchical")
    if variant == "elastic_net":
        return f"flat stack, elastic-net ({missing})"
    if variant == "lr":
        return f"flat stack, logreg ({missing})"
    return f"{variant}/{missing}"


def _full_cohort_clinical_only():
    """Recompute clinical-only on the FULL per-visit clinical test cohort (no MRI filter), reusing the
    fusion harness loaders. Returns agg-style rows (one per visit)."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("t1d_fuse", os.path.join(here, "01_fuse_temporal_T1d.py"))
    fz = importlib.util.module_from_spec(spec); spec.loader.exec_module(fz)
    h = fz.h
    rows = []
    for visit, root in [("bl", fz.CL_BL_ROOT), ("m12", fz.CL_M12_ROOT)]:
        per = {c: [] for c in METRIC_COLS}; ns = []
        for seed in (0, 1, 2):
            d = h.load_clinical(seed, root, "ModernBERT-large", h.CLIN_STRAT, task="T1d_pmci_smci")
            t = d[d["split"] == "test"]
            y = t["y_clin"].astype(int).to_numpy(); P = t[["cp0", "cp1"]].to_numpy(float)
            r = h.metric_row(y, P.argmax(1), P)
            for c in METRIC_COLS:
                per[c].append(r[c])
            ns.append(len(t))
        row = {"variant": "clinical_only", "missing": "-", "clinv": visit, "mri": "—",
               "n_test": float(np.mean(ns))}
        for c in METRIC_COLS:
            row[f"{c}_m"] = float(np.mean(per[c])); row[f"{c}_s"] = float(np.std(per[c], ddof=1))
        rows.append(row)
    return pd.DataFrame(rows)


def render_cleaned(agg, top_n=8):
    # full-cohort clinical-only REPLACES the folder-derived (both-present) clinical-only rows
    agg = agg[agg.variant != "clinical_only"].copy()
    clin_full = _full_cohort_clinical_only()

    # fusion rows ranked by bACC on top; references grouped at the bottom (clinical-only @m12, @bl,
    # then the two MRI-only together with BrainMVP before BrainDINO).
    fusion_sel = agg[~agg.variant.isin(REF_VARIANTS)].sort_values("bacc_m", ascending=False).head(top_n)
    mref = agg[agg.variant == "mri_only"].drop_duplicates("mri").copy()
    mref["_ord"] = mref.mri.map(lambda m: 0 if "BrainMVP" in str(m) else 1)
    mref = mref.sort_values("_ord")                                         # BrainMVP before BrainDINO
    clin_full = clin_full.sort_values("bacc_m", ascending=False)           # @m12 (0.780) then @bl (0.727)
    # fusion ranked on top; then MRI-only block; then clinical-only block (dotted rule before each block)
    sel = pd.concat([fusion_sel, mref, clin_full], ignore_index=True)

    sel["display"] = [cleaned_name(v, cv, (m in ("BrainMVP plus_original", "BrainDINO frozen-none")))
                      for v, cv, m in zip(sel.variant, sel.clinv, sel.mri)]
    sel["aggregation"] = [agg_method(v, mi) for v, mi in zip(sel.variant, sel.missing)]
    sel["mri_disp"] = ["BrainDINO" if "BrainDINO" in str(m) else ("BrainMVP" if "BrainMVP" in str(m) else "—")
                       for m in sel.mri]
    sel.loc[sel.variant.isin(["clinical_only", "clinical_temporal_EN"]), "mri_disp"] = "—"
    sel = sel.reset_index(drop=True)
    # a dotted rule above the FIRST MRI-only row and above the FIRST clinical-only row
    sel["rule_above"] = False
    nf = len(fusion_sel)
    if len(mref):
        sel.loc[nf, "rule_above"] = True
    if len(clin_full):
        sel.loc[nf + len(mref), "rule_above"] = True

    out_cols = (["display", "mri_disp", "aggregation"]
                + [f"{c}_m" for c in METRIC_COLS] + [f"{c}_s" for c in METRIC_COLS] + ["n_test"])
    sel[out_cols].rename(columns={"display": "method", "mri_disp": "MRI", "aggregation": "Aggregation"}) \
        .to_csv(os.path.join(OUT_ROOT, "SUMMARY_top_fusion_by_test_bacc_cleaned.csv"), index=False)

    render_house(sel, os.path.join(OUT_ROOT, "SUMMARY_top_fusion_by_test_bacc_cleaned.png"),
                 N_TEST_SUBTITLE, show_n=True)
    render_house(sel, os.path.join(OUT_ROOT, "SUMMARY_top_fusion_by_test_bacc_cleaned_no_n.png"),
                 N_TEST_SUBTITLE, show_n=False)

    # ---- int_t1d_simple: csv + png + latex (full) and a filtered headline-methods version ----
    def _block(r):
        if str(r["display"]).startswith("MRI-only"):
            return "mri"
        return "clin" if str(r["mri_disp"]) == "—" else "fusion"

    def _with_rules(d):
        d = d.reset_index(drop=True).copy()
        d["block"] = [_block(r) for _, r in d.iterrows()]
        dec = [_t1d_decompose(r["display"], r["mri_disp"]) for _, r in d.iterrows()]   # int_t2_simple cols
        d["Method"] = [x[0] for x in dec]
        d["Clinical (Task)"] = [x[1] for x in dec]
        d["MRI (Task)"] = [x[2] for x in dec]
        d["Aggregation"] = d["aggregation"]
        ra, prev = [], None
        for b in d["block"]:
            ra.append(prev is not None and b != prev); prev = b
        d["rule_above"] = ra
        return d

    KEEP_AGG = {"class-wise, mean", "weighted-avg (val bACC)", "hierarchical, mean"}
    full_sel = _with_rules(sel)
    filt_sel = _with_rules(sel[(sel.apply(_block, axis=1) != "fusion")
                               | sel["aggregation"].isin(KEEP_AGG)])

    def _csv(d, name):
        cols = _T1D_LEAD + [f"{c}_m" for c in METRIC_COLS] + [f"{c}_s" for c in METRIC_COLS] + ["n_test"]
        d[cols].rename(columns={"n_test": "n"}).to_csv(os.path.join(OUT_ROOT, name), index=False)
    _csv(full_sel, "int_t1d_simple.csv")
    _csv(filt_sel, "int_t1d_simple_filtered.csv")
    render_int_t1d(full_sel, os.path.join(OUT_ROOT, "int_t1d_simple.png"), N_TEST_SUBTITLE, show_n=True)
    render_int_t1d(full_sel, os.path.join(OUT_ROOT, "int_t1d_simple_no_n.png"), N_TEST_SUBTITLE, show_n=False)
    render_int_t1d(filt_sel, os.path.join(OUT_ROOT, "int_t1d_simple_filtered.png"), N_TEST_SUBTITLE, show_n=True)
    render_int_t1d(filt_sel, os.path.join(OUT_ROOT, "int_t1d_simple_filtered_no_n.png"), N_TEST_SUBTITLE,
                   show_n=False)
    write_t1d_latex(full_sel, os.path.join(OUT_ROOT, "int_t1d_simple.tex"), filtered=False)
    write_t1d_latex(filt_sel, os.path.join(OUT_ROOT, "int_t1d_simple_filtered.tex"), filtered=True)

    # ---- int_t1d_filtered: the single TOP BrainMVP fusion + its corresponding unimodal references
    #      (clinical-only @m12/@bl + MRI-only BrainMVP); drops the BrainDINO rows. ----
    top_bmvp = (full_sel[(full_sel.block == "fusion") & (full_sel.mri_disp == "BrainMVP")]
                .sort_values("bacc_m", ascending=False).head(1))
    mri_bmvp = full_sel[(full_sel.block == "mri") & (full_sel.mri_disp == "BrainMVP")]
    clin_ref = full_sel[full_sel.block == "clin"]
    min_sel = _with_rules(pd.concat([top_bmvp, mri_bmvp, clin_ref], ignore_index=True))
    _csv(min_sel, "int_t1d_filtered.csv")
    render_int_t1d(min_sel, os.path.join(OUT_ROOT, "int_t1d_filtered.png"), N_TEST_SUBTITLE, show_n=True)
    render_int_t1d(min_sel, os.path.join(OUT_ROOT, "int_t1d_filtered_no_n.png"), N_TEST_SUBTITLE, show_n=False)
    write_t1d_latex(min_sel, os.path.join(OUT_ROOT, "int_t1d_filtered.tex"), filtered=True,
                    label="tab:t1d_late_fusion_min",
                    caption_note=r" Top BrainMVP fusion vs its unimodal (clinical-only and MRI-only) "
                                 r"references; full results in Table~\ref{tab:t1d_late_fusion}.")

    print("=" * 90)
    print(f"  T1d cleaned summary ({len(sel)} rows) — full-cohort clinical-only")
    print("=" * 90)
    print(sel[["mri_disp", "aggregation", "display", "bacc_m", "n_test"]].to_string(index=False))
    print(f"\n  wrote -> SUMMARY_top_fusion_by_test_bacc_cleaned{{,_no_n}}.png (+ .csv)")
    return 0


_T1D_VIS = {"bl+m06+m12": "longitudinal", "bl": "at baseline", "m12": "at 12 months"}
_T1D_STRAT = {"BrainMVP": "ft", "BrainDINO": "frozen"}
_T1D_LEAD = ["Method", "Clinical (Task)", "MRI (Task)", "Aggregation"]


def _t1d_decompose(display, mri_disp):
    """Split a T1d row into (Method plain-English, Clinical (Task), MRI (Task)) — int_t2_simple layout."""
    disp = str(display); model = str(mri_disp)
    strat = _T1D_STRAT.get(model, "")
    mri_tag = (f"{model}-{strat} (T1d)" if strat else f"{model} (T1d)") if model not in ("—", "nan", "") else "—"
    if disp.startswith("MRI-only"):
        return "MRI at 12 months", "—", mri_tag
    m = re.search(r"T1d @(\S+?)(?: \+ MRI @m12)?$", disp)
    phrase = _T1D_VIS.get(m.group(1), m.group(1)) if m else ""
    if "+ MRI @m12" in disp:
        return f"CL {phrase} + MRI at 12 months", "ModernBERT-L-ft (T1d)", mri_tag
    return f"CL {phrase}", "ModernBERT-L-ft (T1d)", "—"


def render_int_t1d(df, out_path, subtitle="", show_n=True):
    """house-style PNG with the int_t2_simple columns: Method | Clinical (Task) | MRI (Task) |
    Aggregation + 5 metrics. Bold-best per metric col; dashed group rules from `block` changes."""
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    headers = _T1D_LEAD + [lbl for _, lbl in CLEAN_METRICS] + (["n"] if show_n else [])
    N_LEAD = len(_T1D_LEAD)
    body, numeric, rules = [], [], []
    prev_b = None
    for _, r in df.iterrows():
        cells = [r["Method"], r["Clinical (Task)"], r["MRI (Task)"], r["Aggregation"]]
        nums = []
        for key, _ in CLEAN_METRICS:
            v = r.get(f"{key}_m"); s = r.get(f"{key}_s")
            if pd.isna(v):
                cells.append("—"); nums.append(np.nan)
            else:
                cells.append(f"{v:.3f}" + (f" ± {s:.3f}" if pd.notna(s) else "")); nums.append(float(v))
        if show_n:
            cells.append(f"{r['n_test']:.0f}")
        body.append(cells); numeric.append(nums)
        rules.append(prev_b is not None and r["block"] != prev_b); prev_b = r["block"]
    numeric = np.array(numeric, float)
    best = [int(np.nanargmax(numeric[:, j])) if not np.all(np.isnan(numeric[:, j])) else -1
            for j in range(len(CLEAN_METRICS))]

    COL_W = [3.50, 3.05, 2.70, 3.45, 1.58, 1.55, 1.55, 1.48, 1.48] + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_H = 0.40 if subtitle else 0.0
    TITLE_H, HEAD_H, ROW_H = 1.25 + SUB_H, 0.40, 0.44
    TOP_PAD, BOT_PAD = 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(COL_W))]

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD
    y_title_top = y; y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) / 2,
            "T1d Late Fusion (sMCI vs pMCI) — clinical ⊕ MRI (simplified)\n"
            "mean ± std across seeds 0/1/2 · fusion ranked by Test balanced accuracy; references below\n"
            "Method = which clinical visits + MRI timepoint are fused\n"
            "Clinical (Task) / MRI (Task) = exact model and task formulation",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    if subtitle:
        ax.text(cx, y + SUB_H / 2, subtitle, ha="center", va="center", fontsize=8.5, fontstyle="italic")
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    LEFT_COLS = {0, 1, 2, 3}
    y_head_top = y; y -= HEAD_H
    ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=9.5, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y; y -= ROW_H
        ymid = (yr_top + y) / 2
        for j in range(len(headers)):
            metric_j = j - N_LEAD
            bold = (0 <= metric_j < len(CLEAN_METRICS) and best[metric_j] == i)
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ymid, cells[j], ha="left", va="center", fontsize=9.0)
            else:
                ax.text(col_cx[j], ymid, cells[j], ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bold else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def write_t1d_latex(df, out_path, filtered=False, label=None, caption_note=None):
    """LaTeX (sidewaystable, booktabs, multicolumn title + group labels, bottom caption with N(test) +
    abbreviations), mirroring the int_t2_simple style. filtered=True drops the F1 sMCI / F1 pMCI columns;
    label / caption_note override the defaults."""
    def esc(s):
        return (str(s).replace("—", "--").replace("&", r"\&").replace("%", r"\%")
                .replace("_", r"\_").replace("#", r"\#"))
    metric_cols = [("bacc", "Test bACC"), ("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1")]
    if not filtered:
        metric_cols += [("f1_sMCI", "F1 sMCI"), ("f1_pMCI", "F1 pMCI")]
    ncol = len(_T1D_LEAD) + len(metric_cols)
    mat = np.array([[r[f"{k}_m"] for k, _ in metric_cols] for _, r in df.iterrows()], float)
    best = [int(np.nanargmax(mat[:, j])) if not np.all(np.isnan(mat[:, j])) else -1
            for j in range(len(metric_cols))]
    BLOCK = {"fusion": r"Clinical $\oplus$ MRI late fusion", "mri": "MRI-only references",
             "clin": "Clinical-only references"}

    def cell(m, s, bold):
        if pd.isna(m):
            return "--"
        txt = f"{m:.3f} \\pm {s:.3f}" if pd.notna(s) else f"{m:.3f}"
        return f"$\\mathbf{{{txt}}}$" if bold else f"${txt}$"

    head = [r"\textbf{Method}", r"\textbf{Clinical (Task)}", r"\textbf{MRI (Task)}",
            r"\textbf{Aggregation}"] + [r"\textbf{" + lbl + "}" for _, lbl in metric_cols]
    L = [r"% T1d Late Fusion -- clinical + MRI late fusion." + (" Filtered (headline methods)." if filtered else ""),
         r"% Requires: \usepackage{rotating, booktabs, graphicx}",
         r"\begin{sidewaystable}[ht]", r"\centering",
         r"\resizebox{\textwidth}{!}{%", r"\normalsize",
         r"\begin{tabular}{llll " + "c" * len(metric_cols) + "}", r"\toprule",
         r"\multicolumn{" + str(ncol) + r"}{c}{\textbf{T1d Late Fusion (sMCI vs pMCI)}} \\", r"\midrule",
         " & ".join(head) + r" \\", r"\midrule"]
    prev = None
    for i, (_, r) in enumerate(df.iterrows()):
        b = r["block"]
        if b != prev:
            if prev is not None:
                L.append(r"\midrule")
            L.append(r"\multicolumn{" + str(ncol) + r"}{l}{\textit{" + BLOCK[b] + r"}} \\")
        cells = [esc(r["Method"]), esc(r["Clinical (Task)"]), esc(r["MRI (Task)"]), esc(r["Aggregation"])]
        cells += [cell(r[f"{k}_m"], r.get(f"{k}_s"), best[j] == i) for j, (k, _) in enumerate(metric_cols)]
        L.append(" & ".join(cells) + r" \\")
        prev = b
    cap = (r"T1d late fusion (sMCI vs pMCI), mean $\pm$ std across seeds 0/1/2, ranked by Test balanced "
           r"accuracy; fusion rows on top, single-modality references below. \textit{Method} indicates "
           r"which clinical visits and the MRI timepoint are fused. \textit{Clinical (Task)} / "
           r"\textit{MRI (Task)} give the exact model and task formulation: T1d = binary sMCI vs pMCI; "
           r"clinical = ModernBERT-large; MRI = BrainMVP or BrainDINO; -ft = full fine-tune, -frozen = "
           r"frozen encoder. N(test), mean per seed: CL baseline $n=33$ (sMCI 21 / pMCI 12); CL 12 months "
           r"$n=30$ (sMCI 19 / pMCI 11); MRI 12 months $n=24$ (sMCI 19 / pMCI 5); all-available $n=32$; "
           r"complete-case $n=20$. Bold = best per column.")
    if caption_note is None:
        caption_note = (r" Headline methods; the full ablation is in Table~\ref{tab:t1d_late_fusion}."
                        if filtered else "")
    cap += caption_note
    lab = label if label is not None else ("tab:t1d_late_fusion_filtered" if filtered
                                           else "tab:t1d_late_fusion")
    L += [r"\bottomrule", r"\end{tabular}%", r"}", r"\caption{" + cap + r"}",
          r"\label{" + lab + r"}", r"\end{sidewaystable}"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  TEX: {out_path}")


def render_house(sel, out_path, subtitle="", show_n=True):
    """Manual house-style table (DejaVu Serif, bordered, ruled, italic headers, dashed dividers,
    bold best-per-column) — mirrors integration/T2_classification/04_summary_top_methods.py."""
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    headers = ["Method", "MRI", "Aggregation"] + [lbl for _, lbl in CLEAN_METRICS] + (["n"] if show_n else [])
    N_LEAD = 3

    body, numeric, rules = [], [], []
    for _, r in sel.iterrows():
        cells = [r["display"], r["mri_disp"], r["aggregation"]]
        nums = []
        for key, _ in CLEAN_METRICS:
            v = r.get(f"{key}_m"); s = r.get(f"{key}_s")
            if pd.isna(v):
                cells.append("—"); nums.append(np.nan)
            else:
                cells.append(f"{v:.3f}" + (f" ± {s:.3f}" if pd.notna(s) else "")); nums.append(float(v))
        if show_n:
            cells.append(f"{r['n_test']:.0f}")
        body.append(cells); numeric.append(nums); rules.append(bool(r.get("rule_above", False)))
    numeric = np.array(numeric, float)
    best = [int(np.nanargmax(numeric[:, j])) if not np.all(np.isnan(numeric[:, j])) else -1
            for j in range(len(CLEAN_METRICS))]

    # every metric cell is "0.XXX ± 0.XXX" -> give each metric column room for the full value
    COL_W = [4.55, 1.05, 3.25, 1.58, 1.55, 1.55, 1.48, 1.48] + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_H = 0.40 if subtitle else 0.0
    TITLE_H, HEAD_H, ROW_H = 1.15 + SUB_H, 0.40, 0.44
    TOP_PAD, BOT_PAD = 0.12, 0.12
    n_rows = len(body)
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * n_rows + BOT_PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(COL_W))]

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD
    y_title_top = y
    y -= TITLE_H
    cx_all = (LEFT + RIGHT) / 2
    ax.text(cx_all, y + SUB_H + (TITLE_H - SUB_H) / 2,
            "T1d Late Fusion (sMCI vs pMCI)\n"
            "mean ± std across seeds 0/1/2 · ranked by Test balanced accuracy\n"
            "clinical = ModernBERT-large (full fine-tune, \"MBERT-L-ft\")\n"
            "MRI = BrainMVP (full fine-tune / plus-original) or BrainDINO (frozen / none)",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    if subtitle:
        ax.text(cx_all, y + SUB_H / 2, subtitle, ha="center", va="center",
                fontsize=8.5, fontstyle="italic")
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    LEFT_COLS = {0, 2}
    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=9.5, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:                                        # dotted rule above MRI-only / clinical-only blocks
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y
        y -= ROW_H
        ymid = (yr_top + y) / 2
        for j in range(len(headers)):
            metric_j = j - N_LEAD
            bold = (0 <= metric_j < len(CLEAN_METRICS) and best[metric_j] == i)
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ymid, cells[j], ha="left", va="center", fontsize=9.0)
            else:
                ax.text(col_cx[j], ymid, cells[j], ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bold else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--cleaned", action="store_true",
                    help="Render a cleaned house-style copy (Method/MRI/Aggregation columns, "
                         "MBERT-L-ft naming, both MRI models, full-cohort clinical-only) -> "
                         "SUMMARY_top_fusion_by_test_bacc_cleaned.{png,csv} (+ _no_n). Leaves the "
                         "existing SUMMARY_top/all outputs untouched.")
    args = ap.parse_args()

    df = collect()
    aggspec = {f"{c}_m": (c, "mean") for c in METRIC_COLS}
    aggspec.update({f"{c}_s": (c, "std") for c in METRIC_COLS})
    aggspec["n_test"] = ("n", "mean")
    agg = df.groupby(["group", "variant", "missing"]).agg(**aggspec).reset_index()
    agg = agg[agg["bacc_m"].notna()].copy()
    agg["method"] = [METHOD.get((v, m), f"{v}/{m}") for v, m in zip(agg.variant, agg.missing)]
    agg["mri"] = agg.group.map(mri_model)
    agg["clinv"] = [clin_visits(g, v) for g, v in zip(agg.group, agg.variant)]

    if args.cleaned:
        return render_cleaned(agg, top_n=max(args.top, 8))

    agg.sort_values("bacc_m", ascending=False).to_csv(
        os.path.join(OUT_ROOT, "SUMMARY_all_fusion_by_test_bacc.csv"), index=False)

    fusion = agg[~agg.variant.isin(REF_VARIANTS)].sort_values("bacc_m", ascending=False)
    top = fusion.head(args.top).copy()
    cref = (agg[agg.variant == "clinical_only"].sort_values("bacc_m", ascending=False)
            .drop_duplicates("clinv"))
    mref = (agg[agg.variant == "mri_only"].sort_values("bacc_m", ascending=False)
            .drop_duplicates("mri"))
    top.to_csv(os.path.join(OUT_ROOT, "SUMMARY_top_fusion_by_test_bacc.csv"), index=False)

    # ---- assemble rows for the figure ----
    def cell(m, s):
        if pd.isna(m):
            return "-"
        return f"{m:.3f} ± {s:.3f}" if not pd.isna(s) else f"{m:.3f}"

    def metric_cells(r):
        return [cell(r["bacc_m"], r["bacc_s"]), cell(r["macro_auc_m"], r["macro_auc_s"]),
                cell(r["macro_f1_m"], r["macro_f1_s"]), cell(r["f1_sMCI_m"], r["f1_sMCI_s"]),
                cell(r["f1_pMCI_m"], r["f1_pMCI_s"]), f"{r['n_test']:.0f}"]

    body, kinds = [], []
    for rank, (_, r) in enumerate(top.iterrows(), 1):
        body.append([str(rank), r["method"], f"ModernBERT-large @ {r['clinv']}",
                     f"{r['mri']} @ m12"] + metric_cells(r))
        kinds.append("fusion")
    for _, r in pd.concat([cref, mref]).iterrows():
        is_clin = r["variant"] == "clinical_only"
        body.append(["ref", r["method"],
                     f"ModernBERT-large @ {r['clinv']}" if is_clin else "—",
                     "—" if is_clin else f"{r['mri']} @ m12"] + metric_cells(r))
        kinds.append("ref")

    col_labels = ["#", "fusion method", "Clinical encoder", "MRI model",
                  "TEST bACC", "macro-AUC", "macro-F1", "F1 sMCI", "F1 pMCI", "n"]
    n_cols = len(col_labels)
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2
                 for j in range(n_cols)]
    # give the wide text columns (method, Clinical encoder, MRI model) extra breathing room
    for j in (1, 2, 3):
        col_chars[j] += 8
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(0.11 * total, 0.6 + 0.36 * (len(body) + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center",
                   colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.4)
    for j in range(n_cols):
        tab[0, j].set_facecolor("#d9d9d9"); tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        for j in (1, 2, 3):
            tab[i, j].set_text_props(ha="left")
        if kinds[i - 1] == "ref":
            for j in range(n_cols):
                tab[i, j].set_facecolor("#f0f0f0")
    tab[1, 4].set_text_props(weight="bold")  # best TEST bACC

    plt.title(f"T1d (pMCI vs sMCI) -- top {args.top} fusion strategies by TEST balanced accuracy\n"
              "clinical encoder = ModernBERT-large (full_ft); fit per-seed on VAL, reported on TEST "
              "(mean ± std across seeds)", pad=8, fontsize=11)
    foot = [
        "Ranked by TEST balanced accuracy (held-out). Cohorts are small (test n ~ 20-32) so std is",
        "non-trivial; equal-weight / reference strategies are the most trustworthy (no val-tuning).",
        "Fusion sources: cross-sectional weighted-avg = CL_<visit> + MRI_m12 (1 scalar w);",
        "  class-wise / temporal = CL_bl + CL_m06 + CL_m12 + MRI_m12 (convex weights, sum=1).",
        "  'equal weights' = 1/k over the present sources; 'val-tuned wts' = grid/simplex search on VAL.",
        "  present-only = drop a missing source and renormalise weights over those present.",
        "References: clinical-only = argmax of the clinical encoder at that visit; MRI-only = argmax MRI_m12.",
    ]
    ax.text(0.0, -0.05, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, family="monospace", linespacing=1.4)
    png = os.path.join(OUT_ROOT, "SUMMARY_top_fusion_by_test_bacc.png")
    fig.savefig(png, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)

    print("=" * 80)
    print(f"  TOP {args.top} T1d fusion strategies by TEST balanced accuracy "
          f"(of {len(fusion)} fusion strategies) + references")
    print("=" * 80)
    for row in body:
        print(f"  {row[0]:>3} {row[1]:<44} {row[4]:>16} AUC={row[5]:<14} n={row[9]}")
    print(f"\n  wrote: {png}")


if __name__ == "__main__":
    main()
