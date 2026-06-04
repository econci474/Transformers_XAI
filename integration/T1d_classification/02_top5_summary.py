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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=5)
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
