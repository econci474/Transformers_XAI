"""
02_summary.py — render the T4 cumulative-horizon fusion table (both cohorts)
============================================================================
Reads fusion_metrics.csv (per (variant, cohort) seed-mean) +
fusion_predictions_all_T4.csv (per-patient, for POOLED bACC), renders a
leaderboard with a COHORT column:
  • both-present : complete-case (only patients with an m12 MRI).
  • present-only : full clinical cohort (single modality proceeds where the
                   other is absent) — the representative, larger test set.
For each method, both cohorts are shown. MRI-only exists only on the
MRI-present subset (both-present). Shows seed-mean AND pooled bACC (the two can
disagree on this tiny cohort; pooled is the trustworthy one).

Run:  conda run -n snp python integration/T4/02_summary.py
Outputs: integration/T4/SUMMARY_T4_fusion.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score

HERE = Path(__file__).resolve().parent
CLASSES = ["h_lt3", "h_3_7", "h_ge7"]
REST_METRICS = [("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1")] + \
               [(f"f1_{c}", f"F1 {c}") for c in CLASSES] + [("val_bacc", "val bACC")]
VAR_LABEL = {
    "clinical_only": "clinical-only [BioClin-L(<3y)+ModernBERT-L(<7y)]",
    "mri_only":      "MRI-only [BrainMVP full_ft plus_original]",
    "equal_0.5":     "equal 0.5  (CL+MRI)",
    "weighted_avg":  "weighted-avg bACC  (CL+MRI)",
    "weighted_avg_J": "weighted-avg Youden  (CL+MRI)",
    "elastic_net":   "elastic-net  (CL+MRI)",
    "lr":            "logistic-reg  (CL+MRI)",
}
ORDER = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg",
         "weighted_avg_J", "elastic_net", "lr"]
COHORTS = [("both_present", "both-present"), ("present_only", "present-only")]


def _cell(m, s):
    if pd.isna(m):
        return "—"
    return f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}"


def main():
    f = HERE / "fusion_metrics.csv"
    if not f.exists():
        raise SystemExit("Run 01_fuse_T4_cumulative.py first.")
    d = pd.read_csv(f)

    pooled = {}
    ap_f = HERE / "fusion_predictions_all_T4.csv"
    if ap_f.exists():
        ap = pd.read_csv(ap_f)
        for (v, c), g in ap.groupby(["variant", "cohort"]):
            gg = g[g["pred"] >= 0]
            pooled[(v, c)] = balanced_accuracy_score(gg["y_true"], gg["pred"]) if len(gg) else np.nan

    body, kinds, bold = [], [], []          # bold = (row_idx, col_idx) to embolden
    for ckey, clabel in COHORTS:
        cd = d[d["cohort"] == ckey]
        if cd.empty:
            continue
        rows_here = [v for v in ORDER if len(cd[cd["variant"] == v])]
        # best (full-cohort variants only) per cohort, for bolding
        full = [v for v in rows_here if v != "mri_only"]
        sm = {v: cd[cd["variant"] == v]["bacc_mean"].iloc[0] for v in full}
        po = {v: pooled.get((v, ckey), np.nan) for v in full}
        best_sm = max(sm, key=lambda v: sm[v]) if sm else None
        best_po = max((v for v in full if not np.isnan(po[v])), key=lambda v: po[v], default=None)
        for v in rows_here:
            r = cd[cd["variant"] == v].iloc[0]
            ri = len(body)
            row = [clabel, VAR_LABEL[v],
                   _cell(r.get("bacc_mean"), r.get("bacc_std")),
                   f"{pooled.get((v, ckey), np.nan):.3f}" if not np.isnan(pooled.get((v, ckey), np.nan)) else "—"]
            row += [_cell(r.get(f"{mk}_mean"), r.get(f"{mk}_std")) for mk, _ in REST_METRICS]
            row += [f"{r['n_mean']:.0f}"]
            body.append(row)
            if v == best_sm:
                bold.append((ri + 1, 2))
            if v == best_po:
                bold.append((ri + 1, 3))

    col_labels = (["Cohort", "Method", "TEST bACC (seed-mean)", "pooled bACC"]
                  + [lbl for _, lbl in REST_METRICS] + ["n"])
    n_cols = len(col_labels)
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2 for j in range(n_cols)]
    col_chars[1] += 8
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(0.115 * total, 0.7 + 0.42 * (len(body) + 1)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(8.5); tab.scale(1.0, 1.5)
    for j in range(n_cols):
        tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        for j in (0, 1):
            tab[i, j].set_text_props(ha="left")
    for (i, j) in bold:
        tab[i, j].set_text_props(weight="bold")

    plt.title("T4 conversion-horizon (<3y / 3-7y / ≥7y) — CUMULATIVE fusion (both cohorts)\n"
              "buckets from cumulative binary T3 probs: clinical@bl(≤3y BioClin-L, ≤7y ModernBERT-L) "
              "⊕ MRI@m12(BrainMVP full_ft plus_original)\nmean ± std across seeds (fit VAL / report TEST)",
              pad=10, fontsize=10.5)
    foot = ["COHORTS:  both-present = complete-case (only patients with an m12 MRI; smaller).",
            "          present-only = FULL clinical cohort; where m12 MRI absent the single present modality",
            "          (clinical) proceeds (weighted-avg falls back to clinical; EN/LR uniform-impute MRI + flag).",
            "MRI-only exists only on the MRI-present subset (both-present). bold = best per cohort among full-cohort",
            "variants (excludes MRI-only — different patients). TEST bACC (seed-mean) is HIGH-VARIANCE (n≈10-16/seed);",
            "pooled bACC (over all test patients across seeds) is the trustworthy metric; the two can disagree.",
            "bucket=[p(<3y), p(<7y)-p(<3y), 1-p(<7y)] per modality (monotone-clipped). Direct 3-class T4 was",
            "degenerate -> this cumulative construction. Chance bACC=0.333. See confusion_T4.png / failure_table_seed*.png."]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    png = HERE / "SUMMARY_T4_fusion.png"
    fig.savefig(png, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"wrote {png}")
    print("pooled bACC:", {f"{v}/{c}": round(pooled[(v, c)], 3)
                           for (v, c) in pooled if not np.isnan(pooled[(v, c)])})


if __name__ == "__main__":
    main()
