"""
02_summary.py — render the T4 cumulative-horizon fusion table
=============================================================
Reads integration/T4/fusion_metrics.csv (per-seed-mean metrics) and
fusion_predictions_all_T4.csv (per-patient preds, for POOLED bACC), and
renders a PNG leaderboard. Run AFTER 01_fuse_T4_cumulative.py.

Shows BOTH:
  - TEST bACC (seed-mean ± std)  — mean of per-seed balanced accuracy. On n≈10/
    seed × 3 classes this is HIGH-VARIANCE and can rank methods misleadingly.
  - pooled bACC                  — balanced accuracy over ALL ~31 test patients
    pooled across seeds. More trustworthy; agrees with confusion_T4.png.
These two can DISAGREE on this tiny cohort (e.g. seed-mean says LR>clinical,
pooled says clinical>LR) — that is the small-n averaging artifact, not a real gain.

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
# metrics after the two bACC columns
REST_METRICS = [("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1")] + \
               [(f"f1_{c}", f"F1 {c}") for c in CLASSES] + [("val_bacc", "val bACC")]
VAR_LABEL = {
    "clinical_only": "clinical-only [BioClin-L(<3y)+ModernBERT-L(<7y) @bl]",
    "mri_only":      "MRI-only [ViT-MAE frozen/none @m12]",
    "equal_0.5":     "equal 0.5  (CL+MRI buckets)",
    "weighted_avg":  "weighted-avg bACC  (CL+MRI)",
    "weighted_avg_J": "weighted-avg Youden  (CL+MRI)",
    "elastic_net":   "elastic-net  (CL+MRI)",
    "lr":            "logistic-reg  (CL+MRI)",
}
ORDER = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg",
         "weighted_avg_J", "elastic_net", "lr"]


def _cell(m, s):
    if pd.isna(m):
        return "—"
    return f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}"


def main():
    f = HERE / "fusion_metrics.csv"
    if not f.exists():
        raise SystemExit("Run 01_fuse_T4_cumulative.py first.")
    d = pd.read_csv(f)
    src = HERE / "mri_source.txt"
    mri_label = src.read_text().strip() if src.exists() else "ViT-MAE frozen/none"
    VAR_LABEL["mri_only"] = f"MRI-only [{mri_label} @m12]"

    # pooled bACC per variant from the per-patient predictions
    pooled = {}
    ap_f = HERE / "fusion_predictions_all_T4.csv"
    if ap_f.exists():
        ap = pd.read_csv(ap_f)
        for v, g in ap.groupby("variant"):
            gg = g[g["pred"] >= 0]
            pooled[v] = balanced_accuracy_score(gg["y_true"], gg["pred"]) if len(gg) else np.nan

    body = []
    for v in ORDER:
        r = d[d["variant"] == v]
        if not len(r):
            continue
        r = r.iloc[0]
        row = [VAR_LABEL[v],
               _cell(r.get("bacc_mean"), r.get("bacc_std")),
               f"{pooled.get(v, np.nan):.3f}" if v in pooled and not np.isnan(pooled[v]) else "—"]
        row += [_cell(r.get(f"{mk}_mean"), r.get(f"{mk}_std")) for mk, _ in REST_METRICS]
        row += [f"{r['n_mean']:.0f}"]
        body.append(row)

    col_labels = (["Method", "TEST bACC (seed-mean)", "pooled bACC"]
                  + [lbl for _, lbl in REST_METRICS] + ["n"])
    n_cols = len(col_labels)
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2 for j in range(n_cols)]
    col_chars[0] += 8
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(0.12 * total, 0.7 + 0.42 * (len(body) + 1)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.5)
    for j in range(n_cols):
        tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        tab[i, 0].set_text_props(ha="left")
    present = [v for v in ORDER if len(d[d["variant"] == v])]   # body rows are in this order
    # bold best seed-mean bACC (col 1)
    sm = {v: d[d["variant"] == v]["bacc_mean"].iloc[0] for v in present}
    if sm:
        best_sm = max(sm, key=lambda v: sm[v])
        tab[present.index(best_sm) + 1, 1].set_text_props(weight="bold")
    # bold best POOLED bACC (col 2; the trustworthy metric)
    if pooled:
        cand = [v for v in present if v in pooled and not np.isnan(pooled[v])]
        if cand:
            best_v = max(cand, key=lambda v: pooled[v])
            tab[present.index(best_v) + 1, 2].set_text_props(weight="bold")

    plt.title("T4 conversion-horizon (<3y / 3-7y / ≥7y) — CUMULATIVE fusion\n"
              "buckets from cumulative binary T3 probs: clinical@bl(≤3y BioClin-L, ≤7y ModernBERT-L) "
              f"⊕ MRI@m12({mri_label})\nmean ± std across seeds (fit VAL / report TEST)",
              pad=10, fontsize=10.5)
    foot = ["bucket = [p(<3y), p(<7y)-p(<3y), 1-p(<7y)] per modality (monotone-clipped, renormalised).",
            "TEST bACC (seed-mean) = mean of per-seed bACC — HIGH-VARIANCE on n≈10/seed × 3 classes.",
            "pooled bACC = bACC over ALL ~31 test patients pooled across seeds — trustworthy; agrees with confusion_T4.png.",
            "The two CAN disagree here (seed-mean ranks LR>clinical; pooled ranks clinical>LR) — a small-n artifact, NOT a real gain.",
            "Direct 3-class T4 models were degenerate (bACC≈chance) → this cumulative construction. Eval = T4 converters,",
            "baseline split (out-of-fold). Chance bACC = 0.333. See lr_weights_T4.png + failure_table_seed*.png + confusion_T4.png."]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    png = HERE / "SUMMARY_T4_fusion.png"
    fig.savefig(png, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"wrote {png}")
    print("pooled bACC:", {v: round(pooled[v], 3) for v in ORDER if v in pooled})


if __name__ == "__main__":
    main()
