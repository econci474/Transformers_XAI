"""
03_lr_weights_and_failures.py — interpret the T4 cumulative LR stacker
======================================================================
Answers "how does logistic regression reach bACC≈0.617 on T4?" with:

  lr_weights_T4.png   — heatmap of the (illustrative, pooled-VAL) multinomial-LR
                        coefficients: 3 horizon classes × 6 bucket features.

  failure_table_seed{0,1,2}.png + failure_summary.png — per-patient colour-coded
                        failure TABLES (rendered by the shared T2 renderer
                        integration/T2_classification/02_render_failure_table.py),
                        columns: fused(LR) / clinical / MRI predicted class
                        (green=correct, red=wrong), + clinical-vs-MRI agreement.

  confusion_T4.png    — SEPARATE plot: true×predicted confusion matrices for
                        clinical-only and LR-fused (pooled over seeds).

Reads lr_coefficients_T4.csv + fusion_predictions_T4.csv (from 01_fuse_T4_cumulative.py).
  conda run -n snp python integration/T4/03_lr_weights_and_failures.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score, confusion_matrix

HERE = Path(__file__).resolve().parent
CLASSES = ["h_lt3", "h_3_7", "h_ge7"]
NAMES = {0: "<3y", 1: "3-7y", 2: "≥7y"}          # readable horizon labels
FEAT_ORDER = ["CL_<3y", "CL_3-7y", "CL_>=7y", "MRI_<3y", "MRI_3-7y", "MRI_>=7y"]


def plot_weights():
    f = HERE / "lr_coefficients_T4.csv"
    if not f.exists():
        print("  [skip] no lr_coefficients_T4.csv"); return
    d = pd.read_csv(f)
    M = (d.pivot_table(index="class_name", columns="feature", values="coef", aggfunc="mean")
         .reindex(index=[c for c in CLASSES if c in d["class_name"].unique()], columns=FEAT_ORDER))
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    fig.subplots_adjust(bottom=0.34, top=0.84)
    vmax = np.nanmax(np.abs(M.values)) or 1.0
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(FEAT_ORDER)))
    ax.set_xticklabels(FEAT_ORDER, rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(list(M.index))
    ax.set_ylabel("LR output class")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > 0.6 * vmax else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="coefficient")
    ax.axvline(2.5, color="k", lw=1.2)                # CL | MRI block divider
    ax.set_title("T4 cumulative LR stacker — coefficients (pooled-VAL, illustrative)\n"
                 "positive ⇒ feature pushes toward that horizon class", fontsize=11)
    fig.text(0.5, 0.04, "Vertical line splits clinical (left 3) from MRI (right 3) bucket probs.   "
             "Tiny VAL (n≈9/seed) ⇒ weights are indicative, high-variance.",
             ha="center", fontsize=8)
    out = HERE / "lr_weights_T4.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.1); plt.close(fig)
    print(f"  wrote {out}")


def write_failure_table():
    """Write a T2-renderer-compatible failure_table.csv with readable class names
    + a 'fused' (LR) modality column, then render it via the shared T2 renderer."""
    f = HERE / "fusion_predictions_T4.csv"
    if not f.exists():
        print("  [skip] no fusion_predictions_T4.csv"); return
    d = pd.read_csv(f)
    # 'fused' = weighted-avg (the best/chosen T4 fusion; LR/EN overfit the tiny val).
    fcol = "wavg_pred" if "wavg_pred" in d.columns else "lr_pred"
    fcor = "wavg_correct" if "wavg_correct" in d.columns else "lr_correct"
    mp = d["mri_present"] if "mri_present" in d.columns else pd.Series(True, index=d.index)
    out = pd.DataFrame({
        "seed": d["seed"], "Patient_ID": d["Patient_ID"],
        "y_true": d["y_true"].map(NAMES),
        "fused_pred": d[fcol].map(NAMES), "fused_correct": d[fcor],
        "clinical_pred": d["clinical_pred"].map(NAMES), "clinical_correct": d["clinical_correct"],
        # present-only: MRI absent -> NaN pred (renders 'absent', grey)
        "mri_pred": d["mri_pred"].map(NAMES), "mri_correct": d["mri_correct"].where(mp, other=np.nan),
    })
    # clinical-vs-MRI complementarity (matches the T2 failure-table 'agreement').
    def agree(r, present):
        if not present:
            return "mri_absent"
        c, m = bool(r["clinical_correct"]), bool(r["mri_correct"])
        return ("both_correct" if c and m else "clinical_only" if c else
                "mri_only" if m else "both_wrong")
    out["agreement"] = [agree(r, bool(p)) for (_, r), p in zip(out.iterrows(), mp)]
    csv = HERE / "failure_table.csv"
    out.to_csv(csv, index=False)

    # reuse the shared T2 renderer (schema-agnostic).
    ft_path = HERE.parent / "T2_classification" / "02_render_failure_table.py"
    spec = importlib.util.spec_from_file_location("_ft", ft_path)
    ft = importlib.util.module_from_spec(spec); sys.modules["_ft"] = ft
    spec.loader.exec_module(ft)
    ft.render_csv(str(csv))


def plot_confusion():
    f = HERE / "fusion_predictions_T4.csv"
    if not f.exists():
        print("  [skip] no fusion_predictions_T4.csv"); return
    d = pd.read_csv(f)
    fcol = "wavg_pred" if "wavg_pred" in d.columns else "lr_pred"
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2))
    for ax, (pred_col, name) in zip(axes, [("clinical_pred", "clinical-only"),
                                           (fcol, "weighted-avg fused")]):
        v = d[d[pred_col] >= 0]
        cm = confusion_matrix(v["y_true"], v[pred_col], labels=[0, 1, 2])
        ba = balanced_accuracy_score(v["y_true"], v[pred_col])
        im = ax.imshow(cm, cmap="Blues")
        ax.set_xticks(range(3)); ax.set_xticklabels([NAMES[i] for i in range(3)])
        ax.set_yticks(range(3)); ax.set_yticklabels([NAMES[i] for i in range(3)])
        ax.set_xlabel("predicted"); ax.set_ylabel("true")
        for i in range(3):
            for j in range(3):
                ax.text(j, i, int(cm[i, j]), ha="center", va="center", fontsize=11,
                        color="white" if cm[i, j] > cm.max() * 0.6 else "black")
        ax.set_title(f"{name}  (pooled bACC={ba:.3f})", fontsize=10)
    fig.suptitle("T4 cumulative fusion — confusion matrices (true × predicted, pooled over seeds; "
                 "n≈31)", fontsize=11)
    out = HERE / "confusion_T4.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.1); plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    plot_weights()
    write_failure_table()
    plot_confusion()
