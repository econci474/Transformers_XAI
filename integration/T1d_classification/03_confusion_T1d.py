"""
03_confusion_T1d.py
===================
Render a 3-panel T1d (sMCI vs pMCI) confusion-matrix figure in the same style as the T2 integration
confusion figure (counts pooled across seeds, row-normalised %): clinical-only @m12, MRI-only @m12
(BrainMVP plus-original), and the best fusion (class-wise mean = convex equal-weight present-only over
CL_bl+CL_m06+CL_m12+MRI_m12). Reuses the temporal-fusion loaders in 01_fuse_temporal_T1d.py.

Out: integration/T1d_classification/outputs/confusion_T1d.png (+ .pdf)
Run: python integration/T1d_classification/03_confusion_T1d.py
"""
import importlib.util
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, balanced_accuracy_score

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t1d_fuse", os.path.join(HERE, "01_fuse_temporal_T1d.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
h = m.h

CLASS_NAMES = ["sMCI", "pMCI"]
SEEDS = (0, 1, 2)
MRI_TMPL = ("D:/ADNI_BIDS_project/derivatives/brainmvp_embeddings/T1d_binary/"
            "aug_plus_original/seed_{seed}/embeddings_seed_{seed}.csv")


def _pooled_single(prefix, present_flag):
    """Pool argmax predictions of one modality's probs over its present cohort, across seeds."""
    ys, ps = [], []
    for s in SEEDS:
        f, _ = m.build_seed_frame(s, MRI_TMPL)
        t = f[f.split == "test"]
        t = t[t[present_flag]]
        P = t[[f"{prefix}0", f"{prefix}1"]].to_numpy(float)
        ys += list(t["y_true"].astype(int)); ps += list(P.argmax(1))
    return np.array(ys), np.array(ps)


def _pooled_fusion():
    """Best fusion = class-wise mean (convex_present / equal) on BrainMVP, pooled across seeds."""
    d = pd.read_csv(os.path.join(HERE, "outputs", "class_wise_fusion", "BrainMVP_plusorig",
                                 "fused_predictions.csv"))
    d = d[(d.variant == "convex_present") & (d.missing == "equal")]
    return d["y_true"].astype(int).to_numpy(), d["pred"].astype(int).to_numpy()


def main():
    panels = [
        ("clinical-only [ModernBERT-large]", "month 12", *_pooled_single("m12", "m12_present")),
        ("MRI-only [BrainMVP plus-original]", "month 12", *_pooled_single("mri", "mri_present")),
        ("best fusion [class-wise, mean]", "CL bl+m06+m12 + MRI m12", *_pooled_fusion()),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(3.6 * 3, 4.3))
    for ax, (label, tp, y, pred) in zip(axes, panels):
        cm = confusion_matrix(y, pred, labels=[0, 1])
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        bacc = balanced_accuracy_score(y, pred) if len(np.unique(y)) > 1 else float("nan")
        n_tot = len(y); n_per_seed = round(n_tot / len(SEEDS))   # mean/seed (matches the table), not floor
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(CLASS_NAMES); ax.set_yticklabels(CLASS_NAMES)
        ax.set_xlabel("predicted"); ax.set_ylabel("true")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]}\n{pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=9, color="white" if pct[i, j] > 55 else "black")
        ax.set_title(f"{label}\n{tp}\n{n_per_seed}/seed TEST (Σ{n_tot})  bACC {bacc:.3f}", fontsize=9)

    fig.suptitle("T1d Multimodal late integration: Confusion matrix", fontsize=13, fontweight="bold",
                 y=0.99)
    fig.text(0.5, 0.925, "counts pooled across seeds, row-normalised %",
             ha="center", va="center", fontsize=9.5, fontstyle="italic")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    out = os.path.join(HERE, "outputs", "confusion_T1d.png")
    fig.savefig(out, dpi=170, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"wrote {out}")
    for label, tp, y, pred in panels:
        b = balanced_accuracy_score(y, pred)
        print(f"  {label:38s} n={len(y)}  bACC={b:.3f}")


if __name__ == "__main__":
    main()
