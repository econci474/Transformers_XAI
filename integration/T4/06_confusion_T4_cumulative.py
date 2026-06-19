"""
06_confusion_T4_cumulative.py
=============================
House-style 3-panel confusion figure (T1e/T3/T4-direct style: counts pooled across seeds, row-normalised
%, dual per-seed-mean + seed-pooled bACC) for the CUMULATIVE-T3 T4 fusion (built by 01_fuse_T4_cumulative
/ summarised by 04_summary_t4_cleaned). Supersedes the old 2-panel counts-only confusion_T4.png.

Panels (complete-case = both-present cohort, same patients):
  clinical-only [BioClin-L ≤3y + MBERT-L ≤7y]   MRI-only [BrainMVP-ft (plus-orig)]   best fusion [weighted-avg (val bACC)]
3-class horizon <3y / 3-7y / ≥7y. Reads the existing per-patient predictions (no re-fit).

Out: integration/T4/confusion_t4_cumulative.{png,pdf,csv}
Run: PYTHONIOENCODING=utf-8 python integration/T4/06_confusion_T4_cumulative.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, balanced_accuracy_score

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = Path(__file__).resolve().parent
SEEDS = [0, 1, 2]
CLASS_DISP = ["<3y", "3–7y", "≥7y"]
COHORT = "both_present"

# (variant, label, sub-line) — best fusion = weighted_avg (val bACC)
PANELS = [
    ("clinical_only", "Clinical-only [BioClin-L ≤3y + MBERT-L ≤7y]", "cumulative T3 buckets @bl"),
    ("mri_only", "MRI-only [BrainMVP-ft (plus-orig)]", "cumulative T3 buckets @m12"),
    ("weighted_avg", "Best fusion [weighted-avg (val bACC)]", "CL @bl ⊕ MRI @m12"),
]


def _bacc(y, p):
    return balanced_accuracy_score(y, p) if len(np.unique(y)) > 1 else float("nan")


def panels():
    ap = pd.read_csv(HERE / "fusion_predictions_all_T4.csv")
    ap = ap[(ap.cohort == COHORT) & (ap.pred >= 0)]
    out = []
    for variant, label, tp in PANELS:
        g = ap[ap.variant == variant]
        ys, ps, per = [], [], []
        for s in SEEDS:
            gs = g[g.seed == s]
            if gs.empty:
                continue
            y = gs.y_true.to_numpy(int); p = gs.pred.to_numpy(int)
            ys += list(y); ps += list(p); per.append(_bacc(y, p))
        out.append((label, tp, np.array(ys), np.array(ps),
                    float(np.nanmean(per)), float(np.nanstd(per, ddof=1))))
    return out


def render(pn, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(3.7 * 3, 4.7))
    csv_rows = []
    for ax, (label, tp, y, pred, mbacc, sbacc) in zip(axes, pn):
        cm = confusion_matrix(y, pred, labels=[0, 1, 2])
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        pooled = _bacc(y, pred)
        n_tot = len(y); n_ps = round(n_tot / len(SEEDS))
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
        ax.set_xticklabels(CLASS_DISP); ax.set_yticklabels(CLASS_DISP)
        ax.tick_params(axis="both", length=0, pad=2)
        ax.set_xlabel("predicted", labelpad=3); ax.set_ylabel("true", labelpad=1)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{cm[i, j]}\n{pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=8.5, color="white" if pct[i, j] > 55 else "black")
                csv_rows.append({"panel": label, "true": CLASS_DISP[i], "pred": CLASS_DISP[j],
                                 "count": int(cm[i, j]), "row_pct": round(float(pct[i, j]), 1),
                                 "n_total": n_tot, "bACC_mean": round(mbacc, 4),
                                 "bACC_std": round(sbacc, 4), "bACC_pooled": round(float(pooled), 4)})
        ax.set_title(f"{label}\n{tp}\n{n_ps}/seed TEST (Σ{n_tot})\n"
                     f"bACC: mean {mbacc:.3f} · pooled {pooled:.3f}", fontsize=9)
    fig.suptitle("T3d (from class-wise T3a/c) multimodal late integration: Confusion matrix",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.text(0.5, 0.93, "complete-case (both modalities present)  ·  counts pooled across seeds, "
             "row-normalised %  ·  bACC reported as per-seed mean and seed-pooled",
             ha="center", va="center", fontsize=8.5, fontstyle="italic")
    fig.tight_layout(rect=(0, 0, 1, 0.85))
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    pd.DataFrame(csv_rows).to_csv(out_path.with_suffix(".csv"), index=False, encoding="utf-8")
    print(f"  PNG: {out_path}")
    for label, tp, y, pred, mbacc, sbacc in pn:
        print(f"  {label:44s} n={len(y)} bACC mean={mbacc:.3f} pooled={_bacc(y, pred):.3f}")


def main():
    render(panels(), HERE / "confusion_t4_cumulative.png")


if __name__ == "__main__":
    main()
