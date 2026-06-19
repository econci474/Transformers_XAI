"""
04_confusion_T1e.py
===================
3-panel confusion-matrix figure (T2/T1d/T3 style: counts pooled across seeds, row-normalised %) for the
top T1e fusion combination — clinical = MBERT-B-ft (ModernBERT-base, full_ft), SNP = meta-PRS (EN),
aggregation = mean (equal 0.5/0.5 blend). Panels: clinical-only (argmax CL), SNP-only (argmax SNP),
fusion (argmax of 0.5·CL + 0.5·SNP). CL and SNP are present for all 22 test subjects per seed, so all
three panels evaluate on the same cohort. Reuses the 01 loaders (no metrics re-fit).

Each panel reports BOTH bACCs: the per-seed **mean** (matches the summary tables) and the **pooled**
bACC (computed on predictions concatenated across seeds — robust to per-seed class imbalance, tells you
whether an apparent gain holds on the aggregated cohort rather than being driven by one lucky seed).

Out: integration/T1e_classification/outputs/baseline/confusion_MBERT_metaPRS.{csv,png,pdf}
Run: python integration/T1e_classification/04_confusion_T1e.py
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
OUT_DIR = os.path.join(HERE, "outputs", "baseline")
spec = importlib.util.spec_from_file_location("t1e_fuse", os.path.join(HERE, "01_fuse_T1e.py"))
t1e = importlib.util.module_from_spec(spec); spec.loader.exec_module(t1e)

SEEDS = [0, 1, 2]
CLASS_NAMES = ["sCN", "pCN"]
CL_KEY = "modernbert_base_full_ft"          # MBERT-B-ft

# (SNP_KEY, snp_short_panel_label, snp_long_panel_label, output_basename)
CONFIG = [
    ("meta_prs",          "meta-PRS",      "meta-PRS (EN)",       "confusion_MBERT_metaPRS"),
    ("meta_prs_filtered", "meta-PRS filt", "meta-PRS (filtered)", "confusion_MBERT_metaPRS_filtered"),
]


def pooled_panels(snp_key, snp_short, snp_long):
    cl = next(v for v in t1e.CLIN_VARIANTS if v["key"] == CL_KEY)
    sn = next(v for v in t1e.SNP_VARIANTS if v["key"] == snp_key)
    pooled = {"clin": [[], []], "snp": [[], []], "fuse": [[], []]}   # each: [y_true, pred] (counts)
    per_seed = {"clin": [], "snp": [], "fuse": []}                    # per-seed bACC (to average)
    for seed in SEEDS:
        clin = t1e.load_clinical(seed, cl["model"], cl["strat"])
        snp = t1e.load_snp(seed, sn["parquet_tmpl"])
        f = clin.merge(snp[["Patient_ID", "sp0", "sp1"]], on="Patient_ID", how="inner")
        f = f[f["split"] == "test"]
        y = f["y_clin"].to_numpy(int)
        cp = f[["cp0", "cp1"]].to_numpy(float)
        sp = f[["sp0", "sp1"]].to_numpy(float)
        preds = {"clin": cp.argmax(1), "snp": sp.argmax(1), "fuse": (0.5 * cp + 0.5 * sp).argmax(1)}
        for k, p in preds.items():
            pooled[k][0] += list(y); pooled[k][1] += list(p)
            per_seed[k].append(balanced_accuracy_score(y, p) if len(np.unique(y)) > 1 else float("nan"))
    mean = {k: float(np.mean(v)) for k, v in per_seed.items()}        # mean-across-seeds (matches tables)
    std = {k: float(np.std(v, ddof=1)) for k, v in per_seed.items()}
    return [
        ("Clinical-only [MBERT-B-ft]", "ModernBERT-base (full_ft) · CL probs",
         *pooled["clin"], mean["clin"], std["clin"]),
        (f"SNP-only [{snp_short}]", f"{snp_long} · SNP probs", *pooled["snp"],
         mean["snp"], std["snp"]),
        ("Best fusion [mean]", "0.5·CL ⊕ 0.5·SNP", *pooled["fuse"], mean["fuse"], std["fuse"]),
    ]


def render(panels, basename):
    fig, axes = plt.subplots(1, 3, figsize=(3.6 * 3, 4.5))
    csv_rows = []
    for ax, (label, tp, y, pred, mbacc, sbacc) in zip(axes, panels):
        y = np.asarray(y, int); pred = np.asarray(pred, int)
        cm = confusion_matrix(y, pred, labels=[0, 1])
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        pooled_bacc = balanced_accuracy_score(y, pred) if len(np.unique(y)) > 1 else float("nan")
        n_tot = len(y); n_ps = round(n_tot / len(SEEDS))
        ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(CLASS_NAMES); ax.set_yticklabels(CLASS_NAMES)
        ax.tick_params(axis="both", length=0, pad=2)
        ax.set_xlabel("predicted", labelpad=3); ax.set_ylabel("true", labelpad=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]}\n{pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=9, color="white" if pct[i, j] > 55 else "black")
                csv_rows.append({"panel": label, "true": CLASS_NAMES[i], "pred": CLASS_NAMES[j],
                                 "count": int(cm[i, j]), "row_pct": round(float(pct[i, j]), 1),
                                 "n_total": n_tot, "bACC_mean": round(float(mbacc), 4),
                                 "bACC_std": round(float(sbacc), 4),
                                 "bACC_pooled": round(float(pooled_bacc), 4)})
        ax.set_title(f"{label}\n{tp}\n{n_ps}/seed TEST (Σ{n_tot})\n"
                     f"bACC: mean {mbacc:.2f} · pooled {pooled_bacc:.2f}", fontsize=9)
    fig.suptitle("T2b multimodal late integration: Confusion matrix",
                 fontsize=13, fontweight="bold", y=0.99)
    fig.text(0.5, 0.93, "counts pooled across seeds, row-normalised %  ·  "
             "bACC reported as per-seed mean and seed-pooled",
             ha="center", va="center", fontsize=9.5, fontstyle="italic")
    fig.tight_layout(rect=(0, 0, 1, 0.85))
    out = os.path.join(OUT_DIR, f"{basename}.png")
    fig.savefig(out, dpi=170, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    pd.DataFrame(csv_rows).to_csv(out.replace(".png", ".csv"), index=False, encoding="utf-8")
    print(f"wrote {out}  (+ .pdf, .csv)")
    for label, tp, y, pred, mbacc, sbacc in panels:
        pooled_bacc = balanced_accuracy_score(y, pred)
        print(f"  {label:30s} n={len(y)} bACC mean={mbacc:.4f} pooled={pooled_bacc:.4f}")


def main():
    for snp_key, snp_short, snp_long, basename in CONFIG:
        print(f"\n=== {basename} ({snp_key}) ===")
        render(pooled_panels(snp_key, snp_short, snp_long), basename)


if __name__ == "__main__":
    main()
