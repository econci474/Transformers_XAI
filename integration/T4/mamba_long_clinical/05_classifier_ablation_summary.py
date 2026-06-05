r"""
05_classifier_ablation_summary.py  (env: clinical or survml — only needs pandas/sklearn/matplotlib)
===================================================================================================
Summarise the T4 3-class horizon CLASSIFIER ablation (04_train_classifier_mamba.py output): per config,
balanced accuracy / macro-F1 / accuracy on the held-out TEST converters (and val), reported BOTH pooled
across the 3 seeds (n=45 test) AND as seed mean +/- std. Renders a standalone table.

This is a SEPARATE table on a DIFFERENT cohort from the survival-derived horizon
(02_survival_comparison.py's derived_horizon_from_survival.csv): here the model is a dedicated weighted-CE
classifier head trained on the T4-aware leakage-safe folds; there it was read off the survival curve.

Run:  python integration/T4/mamba_long_clinical/05_classifier_ablation_summary.py \
        --pred_root D:/ADNI_BIDS_project/derivatives/mamba_classifier/classifier
Out:  <pred_root>/summary/{classifier_metrics.csv, classifier_metrics_summary.csv,
      classifier_ablation_table.png/pdf}
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score, f1_score, accuracy_score

HERE = Path(__file__).resolve().parent
PRED_ROOT = HERE / "outputs" / "classifier"
SEEDS = (0, 1, 2)
LABEL = {"A_default_mamba1_frozen": "mamba1 frozen (default)", "A_q2_finetune": "mamba1 finetune",
         "A_q3_time_prepend": "mamba1 frozen, time prepend", "A_q4_align": "mamba1 frozen + SigLIP align",
         "A_ctrl_meanpool": "mean-pool control", "A_ctrl_gru": "GRU control"}
ORDER = list(LABEL)


def _metrics(y, p):
    return dict(bacc=balanced_accuracy_score(y, p), macro_f1=f1_score(y, p, average="macro", zero_division=0),
                acc=accuracy_score(y, p))


def per_config(name):
    rows, pooled = [], {sp: ([], []) for sp in ("val", "test")}
    for seed in SEEDS:
        f = PRED_ROOT / name / f"seed_{seed}" / "predictions.parquet"
        if not f.exists():
            return None, None
        df = pd.read_parquet(f)
        for sp in ("val", "test"):
            g = df[df.split == sp]
            if len(g) == 0:
                continue
            y, p = g["y_true"].to_numpy(int), g["pred"].to_numpy(int)
            m = _metrics(y, p)
            rows.append({"model": name, "seed": seed, "split": sp, "n": len(g), **m})
            pooled[sp][0].append(y); pooled[sp][1].append(p)
    prows = []
    for sp in ("val", "test"):
        if pooled[sp][0]:
            y = np.concatenate(pooled[sp][0]); p = np.concatenate(pooled[sp][1])
            prows.append({"model": name, "split": sp, "n": len(y), **_metrics(y, p)})
    return rows, prows


def render(summary, pooled, out):
    """summary: per-(model,split) mean/std DataFrame ; pooled: per-(model,split) pooled DataFrame."""
    test = pooled[pooled.split == "test"].set_index("model")
    ssd = summary.xs("test", level="split")
    fig, ax = plt.subplots(figsize=(10.5, 0.6 + 0.45 * len(ORDER)))
    ax.axis("off")
    cols = ["Config", "bACC (pooled)", "bACC (seed mean±sd)", "macro-F1 (pooled)", "macro-F1 (mean±sd)", "acc (pooled)"]
    cell, best_b = [], test["bacc"].max()
    for name in ORDER:
        if name not in test.index:
            continue
        b = test.loc[name, "bacc"]
        cell.append([LABEL[name], f"{b:.3f}", f"{ssd.loc[name, ('bacc','mean')]:.3f}±{ssd.loc[name, ('bacc','std')]:.3f}",
                     f"{test.loc[name,'macro_f1']:.3f}",
                     f"{ssd.loc[name, ('macro_f1','mean')]:.3f}±{ssd.loc[name, ('macro_f1','std')]:.3f}",
                     f"{test.loc[name,'acc']:.3f}"])
    tb = ax.table(cellText=cell, colLabels=cols, loc="center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.5)
    for j in range(len(cols)):
        tb[0, j].set_facecolor("#404040"); tb[0, j].set_text_props(color="white", fontweight="bold")
    for i, name in enumerate([n for n in ORDER if n in test.index], start=1):
        if abs(test.loc[name, "bacc"] - best_b) < 1e-9:
            for j in range(len(cols)):
                tb[i, j].set_facecolor("#d9f0d3")
    ax.set_title("T4 3-class horizon CLASSIFIER head — backbone/training/time/align ablation\n"
                 "held-out TEST converters (T4-aware leakage-safe folds; n=45 pooled / 15 per seed); "
                 "weighted-CE; green = best pooled bACC", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "classifier_ablation_table.png", dpi=150)
    fig.savefig(out / "classifier_ablation_table.pdf"); plt.close(fig)


def main():
    global PRED_ROOT
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_root", type=str, default=str(PRED_ROOT))
    args = ap.parse_args()
    PRED_ROOT = Path(args.pred_root)
    out = PRED_ROOT / "summary"; out.mkdir(parents=True, exist_ok=True)

    configs = [p.name for p in sorted(PRED_ROOT.iterdir()) if p.is_dir() and p.name != "summary"]
    allrows, allpooled = [], []
    for name in configs:
        r, pr = per_config(name)
        if r:
            allrows += r; allpooled += pr
        else:
            print(f"  [skip] {name} (no predictions)")
    if not allrows:
        print("No predictions found."); return
    met = pd.DataFrame(allrows); met.to_csv(out / "classifier_metrics.csv", index=False)
    pooled = pd.DataFrame(allpooled); pooled.to_csv(out / "classifier_metrics_pooled.csv", index=False)
    summ = met.groupby(["model", "split"])[["bacc", "macro_f1", "acc"]].agg(["mean", "std"]).round(4)
    summ.to_csv(out / "classifier_metrics_summary.csv")

    print("\nT4 3-class CLASSIFIER head — TEST (pooled across seeds, n=45):")
    pt = pooled[pooled.split == "test"].copy()
    pt["model"] = pd.Categorical(pt["model"], [c for c in ORDER if c in pt["model"].values], ordered=True)
    print(pt.sort_values("bacc", ascending=False)[["model", "n", "bacc", "macro_f1", "acc"]].round(4).to_string(index=False))
    render(summ, pooled, out)
    print(f"\nDone -> {out}")


if __name__ == "__main__":
    main()
