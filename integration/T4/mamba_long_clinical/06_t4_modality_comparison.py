r"""
06_t4_modality_comparison.py  (env: clinical or survml — pandas/matplotlib only)
================================================================================
Standalone comparison answering "was the longitudinal-MAMBA classifier head worth it?" — it places the
MAMBA T4 3-class horizon classifier (05 output) next to the SIMPLE single-model T4 baselines we already
have, all on the 3-class converter horizon (<3y / 3-7y / >=7y):

  MAMBA longitudinal (clinical bl/m06/m12 embeddings)  -- this arm's 6 configs (04/05)
  Clinical encoder DIRECT head (ModernBERT-large / BioClinical-large full_ft on T4)  -- SAME 116/15/15 cohort
  Clinical tabular (LogReg / RBF-SVM / XGBoost on baseline features)                 -- SAME T4 split
  MRI-only (BrainDINO / BrainMVP / ViT-MAE cross-model T4)                           -- MRI cohort (indicative)

Metric = **mean balanced accuracy over the 3 seeds** (NOT pooled), + macro-F1, mean +/- sd.
NB chance = 0.333 ; test n = 15/seed (clinical+MAMBA share the cohort); MRI is a different cohort/split.

Out:  <pred_root>/summary/{t4_modality_comparison.csv, t4_modality_comparison.png/pdf}

Run:  python integration/T4/mamba_long_clinical/06_t4_modality_comparison.py \
        --pred_root D:/ADNI_BIDS_project/derivatives/mamba_classifier/classifier
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DERIV = Path(r"D:\ADNI_BIDS_project\derivatives")
REPO = Path(r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI")
ENC_T4 = DERIV / "encoder_outputs_no_cdr_post_exclusion_T4"
TAB_T4 = REPO / "clinical_pipeline" / "outputs" / "baseline_no_cdr_post_exclusion" / "results_per_seed_T4.csv"
MRI_T4 = REPO / "mri_pipeline" / "outputs" / "cross_model" / "cross_model_T3T4.csv"
MAMBA_LABEL = {"A_default_mamba1_frozen": "MAMBA: mamba1 frozen (default)", "A_q2_finetune": "MAMBA: mamba1 finetune",
               "A_q3_time_prepend": "MAMBA: mamba1 frozen, time-prepend", "A_q4_align": "MAMBA: mamba1 frozen + SigLIP",
               "A_ctrl_meanpool": "MAMBA: mean-pool control", "A_ctrl_gru": "MAMBA: GRU control"}


def _row(group, approach, bm, bs, fm, fs, cohort):
    return {"group": group, "approach": approach, "bacc_mean": round(bm, 4), "bacc_std": round(bs, 4),
            "f1_mean": round(fm, 4) if fm == fm else np.nan, "f1_std": round(fs, 4) if fs == fs else np.nan,
            "cohort": cohort}


def mamba_rows(pred_root):
    f = Path(pred_root) / "summary" / "classifier_metrics.csv"
    if not f.exists():
        print(f"  [warn] no MAMBA metrics at {f}"); return []
    m = pd.read_csv(f); m = m[m.split == "test"]
    rows = []
    for name, lab in MAMBA_LABEL.items():
        g = m[m.model == name]
        if len(g):
            rows.append(_row("MAMBA longitudinal (clinical emb)", lab, g.bacc.mean(), g.bacc.std(),
                             g.macro_f1.mean(), g.macro_f1.std(), "T4-split 116/15/15"))
    return rows


def encoder_rows():
    rows = []
    for model, lab in [("ModernBERT-large", "Clinical encoder: ModernBERT-large full-ft"),
                       ("BioClinical-ModernBERT-large", "Clinical encoder: BioClinical-large full-ft")]:
        bs, fs = [], []
        for s in range(3):
            p = ENC_T4 / model / "T4_conv_horizon" / f"seed_{s}" / "full_ft" / "metrics.json"
            if p.exists():
                tm = json.load(open(p))["test_metrics"]; bs.append(tm["balanced_acc"]); fs.append(tm["macro_f1"])
        if bs:
            rows.append(_row("Clinical encoder DIRECT head", lab, float(np.mean(bs)), float(np.std(bs)),
                             float(np.mean(fs)), float(np.std(fs)), "T4-split 116/15/15"))
    return rows


def tabular_rows():
    if not TAB_T4.exists():
        print(f"  [warn] no tabular T4 at {TAB_T4}"); return []
    d = pd.read_csv(TAB_T4); d = d[d.Split == "test"]
    rows = []
    for mdl, lab in [("SVM", "Clinical tabular: RBF-SVM"), ("LogReg", "Clinical tabular: LogReg"),
                     ("XGBoost", "Clinical tabular: XGBoost")]:
        g = d[d.Model == mdl]
        if len(g):
            rows.append(_row("Clinical tabular (baseline feats)", lab, g.BalancedAcc.mean(), g.BalancedAcc.std(),
                             g.MacroF1.mean(), g.MacroF1.std(), "T4 split, tabular"))
    return rows


def mri_rows():
    if not MRI_T4.exists():
        print(f"  [warn] no MRI T4 at {MRI_T4}"); return []
    d = pd.read_csv(MRI_T4); d = d[d.task == "T4_conv_horizon"].copy()
    d = d.sort_values("balanced_acc_mean", ascending=False)
    rows = []
    for _, r in d.head(2).iterrows():
        deg = f", {int(r.n_degenerate)}/3 degenerate" if r.n_degenerate else ""
        rows.append(_row("MRI-only (different cohort)", f"MRI: {r.model} {r.variant}/{r.augment}",
                         r.balanced_acc_mean, r.balanced_acc_std, r.f1_mean, r.f1_std, f"MRI cohort{deg}"))
    return rows


GROUP_COLOR = {"MAMBA longitudinal (clinical emb)": "#d9f0d3", "Clinical encoder DIRECT head": "#fde0ef",
               "Clinical tabular (baseline feats)": "#e6f0fb", "MRI-only (different cohort)": "#f2f2f2"}


def render(df, out):
    fig, ax = plt.subplots(figsize=(12.5, 0.7 + 0.42 * len(df))); ax.axis("off")
    cols = ["Approach", "Modality / features", "bACC (mean±sd)", "macro-F1 (mean±sd)", "cohort"]
    best = df["bacc_mean"].max()
    cell, colors = [], []
    for _, r in df.iterrows():
        f1 = "—" if r.f1_mean != r.f1_mean else f"{r.f1_mean:.3f}±{r.f1_std:.3f}"
        cell.append([r.approach, r.group.split("(")[0].strip(), f"{r.bacc_mean:.3f}±{r.bacc_std:.3f}", f1, r.cohort])
        colors.append(GROUP_COLOR.get(r.group, "white"))
    tb = ax.table(cellText=cell, colLabels=cols, loc="center", cellLoc="center",
                  colWidths=[0.30, 0.20, 0.15, 0.15, 0.20])
    tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.5)
    for j in range(len(cols)):
        tb[0, j].set_facecolor("#404040"); tb[0, j].set_text_props(color="white", fontweight="bold")
    for i, (_, r) in enumerate(df.iterrows(), start=1):
        for j in range(len(cols)):
            tb[i, j].set_facecolor(colors[i - 1])
        if abs(r.bacc_mean - best) < 1e-9:
            tb[i, 2].set_text_props(fontweight="bold")
    ax.set_title("T4 AD-conversion horizon (3-class) — longitudinal MAMBA classifier vs single-model baselines\n"
                 "mean balanced accuracy over 3 seeds (chance=0.333; test n=15/seed). "
                 "Clinical encoder+tabular share the MAMBA T4 cohort; MRI is a different cohort (indicative). "
                 "bold = best bACC.", fontsize=9.5)
    fig.tight_layout(); fig.savefig(out / "t4_modality_comparison.png", dpi=150)
    fig.savefig(out / "t4_modality_comparison.pdf"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_root", type=str,
                    default=str(DERIV / "mamba_classifier" / "classifier"))
    args = ap.parse_args()
    out = Path(args.pred_root) / "summary"; out.mkdir(parents=True, exist_ok=True)

    rows = mamba_rows(args.pred_root) + encoder_rows() + tabular_rows() + mri_rows()
    GROUP_ORDER = ["MAMBA longitudinal (clinical emb)", "Clinical encoder DIRECT head",
                   "Clinical tabular (baseline feats)", "MRI-only (different cohort)"]
    df = pd.DataFrame(rows)
    df["__g"] = df["group"].map({g: i for i, g in enumerate(GROUP_ORDER)})
    df = df.sort_values(["__g", "bacc_mean"], ascending=[True, False]).drop(columns="__g").reset_index(drop=True)
    df.to_csv(out / "t4_modality_comparison.csv", index=False)
    print(df[["group", "approach", "bacc_mean", "bacc_std", "f1_mean"]].to_string(index=False))
    render(df, out)
    print(f"\nDone -> {out}  (t4_modality_comparison.csv / .png / .pdf)")


if __name__ == "__main__":
    main()
