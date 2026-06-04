"""
02_summary.py — T3 per-horizon fusion summary (clinical ⊕ MRI)
==============================================================
Reads the engine `fusion_metrics.csv` written by 01_fuse_T3.py under
integration/T3/<horizon>/<arch>/ and renders one ranked table: per horizon,
the clinical-only & MRI-only references plus each MRI arch's fusion variants,
with TEST bACC / macro-AUC / macro-F1 / per-class F1 (mean ± std across seeds).

Run:  conda run -n snp python integration/T3/02_summary.py
Outputs: integration/T3/SUMMARY_T3_fusion.{png,csv}
"""

from __future__ import annotations

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
HORIZON_LABEL = {"3yr": "T3 ≤3y", "5yr": "T3 ≤5y", "7yr": "T3 ≤7y"}
HORIZON_CLIN = {"3yr": "BioClin-L", "5yr": "BioClin-B", "7yr": "ModernBERT-L"}
ARCH_LABEL = {"vit_mae": "ViT-MAE frozen/none", "brainmvp": "BrainMVP frozen/none",
              "brainmvp_ft": "BrainMVP full_ft plus_orig"}
VARIANTS = ["clinical_only", "mri_only", "weighted_avg", "weighted_avg_J", "elastic_net", "lr"]
VAR_LABEL = {"clinical_only": "clinical-only (ref)", "mri_only": "MRI-only (ref)",
             "weighted_avg": "weighted-avg (bACC)", "weighted_avg_J": "weighted-avg (Youden)",
             "elastic_net": "elastic-net", "lr": "logistic-reg"}
METRICS = [("bacc", "TEST bACC"), ("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1"),
           ("f1_conv", "F1 conv")]


def _cell(m, s):
    if pd.isna(m):
        return "—"
    return f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}"


def collect():
    rows = []
    for f in glob.glob(str(HERE / "*" / "*" / "fusion_metrics.csv")):
        horizon = Path(f).parts[-3]
        arch = Path(f).parts[-2]
        if horizon not in HORIZON_LABEL:
            continue
        d = pd.read_csv(f)
        d = d[d["framing"] == "M12"]
        # EN writes 3 identical missing-modes; keep one row per variant.
        d = d.sort_values("missing").drop_duplicates("variant", keep="first")
        for _, r in d.iterrows():
            if r["variant"] not in VARIANTS:
                continue
            rows.append(dict(horizon=horizon, arch=arch, variant=r["variant"],
                             **{f"{mk}_mean": r.get(f"{mk}_mean") for mk, _ in METRICS},
                             **{f"{mk}_std": r.get(f"{mk}_std") for mk, _ in METRICS},
                             n=r.get("n_test_mean")))
    return pd.DataFrame(rows)


def main():
    df = collect()
    if df.empty:
        raise SystemExit("No fusion_metrics.csv under integration/T3/*/* — run 01_fuse_T3.py first.")
    df.to_csv(HERE / "SUMMARY_T3_fusion.csv", index=False)

    # Build display rows grouped by horizon: clinical_only once, then per arch.
    body, kinds, first_of_h = [], [], {}
    for horizon in ["3yr", "5yr", "7yr"]:
        hd = df[df["horizon"] == horizon]
        if hd.empty:
            continue
        # clinical-only (identical across archs) — take first
        cl = hd[hd["variant"] == "clinical_only"]
        if len(cl):
            first_of_h[horizon] = len(body)
            r = cl.iloc[0]
            body.append([HORIZON_LABEL[horizon], HORIZON_CLIN[horizon], "—", VAR_LABEL["clinical_only"]]
                        + [_cell(r[f"{mk}_mean"], r[f"{mk}_std"]) for mk, _ in METRICS] + [f"{r['n']:.0f}"])
            kinds.append("ref")
        for arch in ["vit_mae", "brainmvp", "brainmvp_ft"]:
            ad = hd[hd["arch"] == arch]
            for v in ["mri_only", "weighted_avg", "weighted_avg_J", "elastic_net", "lr"]:
                rr = ad[ad["variant"] == v]
                if not len(rr):
                    continue
                r = rr.iloc[0]
                body.append(["", "—" if v == "mri_only" else HORIZON_CLIN[horizon],
                             ARCH_LABEL[arch], VAR_LABEL[v]]
                            + [_cell(r[f"{mk}_mean"], r[f"{mk}_std"]) for mk, _ in METRICS] + [f"{r['n']:.0f}"])
                kinds.append("ref" if v == "mri_only" else "fusion")

    col_labels = ["Horizon", "Clinical", "MRI", "Method"] + [lbl for _, lbl in METRICS] + ["n"]
    n_cols = len(col_labels)
    col_chars = [max(len(str(x)) for x in [col_labels[j]] + [b[j] for b in body]) + 2 for j in range(n_cols)]
    for j in (0, 1, 2, 3):
        col_chars[j] += 4
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(0.11 * total, 0.6 + 0.4 * (len(body) + 1)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=col_labels, loc="upper center",
                   cellLoc="center", colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.45)
    for j in range(n_cols):
        tab[0, j].set_text_props(weight="bold")
    for i in range(1, len(body) + 1):
        for j in (0, 1, 2, 3):
            tab[i, j].set_text_props(ha="left")
    plt.title("T3 per-horizon late fusion (clinical@bl ⊕ MRI@m12) — mean ± std across seeds\n"
              "clinical-only / MRI-only are references; fusion = weighted-avg / EN / LR",
              pad=10, fontsize=11)
    foot = ["Framing M12: clinical baseline probs + MRI m12 probs; conversion target Label_<Ny>.",
            "Clinical encoders (full_ft): ≤3y BioClinical-L, ≤5y BioClinical-B, ≤7y ModernBERT-L.",
            "MRI cached-head (frozen encoder, augment=none): ViT-MAE / BrainMVP. n = mean test cohort across seeds.",
            "weighted-avg tuned on VAL by bACC / Youden's J; EN/LR stack [clinical,MRI] probs (fit VAL).",
            "NB clinical-only/MRI-only here are on the BOTH-PRESENT m12 intersection cohort (n above), so they",
            "differ slightly from each model's standalone full-test bACC (e.g. ViT-MAE T3a 0.601, BioClin-L ~0.78)."]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.4)
    png = HERE / "SUMMARY_T3_fusion.png"
    fig.savefig(png, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"wrote {png}\nwrote {HERE/'SUMMARY_T3_fusion.csv'}")
    # console
    print("\n", df.sort_values(["horizon", "bacc_mean"], ascending=[True, False])
          [["horizon", "arch", "variant", "bacc_mean", "macro_auc_mean", "n"]].to_string(index=False))


if __name__ == "__main__":
    main()
