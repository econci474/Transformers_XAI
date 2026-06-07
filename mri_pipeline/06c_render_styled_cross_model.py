"""
06c_render_styled_cross_model.py — clinical-style MRI cross-model sub-tables
===========================================================================
Re-renders the MRI cross-model results (06_render_cross_model_table.py outputs)
in the SAME publication aesthetic as the clinical combined table
(clinical_pipeline/04d_plot_combined_table_post_exclusion.py): a bordered box,
grouped task headers with an N(val)/N(test) subtitle per task, dashed model-family
separators, italic metric headers, best-per-column in bold, and NO per-cell "(3)"
seed annotation (seed count stated once in the title).

Produces, under mri_pipeline/outputs/cross_model/:
  mri_t1t1bt2_val.png/.csv   — T1 / T1b / T2, VAL balanced accuracy (only metric tracked on val)
  mri_t1t1bt2_test.png/.csv  — T1 / T1b / T2, TEST Bal.Acc / AUC / F1
  mri_t1d_val.png/.csv       — T1d (pMCI vs sMCI), VAL bACC
  mri_t1d_test.png/.csv      — T1d, TEST Bal.Acc / AUC / F1
  mri_t1t1bt2_{val,test}_aug_hp.png  — the augmentation + per-row HP key, standalone
  mri_t1d_{val,test}_aug_hp.png      — same key (kept beside each table for convenience)

Data sources (already aggregated by 06):
  cross_model_table_val.csv   — val_bacc_mean/std per (model,variant,augment,task)
  cross_model_table_test.csv  — balanced_acc/auc/f1 (+std) per row

Per-task held-out N (scans) is taken from the supervised viscode2-stratified MRI
training metrics (single source the bulk of the rows share); cached-encoder rows
use the larger all-visits cohort — noted in the table footnote.

Usage:  python mri_pipeline/06c_render_styled_cross_model.py
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "outputs" / "cross_model"

# ── Per-task held-out scan counts (all-visits MRI cohort, POST-EXCLUSION). ───
#    Source: diagnostic_coverage_table_cumul_stratified_viscode2_post_exclusion.csv
#    @ window ≤m240 (= all longitudinal visits). Its MRI_val/MRI_test per-class
#    counts match the cached cross-model embeddings exactly. Classes are
#    CN / MCI / AD by PATIENT BASELINE DIAGNOSIS; T1a/T1b groupings derived.
#    Per-seed (0/1/2):  val CN 48/56/59, MCI 73/70/79, AD 13/17/16  (tot 134/143/154)
#                       test CN 56/41/61, MCI 67/72/78, AD 31/18/26 (tot 154/131/165)
#    Reported below = mean across the 3 seeds, rounded. ────────────────────────
COHORT = {
    "val":  {"CN": 54, "MCI": 74, "AD": 15, "total": 143},
    "test": {"CN": 53, "MCI": 72, "AD": 25, "total": 150},
}
# T1d (pMCI vs sMCI) is the MCI-only conversion subset (per-class split not shown).
N_T1D = {"val": 32, "test": 33}

# Conversion-task held-out N + per-class breakdown (Clinical-table style). These
# cohorts differ from the diagnostic tasks (T1e = CN-only stable-vs-progress;
# T3/T4 = AD-conversion subset). Source: dataset_manifest.csv label counts in the
# val/test split, mean across seeds 0/1/2, rounded. Each entry: an ordered list of
# (class_label, count) per split; total is the sum.
N_CONV = {
    "T1e_pcn_vs_scn":  {"val": [("sCN", 15), ("pCN", 6)],
                        "test": [("sCN", 18), ("pCN", 4)]},
    "T3a_conv3y":      {"val": [("non-conv", 35), ("AD≤3yr", 4)],
                        "test": [("non-conv", 34), ("AD≤3yr", 9)]},
    "T3b_conv5y":      {"val": [("non-conv", 23), ("AD≤5yr", 6)],
                        "test": [("non-conv", 23), ("AD≤5yr", 10)]},
    "T3c_conv7y":      {"val": [("non-conv", 17), ("AD≤7yr", 6)],
                        "test": [("non-conv", 14), ("AD≤7yr", 11)]},
    "T4_conv_horizon": {"val": [("<3yr", 6), ("3-7yr", 3), ("≥7yr", 3)],
                        "test": [("<3yr", 6), ("3-7yr", 4), ("≥7yr", 3)]},
}

TASK_LABELS = {
    "T1_binary":     "T1a: CN vs MCI+AD",
    "T1b_binary":    "T1b: CN+MCI vs AD",
    "T1d_binary":    "T1d: pMCI vs sMCI",
    "T2_multiclass": "T2: CN / MCI / AD",
    "T1e_pcn_vs_scn":  "T1e: sCN vs pCN",
    "T3a_conv3y":      "T3a: conv ≤3y",
    "T3b_conv5y":      "T3b: conv ≤5y",
    "T3c_conv7y":      "T3c: conv ≤7y",
    "T4_conv_horizon": "T4: horizon <3 / 3-7 / ≥7y",
}


def task_subtitle(task, split):
    """Two-line subtitle: (N line, per-class line) for the given task/split."""
    sp = "val" if split == "val" else "test"
    if task == "T1d_binary":
        return f"N({sp})={N_T1D[split]}", ""
    if task in N_CONV:
        per_cls = N_CONV[task][split]
        total = sum(c for _, c in per_cls)
        cls_line = ", ".join(f"{lbl}={c}" for lbl, c in per_cls)
        return f"N({sp})={total}", cls_line
    c = COHORT[split]
    n_line = f"N({sp})={c['total']}"
    if task == "T1_binary":
        cls = f"CN={c['CN']}, MCI+AD={c['MCI'] + c['AD']}"
    elif task == "T1b_binary":
        cls = f"CN+MCI={c['CN'] + c['MCI']}, AD={c['AD']}"
    elif task == "T2_multiclass":
        cls = f"CN={c['CN']}, MCI={c['MCI']}, AD={c['AD']}"
    else:
        cls = ""
    return n_line, cls

# ── Row ordering + labelling (mirrors 06_render_cross_model_table.py) ────────
AUG_SUPERSCRIPT = {"none": "⁰", "random": "¹", "stochastic": "²",
                   "plus_original": "³", "flips": "⁴", "flips+strong": "⁵"}
HP_TUNED_MARKER = "†"  # †
MODEL_RANK = {"BrainDINO": 0, "BrainMVP": 2, "ViT-MAE75": 4, "ViT-Base": 6,
              "AG-MS3D": 7, "3D-CNN": 9}
VARIANT_RANK = {"full_ft": 0, "lora": 1, "frozen": 2, "vanilla": 3,
                "separable": 4, "baseline": 5, "scratch": 5, "agms3d": 6, "-": 9}
AUG_RANK = {"plus_original": 0, "stochastic": 1, "random": 1, "none": 2, "-": 9}

# Display renames applied at load: collapse the two AG-MS3D backbones into one
# family (variant carries vanilla/separable), and rename Spasov-CNN -> 3D-CNN,
# ViT-scratch -> ViT-Base. Variant agms3d -> separable (only the AG-MS3D-sep row).
MODEL_DISPLAY = {"Spasov-CNN": "3D-CNN", "AG-MS3D-vanilla": "AG-MS3D",
                 "AG-MS3D-sep": "AG-MS3D", "ViT-scratch": "ViT-Base"}
VARIANT_DISPLAY = {"agms3d": "separable"}


def _aug_clean(a):
    a = "-" if a in (None, "", "nan") or (isinstance(a, float) and np.isnan(a)) else str(a)
    return a[:-len(HP_TUNED_MARKER)] if a.endswith(HP_TUNED_MARKER) else a


def _sort_key(t):
    m, v, a = t
    ac = _aug_clean(a)
    return (MODEL_RANK.get(m, 99), VARIANT_RANK.get(v, 99), AUG_RANK.get(ac, 99), m, v, ac)


def _varaug_label(variant, augment):
    """variant / aug with the Unicode superscript marker + HP-tuned dagger."""
    if augment in ("-", None, "") or (isinstance(augment, float) and np.isnan(augment)):
        return f"{variant}"
    is_hp = isinstance(augment, str) and augment.endswith(HP_TUNED_MARKER)
    ac = _aug_clean(augment)
    sup = AUG_SUPERSCRIPT.get(ac, "")          # superscript keyed on the recipe
    disp = "stochastic" if ac == "random" else ac   # consistent naming (recipe kept via ¹/²)
    return f"{variant} / {disp}{sup}{HP_TUNED_MARKER if is_hp else ''}"


# ── Augmentation + HP key text (verbatim from 06_render_cross_model_table.py) ─
AUG_LEGEND_LINES = [
    "Augmentation key:",
    "⁰ none          = no train-time augmentation (frozen-encoder cached forward, or --augment=none).",
    "¹ stochastic    = ViT trainer default; RandFlip×3 + RandRotate90 + RandScaleIntensity + RandShiftIntensity,",
    "                  each fires stochastically at p ≈ 0.2.",
    "² stochastic    = BrainMVP trainer; RandFlip×3 + RandAffine + RandGaussianNoise + intensity jitter,",
    "                  always-on stochastic.",
    "³ plus_original = originals retained + K=1 deterministic (p=1.0) augmented copy; inner ops match the",
    "                  trainer's stochastic set (ViT: stochastic ops; BrainMVP: stochastic ops).",
    "⁴ flips         = RandFlipd along 3 axes only (AG-MS3D legacy, AG-MS3D rescue1 vanilla, Spasov-CNN).",
    "† HP-tuned      = cached-head sweep (frozen encoder + linear head on cached embeddings; trainer =",
    "                  04_head_finetune_from_embeddings.py). Grid: lr ∈ {1e-3, 1e-4, 1e-5} × drop_rate ∈",
    "                  {0.1, 0.2, 0.3} × label_smoothing ∈ {0.0, 0.1} = 18 HP combos per (task, seed);",
    "                  epochs=50, wd=1e-5, patience=15. HP winner selected per task by mean val_bacc across",
    "                  3 seeds. Non-† rows are single-HP runs.",
]
HP_KEY_LINES = [
    "HP key (training HPs for the runs reported above; verified vs W&B run configs):",
    "  BrainDINO  / full_ft / *           lr=1e-4, wd=1e-5, drop=0.2, ls=0.0, epochs=50,  patience=25, bs=2, Adam.",
    "  BrainDINO  / lora    / *           lr=1e-4, wd=1e-5, drop=0.2, ls=0.0, epochs=8,   patience=25, bs=2, lora_rank=8, Adam.",
    "  BrainDINO  / frozen  / {stoch., plus_orig}  lr=1e-4, wd=1e-5, drop=0.2, ls=0.0, epochs=100, patience=25, bs=2, Adam (encoder frozen).",
    "  BrainDINO  / frozen  / none†       cached-head HP sweep; HP-winner picked per task from the † grid above.",
    "  BrainMVP   / full_ft / *           lr=5e-5, wd=1e-4, drop=0.1, ls=0.0, epochs=200, patience=50, bs=4, AdamW.",
    "  BrainMVP   / frozen  / {stoch., plus_orig}  lr=1e-3, wd=1e-4, drop=0.1, ls=0.0, epochs=100, patience=50, bs=4, AdamW.",
    "  BrainMVP   / frozen  / none†       cached-head HP sweep; HP-winner picked per task from the † grid above.",
    "  ViT-MAE75  / full_ft / *           lr=1e-4, wd=1e-4, drop_path=0.1, attn_drop=0.1, ls=0.0, llrd_gamma=1.0,",
    "                                      epochs=50, patience=10, bs=4 (eff bs=32 via grad_accum=8), AdamW + 10-epoch warmup.",
    "  ViT-MAE75  / frozen  / stochastic  lr=1e-3, wd=1e-4, drop_path=0.1, attn_drop=0.1, ls=0.0, llrd_gamma=0.7",
    "                                      (LLRD is a math no-op for frozen), epochs=70, patience=10, bs=4, grad_accum=8, AdamW.",
    "  ViT-MAE75  / frozen  / plus_orig.  as ViT-MAE75/frozen/stochastic but llrd_gamma=1.0.",
    "  ViT-MAE75  / frozen  / none†       cached-head HP sweep; HP-winner picked per task from the † grid above.",
    "  ViT-Base   / baseline / stochastic lr=1e-3, wd=1e-4, drop_path=0.1, attn_drop=0.1, ls=0.0, llrd_gamma=1.0,",
    "                                      epochs=100, patience=10, bs=4 (eff bs=32 via grad_accum=8), AdamW + warmup, random init.",
    "  AG-MS3D-vanilla / vanilla / flips  lr=3e-3, ls=0.0, epochs=60, bs=8, base_filters=32 (rescue1 vanilla Conv3d, RandFlipd×3).",
    "  AG-MS3D-sep     / agms3d / flips   same HPs as AG-MS3D-vanilla, separable Conv3d backbone (legacy; collapsed 14/15 cells).",
    "  Spasov-CNN / vanilla   / flips     lr=1e-3, wd=1e-4, drop=0.1, ls=0.0, epochs=60, patience=15, bs=8, warmup=5 (n_params≈518k).",
    "  Spasov-CNN / separable / flips     same HPs as Spasov-CNN/vanilla, separable Conv3d backbone (n_params≈26k).",
]

# Metric column sets per split.
VAL_METRICS = [("val_bacc", "Bal.Acc"), ("val_auc", "AUC"), ("val_f1", "F1")]
TEST_METRICS = [("balanced_acc", "Bal.Acc"), ("auc", "AUC"), ("f1", "F1")]

# ── Placeholder rows for runs that did NOT finish ────────────────────────────
# Injected (per split) at chance level and flagged with ‡ so the table is
# complete without implying a measured result. Skipped if a real row already
# exists for that (model, task). TEMPORARY — drop once the run completes.
PLACEHOLDER_ROWS = [
    # ViT-Base T1d trained on Colab (2026-06-06); real row now supersedes the
    # former chance-level placeholder. Keep the mechanism for future gaps.
]
PLACEHOLDER_FILL = {
    "val":  {"val_bacc_mean": 0.5, "val_bacc_std": 0.0},
    "test": {"balanced_acc_mean": 0.5, "balanced_acc_std": 0.0,
             "auc_mean": 0.5, "auc_std": 0.0, "f1_mean": 0.0, "f1_std": 0.0},
}
PLACEHOLDER_NOTE = ("‡ run did not finish — chance-level placeholder "
                    "(0.500), not a measured result.")


def _fmt(mean, std):
    if mean is None or (isinstance(mean, float) and np.isnan(mean)):
        return "—"
    return f"{mean:.3f} ± {std:.3f}" if std is not None and not np.isnan(std) else f"{mean:.3f}"


def _load(split):
    f = OUT_DIR / (f"cross_model_table_{split}.csv")
    df = pd.read_csv(f)
    df["augment"] = df["augment"].fillna("-")
    df["model"] = df["model"].replace(MODEL_DISPLAY)
    df["variant"] = df["variant"].replace(VARIANT_DISPLAY)
    df["is_placeholder"] = False
    # Inject placeholder rows for unfinished runs (only where no real row exists).
    ph = []
    for r in PLACEHOLDER_ROWS:
        exists = ((df["model"] == r["model"]) & (df["task"] == r["task"])).any()
        if not exists:
            ph.append({**r, **PLACEHOLDER_FILL[split], "is_placeholder": True})
    if ph:
        df = pd.concat([df, pd.DataFrame(ph)], ignore_index=True)
        df["is_placeholder"] = df["is_placeholder"].fillna(False)
    return df


def render_table(df, split, tasks, metrics, out_stem, title_lines=None):
    """Clinical-style bordered table for the given tasks/metrics.

    title_lines: optional (line1_bold, line2_italic, line3_italic) override for
    the 3-line title; if None, uses the default diagnostic-table title."""
    # rows present for these tasks (drop all-missing rows)
    metric_mean0 = metrics[0][0] + "_mean"
    present = df[df["task"].isin(tasks)].copy()
    groups = sorted({(r.model, r.variant, r.augment) for r in present.itertuples()
                     if not pd.isna(getattr(r, metric_mean0))}, key=_sort_key)
    if not groups:
        print(f"  [skip] no rows for {out_stem}"); return

    # cell lookup: (model,variant,augment,task) -> row
    idx = {(r.model, r.variant, r.augment, r.task): r for r in present.itertuples()}

    # placeholder groups (unfinished runs) — flagged ‡, excluded from bolding
    ph_keys = {(r.model, r.variant, r.augment) for r in present.itertuples()
               if getattr(r, "is_placeholder", False)}

    n_metrics = len(metrics)
    n_data_cols = len(tasks) * n_metrics

    # Build body strings + numeric matrix (for bolding) + per-task subtitles.
    body, numeric = [], []
    for (m, v, a) in groups:
        cells, nums = [], []
        for t in tasks:
            r = idx.get((m, v, a, t))
            r_is_ph = r is not None and getattr(r, "is_placeholder", False)
            for pref, _ in metrics:
                if r is None or pd.isna(getattr(r, f"{pref}_mean")):
                    cells.append("—"); nums.append(np.nan)
                else:
                    cells.append(_fmt(getattr(r, f"{pref}_mean"), getattr(r, f"{pref}_std")))
                    nums.append(np.nan if r_is_ph else float(getattr(r, f"{pref}_mean")))
        body.append(cells); numeric.append(nums)
    numeric = np.array(numeric, float)
    n_rows = len(groups)

    # best per data column (overall max, ignoring NaN)
    best_mask = np.zeros((n_rows, n_data_cols), bool)
    for j in range(n_data_cols):
        col = numeric[:, j]
        if not np.all(np.isnan(col)):
            best_mask[int(np.nanargmax(col)), j] = True

    # ── Sizing ──────────────────────────────────────────────────────────────
    # Single-metric (val) groups must be wide enough for the "N(val)=…, N(test)=…"
    # subtitle; multi-metric (test) groups are already wide from the 3 sub-columns.
    COL_W = 1.95 if n_metrics == 1 else 1.34
    ROW_H, ROW_H_SUB = 0.32, 0.46   # subtitle band holds two lines (N + per-class)
    MODEL_COL_W, VARAUG_COL_W = 1.45, 2.15
    footnotes = []                 # comments removed per request (key lives in *_aug_hp.png)
    FOOTER_H = 0.20

    title_h = ROW_H * 2.6           # three title lines
    fig_w = 0.3 + MODEL_COL_W + VARAUG_COL_W + n_data_cols * COL_W + 0.3
    fig_h = 0.15 + title_h + ROW_H * (2 + n_rows) + ROW_H_SUB + 0.10 + FOOTER_H
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)

    LEFT = 0.3
    col_widths = [MODEL_COL_W, VARAUG_COL_W] + [COL_W] * n_data_cols
    col_left = [LEFT]
    for w in col_widths:
        col_left.append(col_left[-1] + w)
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(col_widths))]
    RIGHT = col_left[-1]

    TOP = fig_h - 0.15
    row_tops = [TOP]
    y = TOP
    for h in (title_h, ROW_H, ROW_H_SUB, ROW_H):       # title, group, subtitle, header
        y -= h; row_tops.append(y)
    for _ in range(n_rows):
        y -= ROW_H; row_tops.append(y)
    total_visual_rows = 4 + n_rows

    def row_mid(i): return (row_tops[i] + row_tops[i + 1]) / 2

    def hline(yv, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yv, yv], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    hline(row_tops[0], lw=1.5)   # top of title
    hline(row_tops[1], lw=1.5)   # below title
    hline(row_tops[3], lw=0.5, ls=(0, (3, 3)))   # below subtitle band
    hline(row_tops[4], lw=1.2)   # below header row
    hline(row_tops[total_visual_rows], lw=1.5)   # bottom

    # vertical dashed separators between task groups + after the two label cols
    col_cursor = 2
    for ti in range(len(tasks)):
        x = col_left[col_cursor]
        ax.plot([x, x], [row_tops[total_visual_rows], row_tops[1]], color="black",
                linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
        col_cursor += n_metrics
    ax.plot([col_left[1], col_left[1]], [row_tops[total_visual_rows], row_tops[3]],
            color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

    # dashed family separators between model groups
    for ri in range(1, n_rows):
        if groups[ri][0] != groups[ri - 1][0]:
            hline(row_tops[4 + ri], lw=0.6, ls="--")

    # ── Title (three lines: bold headline + two regular sub-lines) ──────────
    if title_lines is not None:
        line1, line2, line3 = title_lines
    else:
        split_word = {"val": "Validation", "test": "Test"}[split]
        metric_word = ("balanced accuracy (early-stop metric)" if split == "val"
                       else "Bal.Acc / AUC / F1")
        line1 = f"MRI models: {split_word} {metric_word} by task"
        line2 = "(mean ± std across 3 seeds)"
        line3 = "Models trained using scans over all longitudinal visits, stratified by patient baseline diagnosis"
    title_fs = 7.5 if n_data_cols <= 3 else 9
    xc = (LEFT + RIGHT) / 2
    if title_lines is not None:
        # Clinical-table style: all three lines bold + centred, one block.
        ax.text(xc, row_tops[0] - title_h * 0.5,
                line1 + "\n" + line2 + "\n" + line3, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", linespacing=1.4)
    else:
        ax.text(xc, row_tops[0] - title_h * 0.30, line1, ha="center", va="center",
                fontsize=title_fs, fontweight="bold")
        ax.text(xc, row_tops[0] - title_h * 0.78, line2 + "\n" + line3, ha="center",
                va="center", fontsize=title_fs - 1.3, fontstyle="italic",
                color="#333333", linespacing=1.25)

    # ── Task group headers + N/per-class subtitle (2 lines) + underline ────
    col_cursor = 2
    for t in tasks:
        x0, x1 = col_left[col_cursor], col_left[col_cursor + n_metrics]
        xm = (x0 + x1) / 2
        ax.text(xm, row_mid(1), TASK_LABELS[t], ha="center", va="center",
                fontsize=8, fontweight="bold")
        n_line, cls_line = task_subtitle(t, split)
        sub_top, sub_bot = row_tops[2], row_tops[3]
        ax.text(xm, sub_top - (sub_top - sub_bot) * 0.34, n_line, ha="center",
                va="center", fontsize=6.4, fontstyle="italic", color="#444444")
        if cls_line:
            ax.text(xm, sub_top - (sub_top - sub_bot) * 0.70, cls_line, ha="center",
                    va="center", fontsize=5.9, fontstyle="italic", color="#666666")
        ax.plot([x0 + 0.05, x1 - 0.05], [row_tops[2], row_tops[2]], color="black",
                linewidth=0.7, zorder=4)
        col_cursor += n_metrics

    # ── Column header row: Model | variant/aug | metric headers ───────────
    ax.text(col_cx[0], row_mid(3), "Model", ha="center", va="center",
            fontsize=7.5, fontstyle="italic")
    ax.text(col_cx[1], row_mid(3), "variant / aug", ha="center", va="center",
            fontsize=7.5, fontstyle="italic")
    col_cursor = 2
    for _t in tasks:
        for _pref, h in metrics:
            ax.text(col_cx[col_cursor], row_mid(3), h, ha="center", va="center",
                    fontsize=7.3, fontstyle="italic")
            col_cursor += 1

    # ── Data rows ─────────────────────────────────────────────────────────
    for ri, (m, v, a) in enumerate(groups):
        row_i = 4 + ri
        ax.text(col_left[0] + 0.06, row_mid(row_i), m, ha="left", va="center",
                fontsize=7.3, fontweight="bold")
        ax.text(col_left[1] + 0.06, row_mid(row_i), _varaug_label(v, a), ha="left",
                va="center", fontsize=7.1)
        for ci in range(n_data_cols):
            fw = "bold" if best_mask[ri][ci] else "normal"
            ax.text(col_cx[2 + ci], row_mid(row_i), body[ri][ci], ha="center",
                    va="center", fontsize=7.0, fontweight=fw)

    ax.add_patch(plt.Rectangle((LEFT, row_tops[total_visual_rows]),
                               RIGHT - LEFT, row_tops[0] - row_tops[total_visual_rows],
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))

    y_fn = row_tops[total_visual_rows] - 0.22
    for line in footnotes:
        ax.text(LEFT, y_fn, line, ha="left", va="top", fontsize=6.3, color="#222222")
        y_fn -= 0.17

    fig.savefig(OUT_DIR / f"{out_stem}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    # CSV mirror of the rendered cells
    col_names = [f"{TASK_LABELS[t]} | {h}" for t in tasks for _, h in metrics]
    recs = []
    for ri, (m, v, a) in enumerate(groups):
        rec = {"Model": m, "variant / aug": _varaug_label(v, a)}
        rec.update(dict(zip(col_names, body[ri])))
        recs.append(rec)
    pd.DataFrame(recs).to_csv(OUT_DIR / f"{out_stem}.csv", index=False)
    print(f"  wrote {out_stem}.png/.csv  ({n_rows} rows)")


def render_aug_hp(out_stem):
    """Standalone augmentation + HP key (the legend that used to sit under the table)."""
    lines = AUG_LEGEND_LINES + [""] + HP_KEY_LINES
    fig_h = 0.22 + 0.165 * len(lines)
    fig, ax = plt.subplots(figsize=(13.0, fig_h))
    ax.axis("off")
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=8, family="monospace", linespacing=1.45)
    fig.savefig(OUT_DIR / f"{out_stem}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  wrote {out_stem}.png")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    val = _load("val")
    test = _load("test")

    TABLES = [
        ("mri_t1t1bt2", ["T1_binary", "T1b_binary", "T2_multiclass"]),
        ("mri_t1d",     ["T1d_binary"]),
    ]
    for stem, tasks in TABLES:
        render_table(val,  "val",  tasks, VAL_METRICS,  f"{stem}_val")
        render_table(test, "test", tasks, TEST_METRICS, f"{stem}_test")
        render_aug_hp(f"{stem}_val_aug_hp")
        render_aug_hp(f"{stem}_test_aug_hp")
    print("Done.")


if __name__ == "__main__":
    main()
