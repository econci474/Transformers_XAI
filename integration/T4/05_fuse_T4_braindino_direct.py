"""
05_fuse_T4_braindino_direct.py
==============================
TRUE / DIRECT 3-class T4 late fusion — a CONTRAST to the cumulative-T3 construction (01/02/04).
Fuses the DIRECT 3-class T4 models that score best on the T4 validation set:
  clinical = MBERT-L  (ModernBERT-large, full_ft, task T4_conv_horizon, prob_0/1/2)
  MRI      = BrainDINO frozen (cached-head HP winner lr1e-03_d0.1_ls0.0, prob_0/1/2 @ m12)

Purpose: show that two weak 3-class models are NOT rescued by fusion. Because BrainDINO did not export
VAL probabilities, only the untuned equal-weight (0.5/0.5) blend is possible (no val-tuned / EN / LR);
equal-weight averaging is exactly the canonical "does naive fusion help?" test.

Cohorts (mirrors the cumulative table): both-present = complete-case (BrainDINO m12 MRI present);
all-available = full clinical T4 test (clinical proceeds where MRI absent). Splits = horizon-stratified
baseline_T4 (clinical T4 anchor). Leakage-safe: clinical-test ∩ BrainDINO-train = ∅.

Out: integration/T4/summary_t4_cleaned_braindino.{csv,png,pdf}
     integration/T4/confusion_t4_braindino.{png,pdf,csv}
Run: PYTHONIOENCODING=utf-8 python integration/T4/05_fuse_T4_braindino_direct.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (balanced_accuracy_score, f1_score, roc_auc_score,
                             precision_recall_fscore_support, confusion_matrix)

from _t4_latex_lib import write_t4_latex

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = Path(__file__).resolve().parent
SEEDS = [0, 1, 2]
CLASSES = ["h_lt3", "h_3_7", "h_ge7"]
CLASS_DISP = ["<3y", "3–7y", "≥7y"]

CLIN_TMPL = ("D:/ADNI_BIDS_project/derivatives/encoder_outputs_no_cdr_post_exclusion_T4/"
             "ModernBERT-large/T4_conv_horizon/seed_{s}/full_ft/embeddings/embeddings.parquet")
BD_DIR = ("D:/ADNI_BIDS_project/derivatives/braindino_outputs/aug_none_hp_tuned/"
          "BrainDINO_vitb16_frozen_cached/T4_conv_horizon/seed_{s}/lr1e-03_d0.1_ls0.0")


# --------------------------------------------------------------------------- #
# Load + align
# --------------------------------------------------------------------------- #
def load_seed(seed):
    """Return a per-patient frame for clinical TEST patients with clinical probs (cp0..2),
    BrainDINO probs (mp0..2, NaN where MRI absent), y_true, and mri_present."""
    cl = pd.read_parquet(CLIN_TMPL.format(s=seed),
                         columns=["Patient_ID", "split", "label", "prob_0", "prob_1", "prob_2"])
    cl = cl[cl["split"] == "test"].rename(columns={"prob_0": "cp0", "prob_1": "cp1", "prob_2": "cp2",
                                                    "label": "y_true"})
    man = pd.read_csv(Path(BD_DIR.format(s=seed)) / "dataset_manifest.csv")
    pr = pd.read_csv(Path(BD_DIR.format(s=seed)) / "test_predictions.csv")
    man_t = man[man["split"] == "test"].reset_index(drop=True)
    assert len(man_t) == len(pr) and (man_t["label"].to_numpy() == pr["y_true"].to_numpy()).all(), \
        f"seed {seed}: BrainDINO manifest/predictions misaligned"
    bd = pd.DataFrame({"Patient_ID": man_t["Patient_ID"].to_numpy(),
                       "bids_ses": man_t["bids_ses"].to_numpy(),
                       "mp0": pr["prob_0"].to_numpy(), "mp1": pr["prob_1"].to_numpy(),
                       "mp2": pr["prob_2"].to_numpy()})
    # one MRI row per patient: BrainDINO (long_mode=all) can have bl+m12; prefer m12 (MRI@m12 convention)
    bd["_m12"] = (bd["bids_ses"] == "m12").astype(int)
    bd = (bd.sort_values("_m12", ascending=False, kind="stable")
            .drop_duplicates("Patient_ID", keep="first")[["Patient_ID", "mp0", "mp1", "mp2"]])
    f = cl.merge(bd, on="Patient_ID", how="left")
    f["mri_present"] = f[["mp0", "mp1", "mp2"]].notna().all(axis=1)
    f["y_true"] = f["y_true"].astype(int)
    return f


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def safe_auc(y, P):
    try:
        if len(np.unique(y)) < 2:
            return np.nan
        return float(roc_auc_score(y, P, multi_class="ovr", average="macro", labels=[0, 1, 2]))
    except Exception:
        return np.nan


def metric_row(y, pred, P):
    _, _, f1c, _ = precision_recall_fscore_support(y, pred, labels=[0, 1, 2], zero_division=0)
    row = {"n": len(y), "bacc": balanced_accuracy_score(y, pred),
           "macro_auc": safe_auc(y, P), "macro_f1": f1_score(y, pred, labels=[0, 1, 2],
                                                              average="macro", zero_division=0)}
    for c, f in zip(CLASSES, f1c):
        row[f"f1_{c}"] = float(f)
    return row


VARIANTS = ["clinical_only", "mri_only", "equal_0.5"]
COHORTS = [("both_present", "both-present"), ("present_only", "all-available")]
VAR_LABEL = {"clinical_only": "clinical-only [MBERT-L]",
             "mri_only": "MRI-only [BrainDINO frozen]",
             "equal_0.5": "weighted-avg (equal)"}
REST_METRICS = [("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1"),
                ("f1_h_lt3", "F1 <3y"), ("f1_h_3_7", "F1 3–7y"), ("f1_h_ge7", "F1 ≥7y")]


def preds_for(f, variant, cohort):
    """Return (y_true, pred, P) for a (variant, cohort) on seed-frame f, or None if N/A."""
    sub = f if cohort == "present_only" else f[f["mri_present"]]
    if len(sub) == 0:
        return None
    y = sub["y_true"].to_numpy(int)
    cp = sub[["cp0", "cp1", "cp2"]].to_numpy(float)
    mp = sub[["mp0", "mp1", "mp2"]].to_numpy(float)
    pres = sub["mri_present"].to_numpy()
    if variant == "clinical_only":
        P = cp
    elif variant == "mri_only":
        if cohort == "present_only":
            return None                      # MRI-only undefined on the full (MRI-absent) cohort
        P = mp
    else:                                    # equal_0.5: blend where present, else clinical
        P = cp.copy()
        if pres.any():
            P[pres] = 0.5 * cp[pres] + 0.5 * mp[pres]
    return y, P.argmax(1), P


def compute():
    frames = {s: load_seed(s) for s in SEEDS}
    per_seed, pooled_yp, breakdown = [], {}, {}
    for ckey, _ in COHORTS:
        for v in VARIANTS:
            ys, ps = [], []
            for s in SEEDS:
                r = preds_for(frames[s], v, ckey)
                if r is None:
                    continue
                y, pred, P = r
                m = metric_row(y, pred, P)
                m.update({"seed": s, "variant": v, "cohort": ckey})
                per_seed.append(m)
                ys += list(y); ps += list(pred)
            if ys:
                pooled_yp[(v, ckey)] = (np.array(ys), np.array(ps))
        # per-seed class breakdown (clinical_only = full coverage for the cohort)
        parts = []
        for s in SEEDS:
            sub = frames[s] if ckey == "present_only" else frames[s][frames[s]["mri_present"]]
            vc = sub["y_true"].value_counts().to_dict()
            parts.append(f"s{s} {vc.get(0, 0)}:{vc.get(1, 0)}:{vc.get(2, 0)}")
        breakdown[ckey] = ", ".join(parts)
    return pd.DataFrame(per_seed), pooled_yp, breakdown, frames


# --------------------------------------------------------------------------- #
# House-style table (mirrors 04_summary_t4_cleaned.py::render)
# --------------------------------------------------------------------------- #
TITLE1 = "T4 Conversion DIRECT Late Fusion (<3y / 3–7y / ≥7y)"
TITLE2 = "direct 3-class T4 models — clinical MBERT-L ⊕ MRI BrainDINO (frozen), equal-weight only"


def render_table(df, subtitle, out_path, show_n=True):
    headers = ["Cohort", "Method", "TEST bACC", "pooled bACC"] + [lbl for _, lbl in REST_METRICS] \
        + (["n"] if show_n else [])
    LEFT_COLS = {0, 1}
    best = {}
    for ckey, g in df.groupby("ckey", sort=False):
        gf = g[g.variant != "mri_only"]
        for key in ("bacc", "pooled", "macro_auc", "macro_f1"):
            best[(ckey, key)] = gf[key].idxmax() if gf[key].notna().any() else -1

    body, rules, bolds, show_c = [], [], [], []
    prev_c = None
    METRIC_KEYS = ["bacc", "pooled", "macro_auc", "macro_f1"] + [mk for mk, _ in REST_METRICS[2:]]
    for idx, r in df.iterrows():
        sc = r.ckey != prev_c
        cells = [r.Cohort if sc else "", r.Method,
                 f"{r.bacc:.3f} ± {r.bstd:.3f}",
                 f"{r.pooled:.3f}" if pd.notna(r.pooled) else "—"]
        for mk, _ in REST_METRICS:
            cells.append(f"{r[mk]:.3f}" if pd.notna(r[mk]) else "—")
        if show_n:
            cells.append(f"{r.n:.0f}")
        body.append(cells)
        bcols = {}
        col_to_key = {2: "bacc", 3: "pooled", 4: "macro_auc", 5: "macro_f1"}
        for cj, key in col_to_key.items():
            bcols[cj] = (best.get((r.ckey, key), -1) == idx)
        bolds.append(bcols)
        rules.append(prev_c is not None and r.ckey != prev_c)
        show_c.append(sc)
        prev_c = r.ckey

    COL_W = [1.35, 3.05, 1.62, 1.20, 1.15, 1.10, 0.95, 1.00, 0.95] + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_LINES = subtitle.count("\n") + 1
    SUB_H = 0.30 * SUB_LINES
    TITLE_H, HEAD_H, ROW_H = 0.72 + SUB_H, 0.40, 0.42
    TOP_PAD, BOT_PAD = 0.12, 0.12
    n_rows = len(body)
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * n_rows + BOT_PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(COL_W))]

    def hline(y, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD
    y_title_top = y
    y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) * 0.70, TITLE1,
            ha="center", va="center", fontsize=12.5, fontweight="bold")
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) * 0.26, TITLE2,
            ha="center", va="center", fontsize=9.5, fontweight="bold")
    ax.text(cx, y + SUB_H / 2, subtitle, ha="center", va="center",
            fontsize=7.8, fontstyle="italic", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=9.0, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.0, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                w = "bold" if (j == 0 and show_c[i]) else "normal"
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center",
                        fontsize=8.6, fontweight=w)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=8.6,
                        fontweight="bold" if bolds[i].get(j, False) else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


# --------------------------------------------------------------------------- #
# 3-panel confusion (both-present cohort; dual mean+pooled bACC)
# --------------------------------------------------------------------------- #
def render_confusion(frames, out_path):
    panels_def = [("Clinical-only [MBERT-L]", "ModernBERT-large · CL probs", "clinical_only"),
                  ("MRI-only [BrainDINO frozen]", "BrainDINO frozen · MRI probs", "mri_only"),
                  ("Best fusion [equal 0.5]", "0.5·CL ⊕ 0.5·MRI", "equal_0.5")]
    panels = []
    for label, tp, variant in panels_def:
        ys, ps, per = [], [], []
        for s in SEEDS:
            r = preds_for(frames[s], variant, "both_present")
            if r is None:
                continue
            y, pred, _ = r
            ys += list(y); ps += list(pred)
            per.append(balanced_accuracy_score(y, pred) if len(np.unique(y)) > 1 else np.nan)
        panels.append((label, tp, np.array(ys), np.array(ps),
                       float(np.nanmean(per)), float(np.nanstd(per, ddof=1))))

    fig, axes = plt.subplots(1, 3, figsize=(3.7 * 3, 4.7))
    csv_rows = []
    for ax, (label, tp, y, pred, mbacc, sbacc) in zip(axes, panels):
        cm = confusion_matrix(y, pred, labels=[0, 1, 2])
        row = cm.sum(1, keepdims=True)
        pct = np.divide(cm, np.where(row == 0, 1, row)) * 100
        pooled_bacc = balanced_accuracy_score(y, pred) if len(np.unique(y)) > 1 else float("nan")
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
                                 "bACC_std": round(sbacc, 4), "bACC_pooled": round(float(pooled_bacc), 4)})
        ax.set_title(f"{label}\n{tp}\n{n_ps}/seed TEST (Σ{n_tot})\n"
                     f"bACC: mean {mbacc:.3f} · pooled {pooled_bacc:.3f}", fontsize=9)
    fig.suptitle("T3d (direct) multimodal late integration: Confusion matrix",
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
    for label, tp, y, pred, mbacc, sbacc in panels:
        print(f"  {label:30s} n={len(y)} bACC mean={mbacc:.3f} pooled={balanced_accuracy_score(y, pred):.3f}")


def main():
    ps_df, pooled_yp, bd, frames = compute()
    # aggregate per (cohort, variant): mean ± std + pooled
    recs = []
    for ckey, clabel in COHORTS:
        for v in VARIANTS:
            g = ps_df[(ps_df.cohort == ckey) & (ps_df.variant == v)]
            if g.empty:
                continue
            yp = pooled_yp.get((v, ckey))
            pooled = balanced_accuracy_score(*yp) if yp is not None else np.nan
            rec = {"ckey": ckey, "Cohort": clabel, "variant": v, "Method": VAR_LABEL[v],
                   "bacc": g.bacc.mean(), "bstd": g.bacc.std(ddof=1), "pooled": pooled, "n": g.n.mean()}
            for mk, _ in REST_METRICS:
                rec[mk] = g[mk].mean()
            recs.append(rec)
    df = pd.DataFrame(recs)

    subtitle = (
        "best single T4 models by validation bACC: clinical MBERT-L (full fine-tune) · MRI BrainDINO "
        "(frozen), fused by unweighted 0.5/0.5 averaging\n"
        "purpose: contrast with the cumulative-T3 construction, two weak 3-class models are not "
        "rescued by fusion\n"
        "both-present = complete-case (BrainDINO m12 MRI present)  ·  all-available = full clinical T4 "
        "test (clinical proceeds where MRI absent)\n"
        "bACC: per-seed mean is high-variance (n ≈ 11–15/seed)  ·  pooled (all test patients) is a more "
        "trustworthy metric  ·  chance = 0.333\n"
        "N(test) per seed on splits stratified by conversion horizon (baseline_T4).  "
        f"both-present: {bd['both_present']}  ·  all-available: {bd['present_only']}   (<3y : 3–7y : ≥7y)\n"
        "no val-tuned / elastic-net / logistic-reg variants — BrainDINO did not export val probabilities")

    out_cols = ["Cohort", "Method", "bacc", "bstd", "pooled"] + [mk for mk, _ in REST_METRICS] + ["n"]
    df[out_cols].rename(columns={"bacc": "TEST_bACC_mean", "bstd": "TEST_bACC_std",
                                 "pooled": "pooled_bACC"}).to_csv(
        HERE / "summary_t4_cleaned_braindino.csv", index=False, encoding="utf-8")
    render_table(df, subtitle, HERE / "summary_t4_cleaned_braindino.png", show_n=True)
    render_table(df, subtitle, HERE / "summary_t4_cleaned_braindino_no_n.png", show_n=False)
    write_t4_latex(df, HERE / "summary_t4_cleaned_braindino.tex", TITLE1,
                   TITLE2 + " " + " ".join(subtitle.split("\n")),
                   "MBERT-L", "BrainDINO (frozen)", "tab:t4_direct_fusion")
    render_confusion(frames, HERE / "confusion_t4_braindino.png")
    print(f"  CSV: summary_t4_cleaned_braindino.csv  ({len(df)} rows)")
    print("  N(test) — both-present:", bd["both_present"], "| all-available:", bd["present_only"])
    print(df[["Cohort", "Method", "bacc", "pooled", "n"]].to_string(index=False))


if __name__ == "__main__":
    main()
