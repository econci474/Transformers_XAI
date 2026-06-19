"""
05_summary_T1e_3way.py
======================
Single summarised house-style table for the T1e (sCN vs pCN) 3-way late fusion (clinical ⊕ SNP ⊕ MRI),
highlighting the 3-way integration: the best combination per 3-way method (avg3 / best-weights3 /
stack-EN3) on top, then the three single-modality references (best CL / SNP / MRI). Reads the existing
aggregated 3-way CSVs (no re-fit); leaves the two 28-row leaderboards (02_fuse_T1e_3way.py) untouched.

present-case == complete-case here (N(test)=22, every test subject has all three modalities), so one
table suffices (cohort = complete).

Out: integration/T1e_classification/outputs/3way/summary_t1e_3way.{csv,png,pdf} (+ _no_n)
Run: PYTHONIOENCODING=utf-8 python integration/T1e_classification/05_summary_T1e_3way.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs", "3way")
COHORT = "complete"

CL_ABBR = {"modernbert_base_full_ft": "MBERT-B-ft", "bioclin_mb_base_frozen": "BioClin-B-frozen"}
SNP_ABBR = {"meta_prs": "meta-PRS (EN)", "kosteridis_prs": "Kosteridis PRS",
            "bmfm_snp_attn1_mlp": "BMFM-attn1", "bmfm_snp_none_xgb": "BMFM-XGB"}
MRI_ABBR = {"vit_mae": "ViT-MAE (frozen)", "braindino": "BrainDINO (frozen)"}
METHOD_LABEL = {"equal_avg3": "avg3 (equal)", "best_weights3": "best-weights3 (val bACC)",
                "stacked_en3": "stack-EN3 (elastic-net)",
                "clinical_only": "CL-only", "snp_only": "SNP-only", "mri_only": "MRI-only"}
THREE_WAY = ["equal_avg3", "best_weights3", "stacked_en3"]
REFS = ["clinical_only", "snp_only", "mri_only"]
CLEAN_METRICS = [("bacc", "Test bACC"), ("macro_f1", "macro-F1"), ("auc", "AUC")]


def _weights():
    per = pd.read_csv(os.path.join(OUT_DIR, "fusion_metrics_per_seed.csv"))
    per = per[(per.split == "test") & (per.cohort == COHORT)]
    g = ["clin_variant", "snp_variant", "mri_variant", "method"]
    return per.groupby(g)[["w_clin", "w_snp", "w_mri"]].mean().round(2).reset_index()


def assemble():
    fm = pd.read_csv(os.path.join(OUT_DIR, "fusion_metrics.csv"))
    t = fm[(fm.split == "test") & (fm.cohort == COHORT)].copy()
    wm = _weights()
    t = t.merge(wm, on=["clin_variant", "snp_variant", "mri_variant", "method"], how="left")

    recs = []
    for blk, methods in (("3-way", THREE_WAY), ("single-modality", REFS)):
        for m in methods:
            sub = t[t.method == m]
            if sub.empty:
                continue
            r = sub.sort_values("bacc_mean", ascending=False, kind="stable").iloc[0]
            ref = m in REFS
            cl = "—" if m in ("snp_only", "mri_only") else CL_ABBR.get(r.clin_variant, r.clin_variant)
            snp = "—" if m in ("clinical_only", "mri_only") else SNP_ABBR.get(r.snp_variant, r.snp_variant)
            mri = "—" if m in ("clinical_only", "snp_only") else MRI_ABBR.get(r.mri_variant, r.mri_variant)
            if ref:
                w = "—"
            elif m == "stacked_en3":
                w = "EN"
            elif any(pd.isna(r.get(k)) for k in ("w_clin", "w_snp", "w_mri")):
                w = "—"
            else:
                w = f"{r.w_clin:.2f}/{r.w_snp:.2f}/{r.w_mri:.2f}"
            rec = {"block": blk, "Integration": blk, "Method": METHOD_LABEL[m],
                   "CL": cl, "SNP": snp, "MRI": mri, "weights": w, "n": r.n_mean}
            for key, _ in CLEAN_METRICS:
                rec[f"{key}_m"] = r[f"{key}_mean"]; rec[f"{key}_s"] = r[f"{key}_std"]
            recs.append(rec)
    df = pd.DataFrame(recs)
    df["rule_above"] = False
    prev = None
    for i in df.index:
        if prev is not None and df.at[i, "block"] != prev:
            df.at[i, "rule_above"] = True
        prev = df.at[i, "block"]
    return df


TITLE = "T1e 3-way Late Fusion (sCN vs pCN)"
SUBTITLE = (
    "clinical ⊕ SNP ⊕ MRI  ·  best combination per method vs single-modality references  ·  "
    "mean ± std across seeds 0/1/2\n"
    "CL = MBERT-B-ft / BioClin-B-frozen  ·  SNP = meta-PRS (EN) / Kosteridis / BMFM  ·  "
    "MRI = ViT-MAE / BrainDINO (frozen)\n"
    "weights (cl/snp/mri) tuned on VAL bACC; avg3 fixed 1/3 each; stack-EN3 = elastic-net over the 6 prob columns\n"
    "N(test) = 22 per seed (present-case = complete-case: every test subject has all three modalities)")


def render(df, out_path, show_n=True):
    headers = (["Integration", "Method", "CL", "SNP", "MRI"]
               + [lbl for _, lbl in CLEAN_METRICS] + ["weights (cl/snp/mri)"] + (["n"] if show_n else []))
    LEFT_COLS = {0, 1, 2, 3, 4, 8}
    METRIC_J0 = 5  # first metric column index
    nm = len(CLEAN_METRICS)
    numeric = np.array([[r[f"{k}_m"] for k, _ in CLEAN_METRICS] for _, r in df.iterrows()], float)
    best = [int(np.nanargmax(numeric[:, j])) for j in range(nm)]

    body, rules, show_g = [], [], []
    prev_g = None
    for _, r in df.iterrows():
        sg = r.block != prev_g
        cells = [r.Integration if sg else "", r.Method, r.CL, r.SNP, r.MRI]
        for k, _ in CLEAN_METRICS:
            cells.append(f"{r[f'{k}_m']:.3f} ± {r[f'{k}_s']:.3f}")
        cells.append(r.weights)
        if show_n:
            cells.append(f"{r.n:.0f}")
        body.append(cells); rules.append(bool(r.rule_above)); show_g.append(sg); prev_g = r.block

    COL_W = [1.75, 2.75, 1.55, 1.65, 1.85, 1.62, 1.62, 1.62, 1.95] + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_LINES = SUBTITLE.count("\n") + 1
    SUB_H = 0.28 * SUB_LINES
    TITLE_H, HEAD_H, ROW_H = 0.58 + SUB_H, 0.40, 0.44
    TOP_PAD, BOT_PAD = 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD

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
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) * 0.55, TITLE,
            ha="center", va="center", fontsize=12.5, fontweight="bold")
    ax.text(cx, y + SUB_H / 2, SUBTITLE, ha="center", va="center",
            fontsize=7.8, fontstyle="italic", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    y_head_top = y; y -= HEAD_H
    ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=8.8, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=8.8, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            mj = j - METRIC_J0
            is_bold = (0 <= mj < nm and best[mj] == i)
            if j in LEFT_COLS:
                w = "bold" if (j == 0 and show_g[i]) else "normal"
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center",
                        fontsize=8.6, fontweight=w)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=8.6,
                        fontweight="bold" if is_bold else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def main():
    df = assemble()
    cols = ["Integration", "Method", "CL", "SNP", "MRI"] \
        + [f"{k}_m" for k, _ in CLEAN_METRICS] + [f"{k}_s" for k, _ in CLEAN_METRICS] + ["weights", "n"]
    df[cols].to_csv(os.path.join(OUT_DIR, "summary_t1e_3way.csv"), index=False, encoding="utf-8")
    render(df, os.path.join(OUT_DIR, "summary_t1e_3way.png"), show_n=True)
    render(df, os.path.join(OUT_DIR, "summary_t1e_3way_no_n.png"), show_n=False)
    print(df[["Integration", "Method", "CL", "SNP", "MRI", "bacc_m", "auc_m", "weights"]].to_string(index=False))


if __name__ == "__main__":
    main()
