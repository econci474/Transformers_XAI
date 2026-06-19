"""
04_summary_t4_cleaned.py
========================
Cleaned house-style copy of the T4 conversion-horizon CUMULATIVE fusion table (the plain table is
`integration/T4/SUMMARY_T4_fusion.png`, built by 02_summary.py). Same DejaVu-Serif / bordered / ruled /
italic-header / dashed-divider / bold-best house style as the T3/T1d/T1e cleaned tables.

Two things this version makes explicit (per request):
  1. T4 is NOT modelled directly — the three horizon buckets (<3y / 3-7y / ≥7y) are CONSTRUCTED by
     combining cumulative binary T3 conversion probabilities; the direct 3-class T4 model was degenerate.
  2. The N(test) horizon breakdown PER SEED, flagged as the OLD baseline-T3 splits (not the
     horizon-stratified baseline_T4 splits) — the classes are tiny/imbalanced (seed 1 has zero ≥7y).

No re-fit: reads fusion_metrics.csv (seed-mean ± std) + fusion_predictions_all_T4.csv (pooled bACC and
the per-seed class breakdown), exactly like 02_summary.py.

Out: integration/T4/summary_t4_cleaned.{csv,png,pdf}
Run: PYTHONIOENCODING=utf-8 python integration/T4/04_summary_t4_cleaned.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score

from _t4_latex_lib import write_t4_latex

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = Path(__file__).resolve().parent
SEEDS = [0, 1, 2]
CLASSES = ["h_lt3", "h_3_7", "h_ge7"]

COHORTS = [("both_present", "both-present"), ("present_only", "all-available")]
ORDER = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg",
         "weighted_avg_J", "elastic_net"]
VAR_LABEL = {
    "clinical_only":  "clinical-only [BioClin-L ≤3y + MBERT-L ≤7y]",
    "mri_only":       "MRI-only [BrainMVP-ft (plus-orig)]",
    "equal_0.5":      "weighted-avg (equal)",
    "weighted_avg":   "weighted-avg (val bACC)",
    "weighted_avg_J": "weighted-avg (Youden J)",
    "elastic_net":    "elastic-net",
    "lr":             "logistic-reg",
}
# rendered metric columns (besides the two bACCs) — val bACC dropped
REST_METRICS = [("macro_auc", "macro-AUC"), ("macro_f1", "macro-F1"),
                ("f1_h_lt3", "F1 <3y"), ("f1_h_3_7", "F1 3–7y"), ("f1_h_ge7", "F1 ≥7y")]


def pooled_bacc():
    """pooled bACC per (variant, cohort) over all test patients across seeds (from the predictions)."""
    ap = pd.read_csv(HERE / "fusion_predictions_all_T4.csv")
    out = {}
    for (v, c), g in ap.groupby(["variant", "cohort"]):
        gg = g[g["pred"] >= 0]
        out[(v, c)] = balanced_accuracy_score(gg["y_true"], gg["pred"]) if len(gg) else np.nan
    return out


def n_breakdown():
    """Per-seed (<3y : 3-7y : ≥7y) TEST counts per cohort, from clinical_only (full coverage)."""
    ap = pd.read_csv(HERE / "fusion_predictions_all_T4.csv")
    bd = {}
    for ckey, _ in COHORTS:
        parts = []
        for s in SEEDS:
            g = ap[(ap.variant == "clinical_only") & (ap.cohort == ckey) & (ap.seed == s)]
            g = g[g.pred >= 0]
            vc = g["y_true"].value_counts().to_dict()
            parts.append(f"s{s} {vc.get(0, 0)}:{vc.get(1, 0)}:{vc.get(2, 0)}")
        bd[ckey] = ", ".join(parts)
    return bd


def aggregate(pooled):
    d = pd.read_csv(HERE / "fusion_metrics.csv")
    recs = []
    for ckey, clabel in COHORTS:
        cd = d[d["cohort"] == ckey]
        for v in ORDER:
            vd = cd[cd["variant"] == v]
            if vd.empty:
                continue
            r = vd.iloc[0]
            rec = {"ckey": ckey, "Cohort": clabel, "variant": v, "Method": VAR_LABEL[v],
                   "bacc": r["bacc_mean"], "bstd": r["bacc_std"],
                   "pooled": pooled.get((v, ckey), np.nan), "n": r["n_mean"]}
            for mk, _ in REST_METRICS:
                rec[mk] = r[f"{mk}_mean"]
            recs.append(rec)
    return pd.DataFrame(recs)


TITLE1 = "T4 Conversion Late Fusion (<3y / 3–7y / ≥7y)"
TITLE2 = "constructed from cumulative binary T3 probabilities, not a direct 3-class T4 model"


def render(df, subtitle, out_path, show_n=True):
    headers = ["Cohort", "Method", "TEST bACC", "pooled bACC"] + [lbl for _, lbl in REST_METRICS] \
        + (["n"] if show_n else [])
    LEFT_COLS = {0, 1}
    METRIC_KEYS = ["bacc", "pooled"] + [mk for mk, _ in REST_METRICS]   # columns 2..

    # best per cohort among FULL-cohort variants (exclude mri_only — different patients)
    best = {}
    for ckey, g in df.groupby("ckey", sort=False):
        gf = g[g.variant != "mri_only"]
        for key in ("bacc", "pooled", "macro_auc", "macro_f1"):
            col = gf[key]
            best[(ckey, key)] = col.idxmax() if col.notna().any() else -1

    body, rules, bolds, show_c = [], [], [], []
    prev_c = None
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
        for cj, key in enumerate(METRIC_KEYS, start=2):
            bkey = key if key in ("bacc", "pooled", "macro_auc", "macro_f1") else None
            bcols[cj] = (bkey is not None and best.get((r.ckey, bkey), -1) == idx)
        bolds.append(bcols)
        rules.append(prev_c is not None and r.ckey != prev_c)
        show_c.append(sc)
        prev_c = r.ckey

    COL_W = [1.35, 4.55, 1.62, 1.20, 1.15, 1.10, 0.95, 1.00, 0.95] + ([0.55] if show_n else [])
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


def main():
    pooled = pooled_bacc()
    df = aggregate(pooled)
    bd = n_breakdown()
    subtitle = (
        "bucket = [p(<3y), p(<7y)−p(<3y), 1−p(<7y)] per modality (monotone-clipped; direct 3-class T4 was degenerate)\n"
        "clinical@bl (≤3y BioClin-L, ≤7y MBERT-L) ⊕ MRI@m12 (BrainMVP-ft, plus-original)  ·  "
        "mean ± std across seeds 0/1/2\n"
        "both-present = complete-case (m12 MRI present)  ·  all-available = full clinical cohort "
        "(single modality proceeds where MRI absent)\n"
        "bACC: per-seed mean is high-variance (n ≈ 10–13/seed)  ·  pooled (all test patients) is a "
        "more trustworthy metric\n"
        "N(test) per seed on splits stratified by baseline diagnosis.  "
        f"both-present: {bd['both_present']}  ·  all-available: {bd['present_only']}   (<3y : 3–7y : ≥7y)")

    out_cols = ["Cohort", "Method", "bacc", "bstd", "pooled"] + [mk for mk, _ in REST_METRICS] + ["n"]
    df[out_cols].rename(columns={"bacc": "TEST_bACC_mean", "bstd": "TEST_bACC_std",
                                 "pooled": "pooled_bACC"}).to_csv(
        HERE / "summary_t4_cleaned.csv", index=False, encoding="utf-8")
    render(df, subtitle, HERE / "summary_t4_cleaned.png", show_n=True)
    render(df, subtitle, HERE / "summary_t4_cleaned_no_n.png", show_n=False)
    write_t4_latex(df, HERE / "summary_t4_cleaned.tex", TITLE1,
                   TITLE2 + " " + " ".join(subtitle.split("\n")),
                   "BioClin-L (≤3y) + MBERT-L (≤7y)", "BrainMVP-ft (plus-orig)",
                   "tab:t4_cumulative_fusion")
    print(f"  CSV: {HERE / 'summary_t4_cleaned.csv'}  ({len(df)} rows)")
    print("  N(test) breakdown — both-present:", bd["both_present"], "| all-available:", bd["present_only"])


if __name__ == "__main__":
    main()
