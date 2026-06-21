"""
07_render_survival_comparison_cleaned.py
========================================
Cleaned house-style copy of the T4 longitudinal-MAMBA survival comparison table (the raw table is
`D:/ADNI_BIDS_project/derivatives/mamba_survival/survival/comparison/survival_comparison_table.png`,
built by 02_survival_comparison.py). DejaVu-Serif / bordered / ruled / italic-header / dashed-divider /
bold-best house style.

Changes vs the raw table (per request):
  - keep ONLY the `valtest` rows (drop val / test);
  - split the cryptic `Model` name into two columns: `Model` (survival head) and `Aggregation`
    (sequence aggregator over the bl/m06/m12 visits);
  - drop the finetune / time-prepend / align ablations (A_q2/q3/q4) and A_q1_mamba2;
  - keep the four tabular baselines (RSF, Cox PH, Weibull AFT, Weibull piecewise) at the bottom (clinical
    baseline vector, no aggregator).

Reads the per-seed master CSV (no re-fit). Saves into the same comparison folder.

If `wald_significance_vs_rsf.csv` (from 08_wald_significance.py) is present, each non-RSF cell is
annotated with a Wald significance marker (* p<0.05, ** p<0.01, *** p<0.001) for the two-sided test of
that model vs RSF on that metric (paired, seed-stratified bootstrap SE). RSF is the reference row. The
matching `survival_comparison_table.tex` is regenerated with the same markers and an honest caption.

Out: D:/ADNI_BIDS_project/derivatives/mamba_survival/survival/comparison/survival_comparison_table_cleaned.{csv,png,pdf}
     integration/T4/mamba_long_clinical/survival_comparison_table.tex
Run: PYTHONIOENCODING=utf-8 python integration/T4/mamba_long_clinical/07_render_survival_comparison_cleaned.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
COMP = Path("D:/ADNI_BIDS_project/derivatives/mamba_survival/survival/comparison")
MASTER = COMP / "master_survival_comparison.csv"
WALD = COMP / "wald_significance_vs_rsf.csv"     # from 08_wald_significance.py (optional)
HERE = Path(__file__).resolve().parent          # also drop copies next to the script
OUT_DIRS = [COMP, HERE]
TEX_OUT = HERE / "survival_comparison_table.tex"


def load_wald():
    """(model_key, metric) -> significance stars, from 08's bootstrap-Wald CSV. {} if absent."""
    if not WALD.exists():
        return {}
    w = pd.read_csv(WALD)
    return {(str(r.model), str(r.metric)): ("" if pd.isna(r.sig) else str(r.sig))
            for r in w.itertuples()}

# model key -> (group, Model = survival head, Aggregation = sequence aggregator)
#   group: "mamba_agg" (aggregation ablation), "mamba_head" (head ablation), "baseline"
ROWS = [
    ("A_ctrl_gru",              "mamba_agg",  "Weibull piecewise", "Gated Recurrent Unit"),
    ("A_ctrl_meanpool",         "mamba_agg",  "Weibull piecewise", "Meanpool"),
    ("A_default_mamba1_frozen", "mamba_agg",  "Weibull piecewise", "MAMBA (frozen)"),
    ("B_cox",                   "mamba_head", "Cox PH",            "MAMBA (frozen)"),
    ("B_soden",                 "mamba_head", "SODEN",             "MAMBA (frozen)"),
    ("B_weibull_aft",           "mamba_head", "Weibull AFT",       "MAMBA (frozen)"),
    ("RSF",                     "baseline",   "RSF",               "—"),
    ("Cox PH",                  "baseline",   "Cox PH",            "—"),
    ("Weibull AFT",             "baseline",   "Weibull AFT",       "—"),
    ("Weibull piecewise",       "baseline",   "Weibull piecewise", "—"),
]
METRICS = [("c_index", "C-index", "max"), ("c_index_ipcw", "IPCW C", "max"), ("ibs", "IBS", "min"),
           ("auc_3y", "AUC@3y", "max"), ("auc_5y", "AUC@5y", "max"),
           ("auc_7y", "AUC@7y", "max"), ("auc_10y", "AUC@10y", "max")]


def aggregate():
    d = pd.read_csv(MASTER)
    d = d[d["split"] == "valtest"]
    recs = []
    for key, group, model, agg in ROWS:
        g = d[d["model"] == key]
        if g.empty:
            raise SystemExit(f"missing valtest rows for model={key}")
        rec = {"key": key, "group": group, "Model": model, "Aggregation": agg,
               "n": g["n"].mean(), "events": g["events"].mean()}
        for m, _, _ in METRICS:
            rec[f"{m}_m"] = g[m].mean(); rec[f"{m}_s"] = g[m].std(ddof=1)
        recs.append(rec)
    return pd.DataFrame(recs)


TITLE = "T4 Survival — longitudinal-MAMBA heads vs tabular baselines (CN+MCI → AD)"
SUBTITLE = (
    "held-out = val+test pooled  ·  mean ± std over 3 seeds  ·  bold = best per column\n"
    "Model = survival head;  Aggregation = sequence aggregator over the bl/m06/m12 visits (risk = 1 − S(t_ref))\n"
    "MAMBA = deep heads on clinical embeddings  ·  baselines = RSF / Cox PH / Weibull AFT / Weibull piecewise "
    "on an 8-feature clinical baseline vector  ·  N(val+test) = 108 (≈25 events)\n"
    "* p<0.05  ** p<0.01  *** p<0.001 : two-sided Wald test vs RSF (reference), paired seed-stratified bootstrap SE")


def render(df, out_path, wald=None):
    wald = wald or {}
    headers = ["Model", "Aggregation", "N (events)"] + [lbl for _, lbl, _ in METRICS]
    LEFT_COLS = {0, 1}
    # best per metric column (max, or min for IBS)
    best = {}
    for m, _, mode in METRICS:
        col = df[f"{m}_m"].to_numpy(float)
        best[m] = int(np.nanargmin(col) if mode == "min" else np.nanargmax(col))

    body, rules, bolds = [], [], []
    prev_g = None
    for i, r in df.iterrows():
        # rule above: strong dashed before baselines; lighter dashed between agg- and head-ablation
        if prev_g is None:
            rule = None
        elif r.group == "baseline" and prev_g != "baseline":
            rule = "strong"
        elif r.group != prev_g:
            rule = "light"
        else:
            rule = None
        cells = [r.Model, r.Aggregation, f"{r['n']:.0f} ({r['events']:.0f})"]
        for m, _, _ in METRICS:
            mark = wald.get((r.key, m), "")            # "" for RSF (reference) and absent files
            cells.append(f"{r[f'{m}_m']:.3f} ± {r[f'{m}_s']:.3f}{mark}")
        body.append(cells)
        bolds.append({3 + j: (best[m] == i) for j, (m, _, _) in enumerate(METRICS)})
        rules.append(rule); prev_g = r.group

    COL_W = [1.95, 2.35, 1.20] + [1.62] * len(METRICS)
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_LINES = SUBTITLE.count("\n") + 1
    SUB_H = 0.28 * SUB_LINES
    TITLE_H, HEAD_H, ROW_H = 0.55 + SUB_H, 0.42, 0.46
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
                    fontsize=9.0, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.0, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i] == "strong":
            hline(y, lw=0.9, ls=(0, (4, 3)))
        elif rules[i] == "light":
            hline(y, lw=0.5, ls=(0, (1, 3)))
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center", fontsize=8.8)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=8.8,
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


# ── LaTeX export (booktabs; Longitudinal / Cross-sectional groups; Wald markers vs RSF) ──
TEX_HEAD = [("c_index", r"C-index $\uparrow$"), ("c_index_ipcw", r"C-idx IPCW $\uparrow$"),
            ("ibs", r"IBS $\downarrow$"), ("auc_3y", r"AUC@3y $\uparrow$"), ("auc_5y", r"AUC@5y $\uparrow$"),
            ("auc_7y", r"AUC@7y $\uparrow$"), ("auc_10y", r"AUC@10y $\uparrow$")]
_TEX_STAR = {"*": "^{*}", "**": "^{**}", "***": "^{***}"}


def write_tex(df, wald, out_path):
    best = {}
    for m, _, mode in METRICS:
        col = df[f"{m}_m"].to_numpy(float)
        best[m] = float(np.nanmin(col) if mode == "min" else np.nanmax(col))

    def cell(r, m):
        v, s = getattr(r, f"{m}_m"), getattr(r, f"{m}_s")
        core = f"{v:.2f} \\pm {s:.2f}"
        if round(v, 2) == round(best[m], 2):
            core = r"\mathbf{" + core + "}"
        return f"${core}{_TEX_STAR.get(wald.get((r.key, m), ''), '')}$"

    def row(r):
        agg = "---" if str(r.Aggregation) in ("—", "-") else str(r.Aggregation)
        return f"{r.Model} & {agg} & " + " & ".join(cell(r, m) for m, _, _ in METRICS) + r" \\"

    longi = df[df.group != "baseline"]; cross = df[df.group == "baseline"]
    lines = [
        r"% Survival model comparison: longitudinal MAMBA aggregation vs cross-sectional baselines",
        r"% Metrics: C-index, C-index IPCW, IBS, AUC@3y/5y/7y/10y (mean $\pm$ std across seeds 0/1/2)",
        r"% Significance: two-sided Wald test vs RSF; SE from paired seed-stratified bootstrap (B=2000). Regenerated by 07_render_survival_comparison_cleaned.py.",
        r"\begin{table}[ht]", r"\centering", r"\resizebox{\textwidth}{!}{%", r"\normalsize",
        r"\begin{tabular}{llccccccc}", r"\toprule",
        r"& & \multicolumn{7}{c}{\textbf{Validation + Test Performance}} \\",
        r"\textbf{Model} & \textbf{Aggregation} & " + " & ".join(rf"\textbf{{{h}}}" for _, h in TEX_HEAD) + r" \\",
        r"\midrule", r"\multicolumn{9}{l}{\textit{Longitudinal}} \\",
        *[row(r) for r in longi.itertuples()],
        r"\midrule", r"\multicolumn{9}{l}{\textit{Cross-sectional}} \\",
        *[row(r) for r in cross.itertuples()],
        r"\bottomrule", r"\end{tabular}%", r"}",
        r"\caption{Survival model performance on the combined validation and test set "
        r"(N=108, mean events$\approx$24.7 per seed). Longitudinal models aggregate clinical visit "
        r"embeddings over baseline, 6-month and 12-month visits; cross-sectional models use clinical "
        r"features from the baseline visit only. Results reported as mean $\pm$ std across seeds 0/1/2; "
        r"bold = best per column. Significance markers ($^{*}p{<}0.05$, $^{**}p{<}0.01$, "
        r"$^{***}p{<}0.001$) denote a two-sided Wald test of each model against the RSF baseline on that "
        r"metric, with the standard error of the difference estimated by a paired, seed-stratified "
        r"bootstrap ($B{=}2000$) of the per-patient predictions; Benjamini--Hochberg FDR-adjusted "
        r"$q$-values are reported in the supplementary \texttt{wald\_significance\_vs\_rsf.csv}. "
        r"IPCW: Inverse Probability of Censoring Weighting. GRU: Gated Recurrent Unit. "
        r"SODEN: Survival Neural ODE. RSF: Random Survival Forest. AFT: Accelerated Failure Time.}",
        r"\label{tab:survival_mamba_val}", r"\end{table}", ""]
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  TeX: {out_path}")


def main():
    df = aggregate()
    wald = load_wald()
    if not wald:
        print("  [note] wald_significance_vs_rsf.csv not found — rendering without significance markers.")
    out_cols = ["Model", "Aggregation", "n", "events"] \
        + [f"{m}_m" for m, _, _ in METRICS] + [f"{m}_s" for m, _, _ in METRICS]
    for d in OUT_DIRS:
        df[out_cols].to_csv(d / "survival_comparison_table_cleaned.csv", index=False, encoding="utf-8")
        render(df, d / "survival_comparison_table_cleaned.png", wald)
        print(f"  CSV: {d / 'survival_comparison_table_cleaned.csv'}")
    write_tex(df, wald, TEX_OUT)
    print(f"  ({len(df)} rows)")
    print(df[["Model", "Aggregation", "c_index_m", "ibs_m"]].to_string(index=False))


if __name__ == "__main__":
    main()
