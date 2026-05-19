"""
19b_render_model_comparison_table.py
=====================================
Aggregate the ``frozen_model_comparison/`` metrics.json tree written by
``19_compare_frozen_models.py`` and render a B&W publication table.

  rows  = (model, aggregation, head)
  cols  = test metrics (AUC, BalAcc, F1, Sens, Spec, Precision, Recall)
          grouped by label-mode cohort, mean ± std over seeds.

Adapts the collect/summarise logic of ``16_aggregate_classification_results.py``
(adding ``model`` and ``head`` group axes) and the B&W house-style renderer of
``16b_render_publication_table.py`` (DejaVu Serif, white fill, thin outline,
no bold).

Outputs (under --root):
  comparison_long.csv, comparison_summary.csv,
  model_comparison_table.png, model_comparison_table.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "DejaVu Serif"

METRICS = [
    ("test_auc",               "AUC"),
    ("test_balanced_accuracy", "BalAcc"),
    ("test_f1",                "F1"),
    ("test_sensitivity",       "Sensitivity"),
    ("test_specificity",       "Specificity"),
    ("test_precision",         "Precision"),
    ("test_recall",            "Recall"),
]
RAW_METRICS = ["auc", "balanced_accuracy", "f1", "precision", "recall",
               "sensitivity", "specificity"]

MODEL_ORDER = ["bmfm_base_ref", "ntv2", "caduceus_ps_8k", "caduceus_ph_8k",
               "caduceus_ps_131k", "caduceus_ph_131k"]
MODEL_PRETTY = {
    "bmfm_base_ref":    "BMFM-DNA-REF (frozen)",
    "ntv2":             "Nucleotide Transformer v2",
    "caduceus_ps_8k":   "Caduceus-PS (8 kb)",
    "caduceus_ph_8k":   "Caduceus-PH (8 kb)",
    "caduceus_ps_131k": "Caduceus-PS (131 kb)",
    "caduceus_ph_131k": "Caduceus-PH (131 kb)",
}
AGG_ORDER = ["chrom_mean_then_mean", "chrom_mean_then_attn"]
AGG_PRETTY = {
    "chrom_mean_then_mean": "chrom-mean → mean",
    "chrom_mean_then_attn": "chrom-mean → attn",
}
HEAD_ORDER = ["mlp2", "mlp3", "rf", "xgb"]
HEAD_PRETTY = {"mlp2": "2-layer MLP", "mlp3": "3-layer MLP",
               "rf": "RandomForest", "xgb": "XGBoost"}
COHORT_PRETTY = {
    "ever_convert": ("Longitudinal conversion  (stable-CN vs ever-MCI/AD)",
                     "Ntest≈62, stable-CN≈15, ever-conv≈47"),
    "bl_multi":     ("Baseline diagnosis  (CN vs MCI+AD)",
                     "Ntest≈62, CN≈21, MCI≈37, AD≈4"),
}


def collect(root: Path) -> pd.DataFrame:
    rows = []
    for mfile in root.rglob("metrics.json"):
        with open(mfile) as f:
            m = json.load(f)
        row = {
            "label_mode": m.get("label_mode", "ever_convert"),
            "model": m.get("model", m.get("upstream_tag")),
            "aggregation": m["aggregation"],
            "head": m.get("head", m.get("clf", "?")),
            "seed": m["seed"],
        }
        for k in RAW_METRICS:
            row[f"val_{k}"] = m["val"].get(k)
            row[f"test_{k}"] = m["test"].get(k)
        row["best_val_auc"] = m.get("fit_info", {}).get("best_val_auc")
        rows.append(row)
    if not rows:
        raise SystemExit(f"No metrics.json found under {root}")
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (lm, model, agg, head), grp in df.groupby(
            ["label_mode", "model", "aggregation", "head"]):
        row = {"label_mode": lm, "model": model, "aggregation": agg,
               "head": head, "n_seeds": len(grp)}
        for k in RAW_METRICS:
            vals = grp[f"test_{k}"].dropna().to_numpy()
            if len(vals):
                row[f"test_{k}"] = f"{np.mean(vals):.3f} ± {np.std(vals, ddof=0):.3f}"
                row[f"test_{k}_mean"] = float(np.mean(vals))
                row[f"test_{k}_std"] = float(np.std(vals, ddof=0))
            else:
                row[f"test_{k}"] = ""
        out.append(row)
    return pd.DataFrame(out)


def _ord(series, order):
    return pd.Categorical(series, categories=order, ordered=True)


def pivot(df: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    df = df.copy()
    df["model"] = _ord(df["model"], MODEL_ORDER)
    df["aggregation"] = _ord(df["aggregation"], AGG_ORDER)
    df["head"] = _ord(df["head"], HEAD_ORDER)
    rows = []
    for (model, agg, head), grp in df.groupby(
            ["model", "aggregation", "head"], observed=True, sort=True):
        r = {"model": model, "aggregation": agg, "head": head}
        for lm in cohorts:
            sub = grp[grp["label_mode"] == lm]
            for mcol, mname in METRICS:
                if len(sub):
                    r[(lm, mname)] = sub.iloc[0].get(mcol, "")
                else:
                    r[(lm, mname)] = ""
        rows.append(r)
    return pd.DataFrame(rows)


def render_png(out: pd.DataFrame, cohorts: list[str], out_png: Path) -> None:
    n_rows, n_metrics, n_coh = len(out), len(METRICS), len(cohorts)
    cell_h_in, header_in, cohort_in, title_in, margin_in = 0.42, 0.55, 0.50, 0.85, 0.30
    fig_w = 6.5 + 1.55 * n_metrics * n_coh
    fig_h = title_in + cohort_in + header_in + n_rows * cell_h_in + 2 * margin_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    left_widths = [0.16, 0.13, 0.11]                       # Model, Aggregation, Head
    left_end = sum(left_widths)
    metric_w = (1.0 - left_end) / (n_metrics * n_coh)
    col_widths = left_widths + [metric_w] * (n_metrics * n_coh)
    col_x = [0.0]
    for w in col_widths:
        col_x.append(col_x[-1] + w)

    fh = fig_h
    title_h, cohort_h, header_h, margin_h = (title_in/fh, cohort_in/fh,
                                             header_in/fh, margin_in/fh)
    header_top = 1.0 - title_h - cohort_h
    body_top = header_top - header_h
    body_h = body_top - margin_h
    cell_h = body_h / max(1, n_rows)

    def box(x, y, w, h, lw=0.4):
        ax.add_patch(plt.Rectangle((x, y), w, h, transform=ax.transAxes,
                                   facecolor="white", edgecolor="black", linewidth=lw))

    def txt(x, y, s, fs=8.5, ha="center"):
        ax.text(x, y, s, ha=ha, va="center", fontsize=fs, transform=ax.transAxes)

    txt(0.5, 1.0 - title_h*0.35,
        "Frozen foundation-model embeddings: classifier-head comparison", fs=14)
    txt(0.5, 1.0 - title_h*0.75,
        "per-window pool = mean;  mean ± std over seeds 0/1/2;  "
        "80/10/10 patient split;  chrom-mean→attn is MLP-only (trainable)",
        fs=9.5)

    for ci, lm in enumerate(cohorts):
        x0 = left_end + ci * n_metrics * metric_w
        w = n_metrics * metric_w
        box(x0, header_top, w, cohort_h, lw=0.7)
        pretty, stats = COHORT_PRETTY.get(lm, (lm, ""))
        txt(x0 + w/2, header_top + cohort_h*0.68, pretty, fs=11)
        txt(x0 + w/2, header_top + cohort_h*0.26, stats, fs=8.5)
    box(0, header_top, left_end, cohort_h, lw=0.7)

    headers = ["Model", "Aggregation", "Head"] + \
              [m[1] for _ in cohorts for m in METRICS]
    for i, h in enumerate(headers):
        box(col_x[i], body_top, col_widths[i], header_h, lw=0.7)
        txt(col_x[i] + col_widths[i]/2, body_top + header_h/2, h, fs=9.5)

    for ri in range(n_rows):
        y0 = body_top - (ri + 1) * cell_h
        row = out.iloc[ri]
        left_vals = [MODEL_PRETTY.get(str(row["model"]), str(row["model"])),
                     AGG_PRETTY.get(str(row["aggregation"]), str(row["aggregation"])),
                     HEAD_PRETTY.get(str(row["head"]), str(row["head"]))]
        for ci, val in enumerate(left_vals):
            box(col_x[ci], y0, col_widths[ci], cell_h)
            txt(col_x[ci] + (0.005 if ci == 0 else col_widths[ci]/2),
                y0 + cell_h/2, val, fs=8.5, ha=("left" if ci == 0 else "center"))
        ci = 3
        for lm in cohorts:
            for _, mname in METRICS:
                box(col_x[ci], y0, col_widths[ci], cell_h)
                txt(col_x[ci] + col_widths[ci]/2, y0 + cell_h/2,
                    str(row[(lm, mname)]), fs=8.0)
                ci += 1
        if ri < n_rows - 1 and str(out.iloc[ri]["model"]) != str(out.iloc[ri+1]["model"]):
            ax.plot([0, 1], [y0, y0], transform=ax.transAxes,
                    linewidth=1.0, color="black")

    fig.savefig(out_png, dpi=150, facecolor="white",
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {out_png}  (figsize={fig_w:.1f} x {fig_h:.1f} in, {n_rows} rows)")


def render_md(out: pd.DataFrame, cohorts: list[str], out_md: Path) -> None:
    lines = ["# Frozen foundation-model embeddings — classifier-head comparison",
             "",
             "Per-window pool = mean. Mean ± std over seeds 0/1/2. "
             "`chrom-mean → attn` is MLP-only (trainable aggregation).", ""]
    for lm in cohorts:
        pretty, stats = COHORT_PRETTY.get(lm, (lm, ""))
        lines.append(f"- **{pretty}** — {stats}")
    lines.append("")
    head = ["Model", "Aggregation", "Head"]
    for lm in cohorts:
        for _, mname in METRICS:
            head.append(f"{lm}:{mname}")
    lines.append("| " + " | ".join(head) + " |")
    lines.append("|" + "|".join(["---"] * len(head)) + "|")
    for ri in range(len(out)):
        row = out.iloc[ri]
        cells = [MODEL_PRETTY.get(str(row["model"]), str(row["model"])),
                 AGG_PRETTY.get(str(row["aggregation"]), str(row["aggregation"])),
                 HEAD_PRETTY.get(str(row["head"]), str(row["head"]))]
        for lm in cohorts:
            for _, mname in METRICS:
                cells.append(str(row[(lm, mname)]))
        lines.append("| " + " | ".join(cells) + " |")
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out_md}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("frozen_model_comparison"))
    args = ap.parse_args()

    long_df = collect(args.root)
    long_df.to_csv(args.root / "comparison_long.csv", index=False)
    print(f"Wrote {args.root/'comparison_long.csv'}  ({len(long_df)} rows)")

    sum_df = summarise(long_df)
    sum_df.to_csv(args.root / "comparison_summary.csv", index=False)
    print(f"Wrote {args.root/'comparison_summary.csv'}  ({len(sum_df)} rows)")

    cohorts = [lm for lm in ("ever_convert", "bl_multi")
               if lm in set(sum_df["label_mode"])]
    piv = pivot(sum_df, cohorts)
    render_png(piv, cohorts, args.root / "model_comparison_table.png")
    render_md(piv, cohorts, args.root / "model_comparison_table.md")

    cols = ["label_mode", "model", "aggregation", "head", "n_seeds",
            "test_auc", "test_balanced_accuracy", "test_f1"]
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print("\n" + "=" * 90)
        print(sum_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
