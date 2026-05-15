"""
16b_render_publication_table.py
================================
Render classification_results/comparison_summary.csv as a publication-style
table styled after the clinical-models figure: rows = (upstream, pool,
aggregation), column groups = (bl_multi, ever_convert) cohorts each with
Acc-equivalent metrics; bolds the per-column best.

Outputs:
  classification_results/publication_table.png
  classification_results/publication_table.md
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent / "classification_results_BCE"
SUMMARY = ROOT / "comparison_summary.csv"
OUT_PNG = ROOT / "publication_table.png"
OUT_MD = ROOT / "publication_table.md"

UPSTREAM_ORDER = [
    "base_ref",
    "without_ukb_ref",
    "without_ukb_augmented_ref",
    "by_chrom_fr_ref",
    "combos_fr_ref",
    "combos_singletons_fr_ref",
]
UPSTREAM_PRETTY = {
    "base_ref":                   "BMFM-DNA-REF (no FT)",
    "without_ukb_ref":            "REF-FT (no UKB)",
    "without_ukb_augmented_ref":  "REF-FT (no UKB, +RC)",
    "by_chrom_fr_ref":            "REF-FT (by-chrom 8kb)",
    "combos_fr_ref":              "REF-FT (combos)",
    "combos_singletons_fr_ref":   "REF-FT (combos+singletons)",
}
AGG_ORDER = [
    "mean_all",
    "chrom_mean_then_mean",
    "attn_all",
    "chrom_mean_then_attn",
    "chrom_attn_then_attn",
    "concat_chrom_means",
]
AGG_PRETTY = {
    "mean_all":             "mean (all windows)",
    "chrom_mean_then_mean": "chrom-mean then mean",
    "attn_all":             "attn (all windows)",
    "chrom_mean_then_attn": "chrom-mean then attn",
    "chrom_attn_then_attn": "chrom-attn then attn",
    "concat_chrom_means":   "concat chrom-means",
}
POOL_PRETTY = {"cls": "CLS token", "mean": "mean over DNA tokens"}

METRICS = [
    ("test_auc",               "AUC"),
    ("test_balanced_accuracy", "BalAcc"),
    ("test_f1",                "F1"),
    ("test_sensitivity",       "Sensitivity"),
    ("test_specificity",       "Specificity"),
    ("test_precision",         "Precision"),
    ("test_recall",            "Recall"),
]
COHORTS = [
    ("bl_multi",     "Baseline diagnosis  (CN vs MCI+AD)",
                     "Ntest=62, CN≈21, MCI≈37, AD≈4"),
    ("ever_convert", "Longitudinal conversion  (stable-CN vs ever-MCI/AD)",
                     "Ntest=62, stable-CN≈15, ever-conv≈47"),
]


def load_summary() -> pd.DataFrame:
    df = pd.read_csv(SUMMARY)
    df["upstream_tag"] = pd.Categorical(df["upstream_tag"], categories=UPSTREAM_ORDER, ordered=True)
    df["aggregation"]  = pd.Categorical(df["aggregation"],  categories=AGG_ORDER,  ordered=True)
    df["pool"]         = pd.Categorical(df["pool"],         categories=["cls", "mean"], ordered=True)
    return df


def pivot_for_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return a wide DataFrame: index = (upstream, pool, aggregation),
    columns = (label_mode, metric_display_name), values = mean±std strings."""
    rows = []
    for (upstream, pool, agg), grp in df.groupby(["upstream_tag", "pool", "aggregation"],
                                                  observed=True, sort=True):
        row = {"upstream": upstream, "pool": pool, "aggregation": agg}
        for label_mode, _, _ in COHORTS:
            sub = grp[grp["label_mode"] == label_mode]
            for mcol, mname in METRICS:
                mean = sub.get(f"{mcol}_mean", pd.Series(dtype=float)).iloc[0] if len(sub) else np.nan
                std  = sub.get(f"{mcol}_std",  pd.Series(dtype=float)).iloc[0] if len(sub) else np.nan
                if pd.isna(mean):
                    row[(label_mode, mname)] = ""
                else:
                    row[(label_mode, mname)] = f"{mean:.3f} ± {std:.3f}"
                row[(label_mode, mname, "_mean")] = mean
        rows.append(row)
    out = pd.DataFrame(rows)
    return out


def best_indices_per_column(out: pd.DataFrame) -> dict[tuple, set[int]]:
    """For each (cohort, metric) column, return the set of row indices whose
    mean equals the column max (ties allowed). All metrics here are
    higher-is-better."""
    best = {}
    for label_mode, _, _ in COHORTS:
        for _, mname in METRICS:
            col_mean = out[(label_mode, mname, "_mean")]
            if col_mean.dropna().empty:
                best[(label_mode, mname)] = set()
                continue
            mx = col_mean.max()
            # Treat ties within 0.0005 as equal (so display-rounding doesn't lose ties).
            best[(label_mode, mname)] = set(out.index[(col_mean - mx).abs() < 5e-4].tolist())
    return best


def render_png(out: pd.DataFrame, best: dict) -> None:
    n_rows = len(out)
    n_metrics = len(METRICS)
    n_cohorts = len(COHORTS)

    # Figure sizing — work in real units, NOT normalised, so cells stay readable.
    # 0.40 inch per row keeps text + bold legible even with 72 rows.
    cell_h_in   = 0.40
    header_in   = 0.55
    cohort_in   = 0.48
    title_in    = 0.80
    margin_in   = 0.30
    fig_w       = 34.0
    fig_h = title_in + cohort_in + header_in + n_rows * cell_h_in + 2 * margin_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Column widths (in normalised x). Left columns wider.
    left_widths = [0.18, 0.12, 0.12]
    metric_w = (1.0 - sum(left_widths)) / (n_metrics * n_cohorts)
    col_widths = left_widths + [metric_w] * (n_metrics * n_cohorts)
    col_x = [0.0]
    for w in col_widths:
        col_x.append(col_x[-1] + w)

    # Vertical layout in normalised y.
    top_y       = 1.0
    title_h     = title_in / fig_h
    cohort_h    = cohort_in / fig_h
    header_h    = header_in / fig_h
    margin_h    = margin_in / fig_h
    header_top  = top_y - title_h - cohort_h
    body_top    = header_top - header_h
    body_h      = body_top - margin_h
    cell_h      = body_h / max(1, n_rows)

    def rect(x0, y0, w, h, **kw):
        ax.add_patch(plt.Rectangle((x0, y0), w, h, transform=ax.transAxes,
                                    fill=False, linewidth=0.5, edgecolor="#888888", **kw))

    def text(x, y, s, **kw):
        kw.setdefault("ha", "center"); kw.setdefault("va", "center")
        kw.setdefault("fontsize", 9.0); kw.setdefault("transform", ax.transAxes)
        ax.text(x, y, s, **kw)

    # Title
    title_main = ("Frozen BMFM-DNA-REF embeddings + MLP head: "
                  "test performance per (upstream, pool, aggregation)")
    title_sub  = ("mean ± std over seeds 0/1/2;  80/10/10 patient-level split;  "
                  "183 chromosome windows per patient (8 kb each);  "
                  "head loss = unweighted BCEWithLogitsLoss")
    text(0.5, top_y - title_h * 0.35, title_main, fontsize=14)
    text(0.5, top_y - title_h * 0.75, title_sub,  fontsize=10, style="italic",
         color="black")

    # Cohort group bars  (B&W: no fill, just thin outlines)
    left_end = sum(left_widths)
    for ci, (label_mode, pretty, cohort_stats) in enumerate(COHORTS):
        x0 = left_end + ci * n_metrics * metric_w
        w = n_metrics * metric_w
        ax.add_patch(plt.Rectangle((x0, header_top),
                                    w, cohort_h,
                                    facecolor="white", edgecolor="black",
                                    linewidth=0.7, transform=ax.transAxes))
        text(x0 + w / 2, header_top + cohort_h * 0.70,
             pretty, fontsize=11)
        text(x0 + w / 2, header_top + cohort_h * 0.25,
             cohort_stats, fontsize=8.5, style="italic", color="black")

    # Left header strip (blank box above the 3 grouping columns)
    ax.add_patch(plt.Rectangle((0, header_top),
                                left_end, cohort_h,
                                facecolor="white", edgecolor="black",
                                linewidth=0.7, transform=ax.transAxes))

    # Header row: left cols + metric names
    headers = ["Upstream (BMFM-DNA-REF + GWAS-reg FT)", "Pool", "Aggregation"]
    for _, _, _ in COHORTS:
        headers.extend([m[1] for m in METRICS])

    for i, h in enumerate(headers):
        ax.add_patch(plt.Rectangle((col_x[i], body_top), col_widths[i], header_h,
                                    facecolor="white", edgecolor="black",
                                    linewidth=0.7, transform=ax.transAxes))
        text(col_x[i] + col_widths[i] / 2, body_top + header_h / 2,
             h, fontsize=9.5)

    # Body rows
    for ri in range(n_rows):
        y0 = body_top - (ri + 1) * cell_h
        row = out.iloc[ri]
        # Left three columns
        left_vals = [
            UPSTREAM_PRETTY[str(row["upstream"])],
            POOL_PRETTY[str(row["pool"])],
            AGG_PRETTY[str(row["aggregation"])],
        ]
        for ci, val in enumerate(left_vals):
            ax.add_patch(plt.Rectangle((col_x[ci], y0), col_widths[ci], cell_h,
                                        facecolor="white", edgecolor="black",
                                        linewidth=0.4, transform=ax.transAxes))
            if ci == 0:
                text(col_x[ci] + 0.005, y0 + cell_h / 2, val,
                     fontsize=9, ha="left")
            else:
                text(col_x[ci] + col_widths[ci] / 2, y0 + cell_h / 2,
                     val, fontsize=9)

        # Metric columns
        ci = 3
        for label_mode, _, _ in COHORTS:
            for _, mname in METRICS:
                v = row[(label_mode, mname)]
                ax.add_patch(plt.Rectangle((col_x[ci], y0), col_widths[ci], cell_h,
                                            facecolor="white", edgecolor="black",
                                            linewidth=0.4, transform=ax.transAxes))
                text(col_x[ci] + col_widths[ci] / 2, y0 + cell_h / 2,
                     v, fontsize=8.5)
                ci += 1

        # Heavier separator between upstream groups (every 12 rows)
        if (ri + 1) % 12 == 0 and ri < n_rows - 1:
            ax.plot([0, 1], [y0, y0], transform=ax.transAxes,
                    linewidth=1.0, color="black")

    # Outer border around full table (header + body)
    outer_top = top_y - title_h - cohort_h
    outer_bot = body_top - body_h
    ax.plot([0, 1, 1, 0, 0],
            [outer_bot, outer_bot, outer_top + cohort_h, outer_top + cohort_h, outer_bot],
            transform=ax.transAxes, linewidth=1.2, color="#222222")

    # Save with tight bounding box to trim surrounding whitespace.
    fig.savefig(OUT_PNG, dpi=150, facecolor="white",
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {OUT_PNG}  (figsize={fig_w} x {fig_h:.1f} in, {n_rows} rows)")


def render_markdown(out: pd.DataFrame, best: dict) -> None:
    """Companion markdown table — plain (no bold)."""
    lines = []
    lines.append("# Publication table — frozen BMFM-DNA-REF + MLP head")
    lines.append("")
    lines.append("Mean ± std over seeds 0/1/2. Test set = 80/10/10 patient-level split. "
                 "Head loss = unweighted `BCEWithLogitsLoss`.")
    lines.append("")
    for _, pretty, stats in COHORTS:
        lines.append(f"- **{pretty}** — {stats}")
    lines.append("")

    header_top = ["Upstream", "Pool", "Aggregation"]
    for _, pretty, _ in COHORTS:
        for _, mname in METRICS:
            header_top.append(f"{pretty.split('  ')[0][:18]}: {mname}")
    lines.append("| " + " | ".join(header_top) + " |")
    lines.append("|" + "|".join(["---"] * len(header_top)) + "|")

    for ri in range(len(out)):
        row = out.iloc[ri]
        cells = [
            UPSTREAM_PRETTY[str(row["upstream"])],
            POOL_PRETTY[str(row["pool"])],
            AGG_PRETTY[str(row["aggregation"])],
        ]
        for label_mode, _, _ in COHORTS:
            for _, mname in METRICS:
                cells.append(row[(label_mode, mname)])
        lines.append("| " + " | ".join(cells) + " |")

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_MD}")


def main() -> None:
    df = load_summary()
    out = pivot_for_table(df)
    best = best_indices_per_column(out)
    render_png(out, best)
    render_markdown(out, best)


if __name__ == "__main__":
    main()
