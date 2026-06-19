r"""
30b_aggregate_corrected_gtex_sweep.py
======================================
Aggregate the 81-run corrected-GTEx func-integration sweep
(9 modes x 3 families x 3 seeds) into a per-config mean+/-std table.

Inputs:
  outputs/corrected_gtex/diff_attn_v3/{head}/{run_name}/metrics.json

Outputs (under outputs/corrected_gtex/sweep_table/):
  - sweep_table.csv           # per-config mean +/- std + per-seed values
  - sweep_table_summary.csv   # one row per (family, mode), mean +/- std only
  - sweep_table.png           # rendered table sorted by val BalAcc desc

Family canonical names:
  - bmfm_mlp  = BMFM-SNP 2L MLP @ 101bp (chrom_hier, w=256, d=0.1, lr=0.001)
  - ntv2_mlp  = NTv2 2L MLP @ 101bp (chrom_hier, w=512, d=0.3, lr=0.0003)
  - bmfm_xgb  = BMFM-SNP XGB @ 1001bp (global_attn, xgb_lr=0.05, depth=4, n=20)

Modes (9 total, ordered by mechanism class):
  none, per_snp_zscored_25, per_snp_summaries_12,
  attn_bias_single, attn_bias_per_modality, attn_bias_single_plus_func_imp_abs,
  attn_bias_per_modality_plus_per_snp_summaries_12,
  value_mult, attn_bias_per_modality_plus_value_mult
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

# Defaults (local Windows). Override via env vars on Colab/HPC.
REPO = Path(__file__).parent.parent.parent
SWEEP_ROOT = Path(os.environ.get("SWEEP_ROOT",
    str(REPO / "Transformers_XAI/snp_pipeline/outputs/corrected_gtex/diff_attn_v3")))
OUT_DIR = Path(os.environ.get("OUT_DIR",
    str(REPO / "Transformers_XAI/snp_pipeline/outputs/corrected_gtex/sweep_table")))

# Family identification rule: (seq_length, model, head) -> family label
FAMILY_RULES = [
    (("101bp",  "bmfm_snp", "mlp2"), "BMFM-SNP 2L MLP (101bp)"),
    (("101bp",  "ntv2",     "mlp2"), "NTv2 2L MLP (101bp)"),
    (("1001bp", "bmfm_snp", "xgb"),  "BMFM-SNP XGB (1001bp)"),
]

MODE_ORDER = [
    "none",
    "per_snp_zscored_25",
    "per_snp_summaries_12",
    "attn_bias_single",
    "attn_bias_per_modality",
    "attn_bias_single_plus_func_imp_abs",
    "attn_bias_per_modality_plus_per_snp_summaries_12",
    "value_mult",
    "attn_bias_per_modality_plus_value_mult",
]


def _family_label(seq_length: str, model: str, head: str) -> str | None:
    for (sl, m, h), lab in FAMILY_RULES:
        if (sl, m, h) == (seq_length, model, head):
            return lab
    return None


def _load_one(metrics_path: Path) -> dict | None:
    try:
        d = json.loads(metrics_path.read_text())
    except Exception as e:
        print(f"  [warn] failed {metrics_path}: {e}")
        return None
    cfg = d.get("config", {})
    seq_length = cfg.get("seq_length", "")
    model = cfg.get("model", "")
    head = cfg.get("head", "")
    mode = cfg.get("func_integration_mode", "")
    seed = int(cfg.get("seed", -1))
    family = _family_label(seq_length, model, head)
    if family is None:
        print(f"  [warn] unknown family {(seq_length, model, head)} -> {metrics_path}")
        return None
    row = {"family": family, "mode": mode, "seed": seed}
    for split in ("val", "test"):
        s = d.get(split, {})
        for metric in ("balanced_accuracy", "roc_auc", "f1", "loss"):
            row[f"{split}_{metric}"] = s.get(metric, np.nan)
    return row


def load_all() -> pd.DataFrame:
    rows = []
    for p in SWEEP_ROOT.rglob("metrics.json"):
        r = _load_one(p)
        if r is not None:
            rows.append(r)
    if not rows:
        raise FileNotFoundError(f"No metrics.json found under {SWEEP_ROOT}")
    df = pd.DataFrame(rows)
    return df


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [c for c in df.columns if c.startswith(("val_", "test_"))]
    agg = (df.groupby(["family", "mode"])[metric_cols]
              .agg(["mean", "std", "count"])
              .reset_index())
    # Flatten MultiIndex cols
    agg.columns = ["_".join(c).rstrip("_") for c in agg.columns]
    # Stable mode ordering
    mode_rank = {m: i for i, m in enumerate(MODE_ORDER)}
    agg["__mode_rank"] = agg["mode"].map(mode_rank).fillna(99)
    agg = agg.sort_values(["family", "__mode_rank"]).drop(columns="__mode_rank")
    return agg.reset_index(drop=True)


def _fmt_pm(mean: float, std: float, count: int) -> str:
    if np.isnan(mean):
        return "—"
    if count == 1 or np.isnan(std):
        return f"{mean:.3f}"
    return f"{mean:.3f} ± {std:.3f}"


# Family display order (best-val-BalAcc family first)
def _family_block_order(agg: pd.DataFrame) -> list:
    fam_best = (agg.groupby("family")["val_balanced_accuracy_mean"]
                   .max().sort_values(ascending=False))
    return fam_best.index.tolist()


def render_png(agg: pd.DataFrame, out_path: Path):
    # Sort within each family by val BalAcc desc; families ordered by best val
    fam_order = _family_block_order(agg)
    blocks = []
    for fam in fam_order:
        block = (agg[agg["family"] == fam]
                  .sort_values("val_balanced_accuracy_mean", ascending=False))
        blocks.append(block)
    plot_df = pd.concat(blocks, ignore_index=True)

    # Find best value per metric column for bolding
    metric_keys = [
        ("val_balanced_accuracy", "max"),
        ("val_roc_auc",           "max"),
        ("val_f1",                "max"),
        ("test_balanced_accuracy", "max"),
        ("test_roc_auc",          "max"),
        ("test_f1",               "max"),
    ]
    best_idx = {}
    for key, mode in metric_keys:
        col = f"{key}_mean"
        if mode == "max":
            best_idx[key] = plot_df[col].idxmax()
        else:
            best_idx[key] = plot_df[col].idxmin()

    header = ["Family", "Mode",
              "Val bACC", "Val AUC", "Val F1",
              "Test bACC", "Test AUC", "Test F1"]
    cells = []
    row_idx_seq = []   # parallel list — index of (family, mode) in plot_df
    for idx, r in plot_df.iterrows():
        n = int(r["val_balanced_accuracy_count"]) if not np.isnan(r["val_balanced_accuracy_count"]) else 0
        row = [
            r["family"],
            r["mode"],
            _fmt_pm(r["val_balanced_accuracy_mean"], r["val_balanced_accuracy_std"], n),
            _fmt_pm(r["val_roc_auc_mean"],          r["val_roc_auc_std"],          n),
            _fmt_pm(r["val_f1_mean"],               r["val_f1_std"],               n),
            _fmt_pm(r["test_balanced_accuracy_mean"], r["test_balanced_accuracy_std"], n),
            _fmt_pm(r["test_roc_auc_mean"],         r["test_roc_auc_std"],         n),
            _fmt_pm(r["test_f1_mean"],              r["test_f1_std"],              n),
        ]
        cells.append(row)
        row_idx_seq.append(idx)

    # Family group separators: dashed line between family blocks.
    # block_boundaries holds plot_df row indices of the FIRST row of each
    # non-leading block. cells_dict uses (table_row, col) where table_row=0
    # is the header and table_row=k is plot_df row k-1, so the corresponding
    # cells_dict row is plot_df_idx + 1.
    family_seq = plot_df["family"].tolist()
    block_boundaries = []  # plot_df indices of next-block-first-row
    for i in range(len(family_seq) - 1):
        if family_seq[i] != family_seq[i + 1]:
            block_boundaries.append(i + 1)

    # Map column key -> table column index (2-based since cols 0,1 are family/mode)
    col_idx = {
        "val_balanced_accuracy":  2, "val_roc_auc":  3, "val_f1": 4,
        "test_balanced_accuracy": 5, "test_roc_auc": 6, "test_f1": 7,
    }

    n_rows = len(cells)
    fig_h = 0.40 * (n_rows + 1) + 1.4
    fig, ax = plt.subplots(figsize=(15.5, fig_h))
    ax.set_axis_off()

    # Column widths (relative; sum to ~1)
    col_widths = [0.22, 0.30,   # family + mode (wide)
                   0.080, 0.080, 0.080,   # val 3
                   0.080, 0.080, 0.080]   # test 3
    tbl = ax.table(cellText=cells, colLabels=header,
                    colWidths=col_widths,
                    loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.35)

    # Style every cell — match T1d reference (white bg, italic header, hide
    # default cell borders so we can draw custom dashed inner separators).
    for (i, j), cell in tbl.get_celld().items():
        cell.set_facecolor("#ffffff")
        cell.set_edgecolor("#ffffff")  # hide default solid borders
        if i == 0:
            cell.set_text_props(style="italic", weight="normal",
                                  color="#222")
            cell.set_height(cell.get_height() * 1.05)
        else:
            cell.set_text_props(color="#222")
            # Left-align family + mode columns; centre numeric columns
            if j in (0, 1):
                cell.set_text_props(ha="left",
                                      style="italic" if j == 1 else "normal")
                cell._loc = "left"
                # Padding-left via x-margin (matplotlib quirk: do via PAD)
                cell.PAD = 0.04

    # Bold the global-best cell per metric column
    for key, idx in best_idx.items():
        # idx is the row-index in plot_df; the table rows are in the same order
        # as plot_df (we iterated with iterrows over its sorted form).
        # Map plot_df index back to table row position
        try:
            tbl_row = plot_df.index.get_loc(idx) + 1   # +1 for header
        except KeyError:
            continue
        col = col_idx[key]
        c = tbl.get_celld()[(tbl_row, col)]
        c.set_text_props(weight="bold")

    # --- custom borders (drawn over the table) ---
    # Get bbox of header + body to draw outer frame + lines
    fig.canvas.draw()
    cells_dict = tbl.get_celld()
    # outer corners
    top_left = cells_dict[(0, 0)]
    top_right = cells_dict[(0, len(header) - 1)]
    bot_left = cells_dict[(n_rows, 0)]
    bot_right = cells_dict[(n_rows, len(header) - 1)]
    x0 = top_left.get_x()
    x1 = top_right.get_x() + top_right.get_width()
    y_top = top_left.get_y() + top_left.get_height()
    y_below_header = top_left.get_y()
    y_bot = bot_left.get_y()

    # Outer frame
    rect = plt.Rectangle((x0, y_bot), x1 - x0, y_top - y_bot,
                          fill=False, edgecolor="#222", lw=1.2,
                          transform=ax.transAxes, zorder=5)
    ax.add_patch(rect)
    # Solid line below header
    ax.plot([x0, x1], [y_below_header, y_below_header],
            color="#222", lw=1.0, transform=ax.transAxes, zorder=6)

    # Dashed family separators — line at TOP of next-block's first cell,
    # which equals the BOTTOM of the previous block's last cell (the
    # boundary between the two blocks).
    for r_pdf in block_boundaries:
        cell = cells_dict[(r_pdf + 1, 0)]   # +1 to skip header row
        y_sep = cell.get_y() + cell.get_height()
        ax.plot([x0, x1], [y_sep, y_sep], color="#666", lw=0.8,
                 linestyle=(0, (4, 3)), transform=ax.transAxes, zorder=6)

    # Dashed vertical inter-column separators
    for j in range(len(header) - 1):
        cell = cells_dict[(0, j)]
        x_sep = cell.get_x() + cell.get_width()
        ax.plot([x_sep, x_sep], [y_bot, y_top], color="#888", lw=0.6,
                 linestyle=(0, (3, 3)), transform=ax.transAxes, zorder=6)

    # Title (3 lines, T1d-style) — placed close to table top
    title_text = (
        "Corrected-GTEx AlphaGenome func-integration sweep "
        "(ADNI ad_vs_cn binary)\n"
        "mean ± std across seeds 0/1/2  ·  ranked by Val balanced accuracy "
        "within each family\n"
        "Families: BMFM-SNP 2L MLP (101bp, chrom_hier, w=256, d=0.1, lr=0.001)  ·  "
        "NTv2 2L MLP (101bp, chrom_hier, w=512, d=0.3, lr=0.0003)  ·  "
        "BMFM-SNP XGB (1001bp, global_attn, lr=0.05, depth=4, n_est=20)")
    n_train_text = ("N per seed (typical): train 265-270 (pos 146-150 / neg 119-120) "
                    "·  val 30-32 (pos 16-17 / neg 14-15) "
                    "·  test 34-37 (pos 16-19 / neg 14-18)")
    # Place title + subtitle just above the table (in axes coords)
    ax.text(0.5, y_top + 0.060, title_text,
             ha="center", va="bottom", transform=ax.transAxes,
             fontsize=11, color="#222")
    ax.text(0.5, y_top + 0.012, n_train_text,
             ha="center", va="bottom", transform=ax.transAxes,
             fontsize=9, style="italic", color="#444")

    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_all()
    print(f"[load] {len(df)} runs from {SWEEP_ROOT}")
    print(f"[load] families: {df['family'].value_counts().to_dict()}")
    print(f"[load] modes:    {df['mode'].value_counts().to_dict()}")
    print(f"[load] seeds:    {sorted(df['seed'].unique())}")

    # Per-seed long-form CSV
    df.to_csv(OUT_DIR / "sweep_table_per_seed.csv", index=False)

    # Aggregated summary CSV
    agg = aggregate(df)
    agg.to_csv(OUT_DIR / "sweep_table_summary.csv", index=False)

    # Wide PNG sorted by val BalAcc desc
    png_path = OUT_DIR / "sweep_table.png"
    render_png(agg, png_path)

    n_complete = len(agg[agg["val_balanced_accuracy_count"] == 3])
    print(f"\n[summary] {len(agg)} (family, mode) groups, "
          f"{n_complete} fully complete (3 seeds), "
          f"{len(agg) - n_complete} partial")
    print(f"\n[out] {OUT_DIR / 'sweep_table_per_seed.csv'}")
    print(f"[out] {OUT_DIR / 'sweep_table_summary.csv'}")
    print(f"[out] {png_path}")


if __name__ == "__main__":
    main()
