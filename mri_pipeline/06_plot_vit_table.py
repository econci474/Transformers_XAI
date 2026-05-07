"""
06_plot_vit_table.py
====================
Renders a publication-quality table of ViT_B_mae75 test performance,
mirroring the style of clinical_pipeline/04b_plot_combined_table.py.

Layout
------
  Title row
  Group-header row    (4 task groups: T1 binary, T2 multiclass, T3a, T3b)
  Sub-header row      (Acc, Bal.Acc, AUC, F1/MacroF1)
  Data rows grouped by Cohort, then Strategy:
      bl         | frozen
      bl         | full fine-tune
      bl + m12   | frozen
      bl + m12   | full fine-tune

Cells with no completed runs show '---'. Within each metric column the best
mean is bolded.

Outputs (all under outputs/ViT_B_mae75/):
  vit_combined_table.png
  vit_combined_table.pdf
  vit_combined_table.tex

Usage:
  python mri_pipeline/06_plot_vit_table.py
  python mri_pipeline/06_plot_vit_table.py --out_dir /alt/outputs/root
"""

import argparse
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.rcParams["font.family"] = "DejaVu Serif"

THIS_DIR = Path(__file__).resolve().parent

# -- CLI ----------------------------------------------------------------------
p = argparse.ArgumentParser()
p.add_argument("--out_dir", type=str,
               default=str(THIS_DIR / "outputs"),
               help="Root scanned recursively for metrics.json. Default "
                    "covers all three sweep dirs (ViT_B_mae75/, "
                    "vit_outputs_bl/, vit_outputs_bl_m12/) and writes the "
                    "PNG/PDF/LaTeX into the same root.")
p.add_argument("--exclude", nargs="*", default=["smoke_test"],
               help="Path components to skip (default: smoke_test, since the "
                    "proper bl-only T1 sweep under vit_outputs_bl/ supersedes "
                    "the early 2-epoch smoke run).")
args = p.parse_args()

OUT_DIR = Path(args.out_dir)


def cohort_from_cfg(cfg: dict) -> str:
    s = cfg.get("session")
    if s is not None:
        return str(s)
    sp = cfg.get("session_policy", "current")
    if sp == "baseline_only":
        return "bl"
    if cfg.get("long_mode") == "cutoff" and cfg.get("max_months") == 12:
        return "bl + m12"
    return "bl"


# -- Collect runs -------------------------------------------------------------
rows = []
for jf in OUT_DIR.rglob("metrics.json"):
    if any(ex in jf.parts for ex in args.exclude):
        continue
    with open(jf) as f:
        data = json.load(f)
    cfg = data["config"]
    met = data["test_metrics"]
    n_total = (cfg.get("n_train", 0) + cfg.get("n_val", 0)
               + cfg.get("n_test", 0))
    rows.append({
        "task":     cfg["task"],
        "strategy": cfg["strategy"],
        "cohort":   cohort_from_cfg(cfg),
        "seed":     cfg["seed"],
        "n_total":  n_total,
        **met,
    })
raw_df = pd.DataFrame(rows)
if raw_df.empty:
    raise SystemExit(f"No runs found under {OUT_DIR}")

# Per-cohort representative N (largest n_train+n_val+n_test seen across runs).
# T1/T2 don't filter AD subjects, so they yield the cohort-wide scan count.
cohort_N = {c: int(raw_df.loc[raw_df["cohort"] == c, "n_total"].max())
            for c in raw_df["cohort"].unique()}


def cohort_class_breakdown():
    """For each cohort, count CN/MCI/AD by reading a T2_multiclass manifest.
    The bl+m12 T2 manifest contains both bl and m12 rows, so we can derive
    both cohort breakdowns from a single manifest file. Returns
    {cohort -> {0: n_CN, 1: n_MCI, 2: n_AD}} or {} if no manifest is found."""
    out = {}
    candidates = []
    for jf in OUT_DIR.rglob("metrics.json"):
        if "T2_multiclass" not in jf.parts:
            continue
        with open(jf) as f:
            cfg = json.load(f)["config"]
        candidates.append((cohort_from_cfg(cfg),
                           jf.parent / "dataset_manifest.csv"))
    if not candidates:
        return out
    # Prefer bl + m12 (superset cohort) so we can derive both bl and bl+m12
    candidates.sort(key=lambda x: 0 if x[0] == "bl + m12" else 1)
    _, manifest_path = candidates[0]
    if not manifest_path.exists():
        return out
    m = pd.read_csv(manifest_path)
    # Drop duplicates so the same scan isn't double-counted across splits
    m = m.drop_duplicates(subset=["Patient_ID", "bids_ses"])
    for cohort, mask in [("bl", m["bids_ses"] == "bl"),
                         ("bl + m12", m["bids_ses"].isin(["bl", "m12"]))]:
        sub = m[mask]
        if sub.empty:
            continue
        counts = sub["label"].value_counts().to_dict()
        out[cohort] = {int(k): int(v) for k, v in counts.items()}
    return out


class_breakdown = cohort_class_breakdown()


def fmt_breakdown(cohort):
    bd = class_breakdown.get(cohort)
    if not bd:
        return None
    return f"CN={bd.get(0,0)}  MCI={bd.get(1,0)}  AD={bd.get(2,0)}"

# Aggregate mean +/- std per (cohort, strategy, task)
metric_cols = [c for c in raw_df.columns
               if c not in {"task", "strategy", "cohort", "seed"}]
summary = {}
for (cohort, strategy, task), grp in raw_df.groupby(["cohort", "strategy", "task"]):
    entry = {}
    for m in metric_cols:
        vals = grp[m].dropna()
        if len(vals) > 0:
            mu = vals.mean()
            sd = vals.std() if len(vals) > 1 else 0.0
            entry[m] = f"{mu:.4f} ± {sd:.4f}"
        else:
            entry[m] = "---"
    summary[(cohort, strategy, task)] = entry


# -- Table layout -------------------------------------------------------------
GROUPS = [
    {
        "title":   "CN vs MCI+AD (binary)",
        "task":    "T1_binary",
        "metrics": ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":   "CN / MCI / AD (3-class)",
        "task":    "T2_multiclass",
        "metrics": ["accuracy", "balanced_acc", "auc_roc_ovr", "macro_f1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "MacroF1"],
    },
    {
        "title":   "AD conversion ≤ 3 yrs",
        "task":    "T3a_conv3y",
        "metrics": ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
    },
    {
        "title":   "AD conversion ≤ 5 yrs",
        "task":    "T3b_conv5y",
        "metrics": ["accuracy", "balanced_acc", "auc_roc", "f1"],
        "headers": ["Acc", "Bal.Acc", "AUC", "F1"],
    },
]

# Row order: cohort sections in fixed order; strategies in fixed order within each cohort.
COHORT_ORDER = ["bl", "bl + m12"]
STRATEGY_ORDER = [("frozen", "frozen"), ("full_ft", "full fine-tune")]

ROW_ENTRIES = []   # list of (cohort, strategy_id, strategy_label)
for cohort in COHORT_ORDER:
    if not any(c == cohort for (c, _, _) in [(c, s, t) for (c, s, t) in summary.keys()]):
        continue
    for sid, slabel in STRATEGY_ORDER:
        ROW_ENTRIES.append((cohort, sid, slabel))

n_rows = len(ROW_ENTRIES)


# -- Pre-compute cells -------------------------------------------------------
def cell_for(cohort, strategy, group):
    entry = summary.get((cohort, strategy, group["task"]), {})
    return [entry.get(m, "---") for m in group["metrics"]]


def parse_mean(s):
    s = str(s).strip()
    if s in ("", "---", "nan", "NaN"):
        return float("-inf")
    try:
        return float(s.split("±")[0].strip())
    except ValueError:
        return float("-inf")


cell_grid = []
for cohort, sid, slabel in ROW_ENTRIES:
    row = []
    for grp in GROUPS:
        row.extend(cell_for(cohort, sid, grp))
    cell_grid.append(row)

n_total_cols = sum(len(g["headers"]) for g in GROUPS)
best_mask = [[False] * n_total_cols for _ in range(n_rows)]
for col_idx in range(n_total_cols):
    means = [parse_mean(cell_grid[r][col_idx]) for r in range(n_rows)]
    best_val = max(means)
    if best_val > float("-inf"):
        for r in range(n_rows):
            if means[r] == best_val:
                best_mask[r][col_idx] = True


# -- Sizing -------------------------------------------------------------------
COL_W        = 1.45
ROW_H        = 0.32
COHORT_COL_W = 1.55
STRAT_COL_W  = 1.10

n_data_cols = sum(len(g["headers"]) for g in GROUPS)

fig_w = 0.3 + COHORT_COL_W + STRAT_COL_W + n_data_cols * COL_W + 0.3
fig_h = 0.15 + ROW_H * (3 + n_rows) + 0.10

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.axis("off")
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)

LEFT = 0.3
col_widths = [COHORT_COL_W, STRAT_COL_W] + [COL_W] * n_data_cols
col_left = [LEFT]
for w in col_widths:
    col_left.append(col_left[-1] + w)
col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(col_widths))]
RIGHT = col_left[-1]

TOP = fig_h - 0.15
header_rows = 3
total_visual_rows = header_rows + n_rows
row_tops = [TOP - i * ROW_H for i in range(total_visual_rows + 1)]


def row_mid(i):
    return (row_tops[i] + row_tops[i + 1]) / 2


def hline(y, lw=1.0, ls="-"):
    ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
            solid_capstyle="butt", zorder=3)


# Top, below-title, below-group, below-subheader, between-cohort sections, bottom
hline(row_tops[0], lw=1.5)
hline(row_tops[1], lw=1.5)
hline(row_tops[2], lw=0.8)
hline(row_tops[3], lw=1.2)

# Section divider between cohort blocks (cohort changes -> dashed line)
prev_cohort = None
for r_idx, (cohort, _, _) in enumerate(ROW_ENTRIES):
    if prev_cohort is not None and cohort != prev_cohort:
        hline(row_tops[3 + r_idx], lw=0.6, ls="--")
    prev_cohort = cohort
hline(row_tops[total_visual_rows], lw=1.5)

# Dotted vertical separators between task groups
col_cursor = 2
for grp in GROUPS:
    x = col_left[col_cursor]
    ax.plot([x, x], [row_tops[total_visual_rows], row_tops[1]],
            color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)
    col_cursor += len(grp["headers"])

# Dotted between cohort and strategy columns
x_strat = col_left[1]
ax.plot([x_strat, x_strat], [row_tops[total_visual_rows], row_tops[2]],
        color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)
x_data = col_left[2]
ax.plot([x_data, x_data], [row_tops[total_visual_rows], row_tops[2]],
        color="black", linewidth=0.5, linestyle=(0, (3, 3)), zorder=2)

# Title
ax.text((LEFT + RIGHT) / 2, row_mid(0),
        "ViT_B_mae75: Test Performance   "
        "(mean ± std, seeds 0/1/2, 80/10/10 split)",
        ha="center", va="center", fontsize=8.5, fontweight="bold", color="black")

# Group headers
col_cursor = 2
for grp in GROUPS:
    n = len(grp["headers"])
    x0 = col_left[col_cursor]
    x1 = col_left[col_cursor + n]
    xm = (x0 + x1) / 2
    ax.text(xm, row_mid(1), grp["title"],
            ha="center", va="center", fontsize=8, fontweight="bold", color="black")
    pad = 0.05
    ax.plot([x0 + pad, x1 - pad], [row_tops[2], row_tops[2]],
            color="black", linewidth=0.7, zorder=4)
    col_cursor += n

# Sub-headers
ax.text(col_cx[0], row_mid(2), "Cohort",
        ha="center", va="center", fontsize=7.5, fontstyle="italic", color="black")
ax.text(col_cx[1], row_mid(2), "Strategy",
        ha="center", va="center", fontsize=7.5, fontstyle="italic", color="black")
col_cursor = 2
for grp in GROUPS:
    for h in grp["headers"]:
        ax.text(col_cx[col_cursor], row_mid(2), h,
                ha="center", va="center", fontsize=7.5,
                fontstyle="italic", color="black")
        col_cursor += 1

# Data rows
prev_cohort = None
for r_idx, (cohort, sid, slabel) in enumerate(ROW_ENTRIES):
    row_i = 3 + r_idx
    # Show cohort label only on first row of its block; vertically center
    # over the whole block (frozen + full_ft) and stack a class-breakdown
    # line beneath the N count so the underpowered-3-class case is obvious.
    if cohort != prev_cohort:
        n = cohort_N.get(cohort)
        n_str = f"  (N={n})" if n is not None else ""
        bd_str = fmt_breakdown(cohort)
        # Vertical center of this cohort block (count rows that share cohort)
        block_size = sum(1 for c, _, _ in ROW_ENTRIES[r_idx:] if c == cohort)
        block_top    = row_tops[row_i]
        block_bottom = row_tops[row_i + block_size]
        block_mid    = (block_top + block_bottom) / 2
        if bd_str:
            line_off = ROW_H * 0.25
            ax.text(col_left[0] + 0.06, block_mid + line_off,
                    f"{cohort}{n_str}",
                    ha="left", va="center", fontsize=7.5,
                    fontweight="bold", color="black")
            ax.text(col_left[0] + 0.06, block_mid - line_off,
                    bd_str,
                    ha="left", va="center", fontsize=6.3, color="black")
        else:
            ax.text(col_left[0] + 0.06, block_mid,
                    f"{cohort}{n_str}",
                    ha="left", va="center", fontsize=7.5,
                    fontweight="bold", color="black")
    ax.text(col_cx[1], row_mid(row_i), slabel,
            ha="center", va="center", fontsize=7, color="black")

    col_cursor = 2
    for ci, val in enumerate(cell_grid[r_idx]):
        fw = "bold" if best_mask[r_idx][ci] else "normal"
        ax.text(col_cx[col_cursor], row_mid(row_i), val,
                ha="center", va="center", fontsize=7.3,
                color="black", fontweight=fw)
        col_cursor += 1
    prev_cohort = cohort

rect = plt.Rectangle((LEFT, row_tops[total_visual_rows]),
                     RIGHT - LEFT, row_tops[0] - row_tops[total_visual_rows],
                     facecolor="none", edgecolor="black", linewidth=1.5, zorder=5)
ax.add_patch(rect)

plt.tight_layout(pad=0.1)
OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_DIR / "vit_combined_table.pdf", bbox_inches="tight", dpi=300)
fig.savefig(OUT_DIR / "vit_combined_table.png", bbox_inches="tight", dpi=300)
print(f"Saved PDF -> {OUT_DIR / 'vit_combined_table.pdf'}")
print(f"Saved PNG -> {OUT_DIR / 'vit_combined_table.png'}")
plt.close()


# -- LaTeX (booktabs) ---------------------------------------------------------
def latex_val(s):
    return str(s).replace("±", r"$\pm$").replace("---", r"---")


col_spec = "ll" + "".join(
    "@{\\hspace{3pt}}|@{\\hspace{3pt}}" + "c" * len(g["headers"])
    for g in GROUPS
)
lines = [
    r"\begin{table}[ht]",
    r"\centering",
    r"\caption{ViT-B (MAE-75) test performance on ADNI MRI "
    r"(mean $\pm$ std across seeds 0/1/2, 80/10/10 split). "
    r"Cohort `bl' uses baseline scans only; `bl + m12' uses baseline plus "
    r"month-12 scans, paired with the same baseline-anchored label.}",
    r"\label{tab:vit_combined}",
    r"\small\setlength{\tabcolsep}{4pt}",
    r"\begin{tabular}{" + col_spec + "}",
    r"\toprule",
]
grp_row = " & & " + " & ".join(
    r"\multicolumn{" + str(len(g["headers"])) + r"}{c}{\textbf{" + g["title"] + r"}}"
    for g in GROUPS
) + r" \\"
lines.append(grp_row)
cursor = 3
cmidrules = ""
for g in GROUPS:
    n = len(g["headers"])
    cmidrules += r"\cmidrule(lr){" + f"{cursor}-{cursor+n-1}" + "}"
    cursor += n
lines.append(cmidrules)
lines.append("Cohort & Strategy & " + " & ".join(
    r"\textit{" + h + r"}" for g in GROUPS for h in g["headers"]
) + r" \\")
lines.append(r"\midrule")

prev_cohort = None
for r_idx, (cohort, sid, slabel) in enumerate(ROW_ENTRIES):
    if prev_cohort is not None and cohort != prev_cohort:
        lines.append(r"\midrule")
    if cohort != prev_cohort:
        n = cohort_N.get(cohort)
        n_suffix = f" (N={n})" if n is not None else ""
        bd_str = fmt_breakdown(cohort)
        if bd_str:
            # Requires \usepackage{makecell} in the LaTeX preamble
            cohort_cell = (r"\makecell[l]{\textbf{" + cohort + n_suffix
                           + r"} \\ \footnotesize " + bd_str + r"}")
        else:
            cohort_cell = r"\textbf{" + cohort + n_suffix + r"}"
    else:
        cohort_cell = ""
    parts = [cohort_cell, slabel]
    for ci, val in enumerate(cell_grid[r_idx]):
        lv = latex_val(val)
        if best_mask[r_idx][ci]:
            lv = r"\textbf{" + lv + "}"
        parts.append(lv)
    lines.append(" & ".join(parts) + r" \\")
    prev_cohort = cohort

lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
tex_path = OUT_DIR / "vit_combined_table.tex"
tex_path.write_text("\n".join(lines), encoding="utf-8")
print(f"Saved LaTeX -> {tex_path}")
