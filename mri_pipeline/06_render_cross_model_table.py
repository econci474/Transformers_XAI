"""
06_render_cross_model_table.py — publication-style cross-model table
=====================================================================
Builds a single per-model x per-task comparison table at
`mri_pipeline/outputs/cross_model/`, mirroring the publication style of
`outputs/ViT_B_mae75/vit_combined_table.*`.

Rows are `model / variant / augment` groups (each comparable unit is one
row); columns are the four diagnostic + one prognostic-binary tasks the
thesis compares across:
  T1_binary, T1b_binary, T1c_binary, T1d_binary, T2_multiclass.

Conversion tasks (T3a, T3b) are intentionally dropped per user request;
the existing per-task aggregator at `05_aggregate_mri_results.py` keeps
them for the full ablation.

Each cell is `bal_acc_mean ± std (n_seeds)` — Bal Acc first, plain
accuracy dropped (uninformative under class imbalance).

Outputs:
  cross_model/cross_model_table.csv  — long-form per-cell numbers
                                       (model, variant, augment, task,
                                        bal_acc_mean, bal_acc_std, auc...)
  cross_model/cross_model_table.png  — rendered table figure
  cross_model/cross_model_table.tex  — LaTeX booktabs table

Usage
-----
  python mri_pipeline/06_render_cross_model_table.py
  python mri_pipeline/06_render_cross_model_table.py --root <derivatives>
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import io
import os
import sys

import numpy as np
import pandas as pd

# Windows cp1252 stdout can't print the Unicode superscripts used in the
# augment legend. Re-bind stdout to a utf-8 buffer so the sanity-check dump
# doesn't crash. (PNG/CSV/TeX renderers always write utf-8 to disk and are
# unaffected by this.)
if hasattr(sys.stdout, "buffer"):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                      errors="replace", line_buffering=True)
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Import helpers from the existing per-task aggregator (single source of
# truth for the metrics.json schema + degeneracy logic).
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_05_PATH = os.path.join(_HERE, "05_aggregate_mri_results.py")
_spec = importlib.util.spec_from_file_location("_05_aggregate", _05_PATH)
_05 = importlib.util.module_from_spec(_spec)
sys.modules["_05_aggregate"] = _05
_spec.loader.exec_module(_05)

read_run = _05.read_run
_is_degenerate = _05._is_degenerate

DEFAULT_ROOT = r"D:/ADNI_BIDS_project/derivatives"
DEFAULT_OUT = os.path.join(
    _HERE, "outputs", "cross_model"
)

# Reuse the upstream MODEL_TREES directly — it already names the
# AG-MS3D-sep + AG-MS3D-vanilla + BrainDINO rows (added when Stage B+C
# were wired in). No local extension needed.
MODEL_TREES = list(_05.MODEL_TREES)

# Order tasks columns left-to-right, hardest-on-the-right.
TASK_ORDER = [
    "T1_binary",
    "T1b_binary",
    "T1c_binary",
    "T1d_binary",
    "T2_multiclass",
]

TASK_LABELS = {
    "T1_binary":     "T1: CN vs MCI+AD",
    "T1b_binary":    "T1b: CN+MCI vs AD",
    "T1c_binary":    "T1c: CN vs AD",
    "T1d_binary":    "T1d: pMCI vs sMCI",
    "T2_multiclass": "T2: CN / MCI / AD",
}

CHANCE = {
    "T1_binary": 0.5,  "T1b_binary": 0.5,
    "T1c_binary": 0.5, "T1d_binary": 0.5,
    "T2_multiclass": 1 / 3,
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="Directory containing the model output trees.")
    p.add_argument("--out",  default=DEFAULT_OUT,
                   help="Output dir for the CSV/PNG/TeX (default: "
                        "<repo>/mri_pipeline/outputs/cross_model).")
    return p.parse_args()


def _group_key(row) -> tuple:
    """Hashable (model, variant, augment) identifier for the comparison
    unit. Augment is only meaningful where a model actually varied it."""
    m, v, a = row["model"], row["variant"], row["augment"]
    return (m, v if v else "-", a if a not in (None, "", "-") else "-")


def _row_label(model, variant, augment) -> str:
    """Render the row header. Augment carries a Unicode-superscript marker
    (see AUG_LEGEND_LINES under the table) so a reader can decode the
    augmentation policy without checking source.

    If the augment value ends in HP_TUNED_MARKER ('†'), the row was produced
    by the cached-head HP-tuning sweep -- the marker is stripped for the
    superscript lookup and re-appended as the final label-tail token."""
    lab = f"{model} / {variant}"
    if augment in ("-", None, ""):
        return lab
    is_hp = isinstance(augment, str) and augment.endswith(HP_TUNED_MARKER)
    aug_clean = augment[:-len(HP_TUNED_MARKER)] if is_hp else augment
    sup = AUG_SUPERSCRIPT.get(aug_clean, "")
    suffix = HP_TUNED_MARKER if is_hp else ""
    return f"{lab} / {aug_clean}{sup}{suffix}"


# ---------------------------------------------------------------------------
# Augment legend: each label below maps to a concrete MONAI pipeline. The
# table cells get a Unicode-superscript marker after the aug name and a
# matching footnote block under the PNG / a caption note in the TeX so the
# reader can decode each label without having to read source.
AUG_SUPERSCRIPT = {
    "none":         "⁰",   # 0
    "random":       "¹",   # 1
    "stochastic":   "²",   # 2
    "plus_original": "³",  # 3
    "flips":        "⁴",   # 4
    "flips+strong": "⁵",   # 5
}
# Marker postfixed to an augment value in summ_df (e.g. "none†") to flag rows
# that were produced by the cached-head HP-tuning sweep rather than a single-HP
# trainer run. _row_label strips this for the superscript lookup and re-appends
# it as the final label-tail marker; the legend block below explains the sweep.
HP_TUNED_MARKER = "†"
AUG_LEGEND_LINES = [
    "Augmentation key:",
    "⁰ none = no train-time augmentation (frozen-encoder cached forward, "
    "or --augment=none).",
    "¹ random = ViT trainer default: RandFlip×3 + RandRotate90 + "
    "RandScaleIntensity + RandShiftIntensity (each fires stochastically at "
    "p≈0.2).",
    "² stochastic = BrainMVP trainer: RandFlip×3 + RandAffine + "
    "RandGaussianNoise + intensity jitter (BrainMVP-internal pipeline, "
    "always-on stochastic).",
    "³ plus_original = originals retained + K=1 deterministic augmented "
    "copy (p=1.0); the inner aug ops are the SAME as the trainer's stochastic "
    "set (ViT: random ops; BrainMVP: stochastic ops).",
    "⁴ flips = RandFlipd along 3 axes only (AG-MS3D legacy, AG-MS3D "
    "rescue1 vanilla, Spasov-CNN).",
    "⁵ flips+strong = flips + RandAffine + RandGaussianNoise + "
    "RandBiasField + RandAdjustContrast (AG-MS3D rescue2, --strong_aug).",
    "† HP-tuned via the cached-head sweep (frozen encoder + linear head "
    "on cached embeddings; trainer = 04_head_finetune_from_embeddings.py). "
    "Grid: lr ∈ {1e-3, 1e-4, 1e-5} × drop_rate ∈ {0.1, 0.2, 0.3} × "
    "label_smoothing ∈ {0.0, 0.1} = 18 HP combos per (task, seed); "
    "epochs=50, wd=1e-5, patience=15, batch=full. HP winner selected per "
    "task by mean val_bacc across 3 seeds (n=3 in the cell = 3 seeds at "
    "the winning HP combo). Non-† rows are single-HP runs from the "
    "on-the-fly trainer.",
]

# Augment back-fill: the per-pipeline trainers were written independently and
# do not all record an `augment` field in metrics.json. The values below are
# derived from each trainer's own _train_transform (verified by inspecting
# the trainer source):
#   - ViT (04_supervised_finetuning_ViT.py): --augment {none,random,plus_original};
#     "random" = RandFlip x3 + RandRotate90 + RandScaleIntensity + RandShiftIntensity,
#     all on-the-fly stochastic. Recorded explicitly in metrics.json.
#   - BrainMVP (brain_mvp/04_supervised_finetuning_BrainMVP.py): --augment
#     {none,stochastic,plus_original}. Recorded explicitly in metrics.json.
#   - BrainDINO (brain_dino/02_supervised_finetuning_BrainDINO.py): cached-head
#     pipeline used to date is encoder-frozen + deterministic forward, so no
#     train-time aug -> "none".
#   - AG-MS3D legacy + rescue1 (3d_cnn_vit/train_agms3d.py at args.strong_aug=False):
#     RandFlipd x3 (3 axes), nothing else. Not recorded in metrics.json -> "flips".
#   - AG-MS3D rescue2 (--strong_aug): RandFlip x3 + RandAffine + RandGaussianNoise
#     + RandBiasField + RandAdjustContrast. Recorded as `strong_aug: True` in
#     metrics.json -> "flips+strong".
#   - Spasov-CNN (3d_conv_net/train_3dcnn.py): RandFlipd x3 only. Not recorded
#     in metrics.json -> "flips".
# When a model adds new aug options later, extend this mapping.
_DEFAULT_AUG_BY_MODEL = {
    "AG-MS3D-sep":      "flips",
    "AG-MS3D-vanilla":  "flips",
    "AG-MS3D-r2":       "flips+strong",
    "Spasov-CNN":       "flips",
    # ViT family records `augment` in metrics.json; if for some reason it's
    # missing (very old runs), fall back to the trainer default "random".
    "ViT-MAE75":        "random",
    "ViT-Tiny":         "random",
    "ViT-scratch":      "random",
}


def _backfill_augment(df):
    """For every row whose `augment` is missing/blank/'-', substitute the
    trainer-derived default (see _DEFAULT_AUG_BY_MODEL). Makes every cross-
    model row self-documenting about the augmentation policy used."""
    if df.empty:
        return df
    df = df.copy()
    blank = df["augment"].isin([None, "", "-"]) | df["augment"].isna()
    df.loc[blank, "augment"] = df.loc[blank, "model"].map(
        lambda m: _DEFAULT_AUG_BY_MODEL.get(m, "-"))
    return df


def _filter_cached_to_hp_winners(df):
    """For models whose name ends in '-cached' (cached-head HP sweeps with
    many HP-leaves per task), filter rows to only the HP-winner per
    (model, task). Winner = HP-combo with highest mean best_val_balanced_acc
    across seeds. This makes cached rows comparable to single-HP rows in
    the cross-model table (otherwise their mean is diluted by the full
    HP grid, including known-bad cells).

    The HP combo is identified by the path leaf (`<task>/seed_<n>/<hp>/
    metrics.json`)."""
    if df.empty or "path" not in df.columns:
        return df
    cached_mask = df["model"].str.endswith("-cached")
    if not cached_mask.any():
        return df
    cached   = df[cached_mask].copy()
    other    = df[~cached_mask].copy()
    cached["hp_leaf"] = cached["path"].apply(lambda p: os.path.basename(os.path.dirname(p)))

    val = pd.to_numeric(cached["best_val_bacc"], errors="coerce")
    cached["_val"] = val

    # Mean best_val_bacc across seeds for each (model, task, hp_leaf).
    grp = (cached.groupby(["model", "task", "hp_leaf"], dropna=False)["_val"]
                 .mean().reset_index())
    # Winner per (model, task): hp_leaf with the highest mean.
    winners = grp.loc[grp.groupby(
        ["model", "task"], dropna=False)["_val"].idxmax()]
    keep = set(zip(winners["model"], winners["task"], winners["hp_leaf"]))

    cached = cached[cached.apply(
        lambda r: (r["model"], r["task"], r["hp_leaf"]) in keep, axis=1)]
    cached = cached.drop(columns=["hp_leaf", "_val"])

    n_dropped = int(cached_mask.sum()) - len(cached)
    if n_dropped:
        print(f"  HP-winner filter: kept {len(cached)} cached rows, "
              f"dropped {n_dropped} non-winning HP cells.")
    return pd.concat([other, cached], ignore_index=True)


def _collect_runs(root):
    """Read every metrics.json under the model trees and return a long
    DataFrame (one row per training run).

    Tolerates two layouts on disk:
      flat:           <root>/<tree>/...
      double-nested:  <root>/<tree>/<tree>/...
    The second arises when local downloads were rsync'd from CSD3 with
    the source dir as the leaf (e.g. `rsync -a .../cnn3d_outputs/ <root>/cnn3d_outputs/`
    sometimes producing `<root>/cnn3d_outputs/cnn3d_outputs/<contents>`)."""
    rows = []
    for model, rel in MODEL_TREES:
        tree = rel.split("/")[0]
        # Try the flat layout first.
        files = sorted(glob.glob(os.path.join(root, rel)))
        # Fallback: double-nested (rsync side-effect).
        if not files:
            files = sorted(glob.glob(os.path.join(root, tree, rel)))
        if not files:
            print(f"  [skip] {model:18s} — no metrics.json under {tree}/")
            continue
        n_ok = 0
        for f in files:
            try:
                rows.append(read_run(f, model))
                n_ok += 1
            except Exception as exc:
                print(f"  [WARN] unreadable: {f}  ({exc})")
        print(f"  [ok]   {model:18s} — {n_ok} runs ({tree}/)")
    return pd.DataFrame(rows)


def _summarise(df):
    """Mean+std over seeds, grouped by (model, variant, augment, task).
    Output is a tidy DataFrame keyed by group tuple + task.

    Augment is always preserved verbatim now (no `aug_varies` blanking) --
    every cross-model row in the publication table explicitly carries the
    augmentation policy used. Missing/blank augments are back-filled
    upstream via _backfill_augment, so by the time _summarise sees a row
    the augment column is always a real label."""
    if df.empty:
        return df

    summ = []
    for keys, sub in df.groupby(
        ["model", "variant", "augment", "task"], dropna=False
    ):
        rec = dict(zip(["model", "variant", "augment", "task"], keys))
        rec["n_seeds"] = len(sub)
        rec["n_degenerate"] = int(sub["degenerate"].sum())
        # Test metrics (from metrics.json `test_metrics` block).
        for m in ["balanced_acc", "auc", "f1", "accuracy"]:
            v = pd.to_numeric(sub[m], errors="coerce").dropna().to_numpy()
            rec[f"{m}_mean"] = float(v.mean()) if len(v) else None
            rec[f"{m}_std"]  = float(v.std())  if len(v) else None
        # Val balanced acc (from metrics.json `config.best_val_balanced_acc`,
        # the early-stop / checkpoint-selection metric the trainers track).
        v = pd.to_numeric(sub["best_val_bacc"], errors="coerce").dropna().to_numpy()
        rec["val_bacc_mean"] = float(v.mean()) if len(v) else None
        rec["val_bacc_std"]  = float(v.std())  if len(v) else None
        summ.append(rec)
    return pd.DataFrame(summ)


def _cell_str(row, metric: str = "balanced_acc"):
    """Format the metric cell shown in PNG / LaTeX: '0.784 ± 0.052 (3)'."""
    if row is None or row.empty:
        return "—"
    m, s, n = row[f"{metric}_mean"], row[f"{metric}_std"], int(row["n_seeds"])
    if m is None or np.isnan(m):
        return "—"
    return f"{m:.3f} ± {s:.3f} ({n})"


def _build_pivot(summ_df, metric: str = "balanced_acc"):
    """Pivot: rows = (model, variant, augment), cols = tasks.
    Returns the pivot DataFrame whose cells are pre-formatted strings,
    plus a parallel DataFrame of numeric metric means (for highlighting).

    `metric` selects which mean/std columns the pivot uses --
    e.g. 'balanced_acc' for test bacc or 'val_bacc' for the early-stop /
    checkpoint-selection metric."""
    # Establish a stable row ordering: model (paper-quality first), then
    # variant, then augment. The order list is fixed; new models can be
    # appended without disturbing the existing rows.
    # Row ordering (user-specified): BrainDINO -> BrainMVP -> ViT-MAE75 ->
    # AG-MS3D -> Spasov-CNN. Cached-head sweeps sit immediately after their
    # parent model; ViT-Tiny / ViT-scratch sit under the ViT block (they
    # share the ViT-B/16 architecture, just different sizes / inits).
    MODEL_RANK = {
        "BrainDINO":         0,
        "BrainDINO-cached":  1,
        "BrainMVP":          2,
        "BrainMVP-cached":   3,
        "ViT-MAE75":         4,
        "ViT-MAE-cached":    5,
        "ViT-Tiny":          6,
        "ViT-scratch":       7,
        "AG-MS3D-r2":        8,
        "AG-MS3D-vanilla":   9,
        "AG-MS3D-sep":      10,
        "Spasov-CNN":       11,
    }
    VARIANT_RANK = {
        "full_ft": 0, "lora": 1, "frozen": 2,
        "vanilla": 3, "separable": 4,
        "scratch": 5, "agms3d": 6, "-": 9,
    }
    AUG_RANK = {
        "plus_original": 0, "stochastic": 1, "random": 1, "none": 2, "-": 9,
    }

    def _sort_key(t):
        m, v, a = t
        return (
            MODEL_RANK.get(m, 99),
            VARIANT_RANK.get(v, 99),
            AUG_RANK.get(a, 99),
            m, v, a,
        )

    groups = sorted(
        {(r.model, r.variant, r.augment) for r in summ_df.itertuples()},
        key=_sort_key,
    )

    fmt = pd.DataFrame(index=pd.MultiIndex.from_tuples(
        groups, names=["model", "variant", "augment"]
    ), columns=TASK_ORDER, dtype=object)
    num = pd.DataFrame(index=fmt.index, columns=TASK_ORDER, dtype=float)

    for (m, v, a), sub in summ_df.groupby(["model", "variant", "augment"]):
        for task in TASK_ORDER:
            cell = sub[sub["task"] == task]
            if cell.empty:
                fmt.at[(m, v, a), task] = "—"
                num.at[(m, v, a), task] = np.nan
            else:
                r = cell.iloc[0]
                fmt.at[(m, v, a), task] = _cell_str(r, metric=metric)
                num.at[(m, v, a), task] = r[f"{metric}_mean"] or np.nan

    return fmt, num


def _render_csv(summ_df, path):
    """Long-form CSV with all metrics (test + val), one row per (model,
    variant, aug, task)."""
    cols = [
        "model", "variant", "augment", "task", "n_seeds", "n_degenerate",
        "balanced_acc_mean", "balanced_acc_std",
        "auc_mean", "auc_std",
        "f1_mean", "f1_std",
        "accuracy_mean", "accuracy_std",
        "val_bacc_mean", "val_bacc_std",
    ]
    summ_df[cols].to_csv(path, index=False, float_format="%.4f")
    print(f"  CSV       : {path}")


def _render_png(fmt, num, path, title_metric: str = "test balanced accuracy"):
    """Render a matplotlib `ax.table` mirroring the look of
    `outputs/ViT_B_mae75/vit_combined_table.png`. Best-per-column cells
    are bolded. `title_metric` shows which metric the cells represent."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed — skipping PNG.")
        return

    n_rows = len(fmt)
    n_cols = len(TASK_ORDER) + 1  # +1 for the "Model / variant / augment" col

    fig, ax = plt.subplots(
        figsize=(2.2 + 2.4 * len(TASK_ORDER), 0.55 + 0.45 * (n_rows + 2))
    )
    ax.axis("off")

    # ── Header rows ─────────────────────────────────────────────────────
    col_labels = ["Model / variant / augment"] + [TASK_LABELS[t] for t in TASK_ORDER]

    # ── Body rows ───────────────────────────────────────────────────────
    body = []
    for (m, v, a), row in fmt.iterrows():
        row_header = _row_label(m, v, a)
        body.append([row_header] + [row[t] for t in TASK_ORDER])

    tab = ax.table(
        cellText=body,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1.0, 1.35)

    # Bold the best (highest bal_acc) cell in each task column
    for j, task in enumerate(TASK_ORDER, start=1):
        col_vals = num[task].to_numpy()
        if np.all(np.isnan(col_vals)):
            continue
        best_idx = int(np.nanargmax(col_vals))
        # +1 because mpl tables index row 0 as the header
        cell = tab[best_idx + 1, j]
        cell.set_text_props(weight="bold")

    # Style header row
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC")
        tab[0, j].set_text_props(weight="bold")

    # Left-align the row-header column
    for i in range(n_rows + 1):
        tab[i, 0].set_text_props(ha="left")
        tab[i, 0].PAD = 0.04

    plt.title(f"MRI models — {title_metric} by task "
              "(mean ± std across seeds, n)", pad=14, fontsize=11)

    # ── Footnote: augment legend ────────────────────────────────────────
    # Placed below the table via fig.text so it doesn't compete with the
    # cell grid. Reserves room at the bottom of the figure via subplots_
    # adjust + figsize bump so bbox_inches='tight' doesn't clip it.
    legend_text = "\n".join(AUG_LEGEND_LINES)
    fig.text(0.02, -0.005, legend_text, ha="left", va="top",
             fontsize=7, family="monospace",
             linespacing=1.4)

    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  PNG       : {path}")


def _render_tex(fmt, num, path, title_metric: str = "test balanced accuracy"):
    """Booktabs LaTeX table — drop into the thesis source tree directly."""
    n_tasks = len(TASK_ORDER)
    col_spec = "l" + "c" * n_tasks

    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{MRI models -- " + title_metric + r" by task "
        r"(mean $\pm$ std across seeds; n in parentheses). Best per task "
        r"in \textbf{bold}.}",
        r"\label{tab:mri_cross_model}",
        r"\small\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{" + col_spec + r"}",
        r"\toprule",
        "Model / variant / aug & "
        + " & ".join(
            r"\textbf{" + TASK_LABELS[t].replace("&", r"\&") + r"}"
            for t in TASK_ORDER
        ) + r" \\",
        r"\midrule",
    ]

    # Identify column maxima for bold rendering
    bests = {}
    for t in TASK_ORDER:
        col_vals = num[t].to_numpy()
        bests[t] = int(np.nanargmax(col_vals)) if not np.all(np.isnan(col_vals)) else None

    for i, ((m, v, a), row) in enumerate(fmt.iterrows()):
        cells = []
        for t in TASK_ORDER:
            raw = row[t]
            txt = raw.replace("±", r"$\pm$") if isinstance(raw, str) else "—"
            if bests[t] == i:
                txt = r"\textbf{" + txt + "}"
            cells.append(txt)
        header = _row_label(m, v, a).replace("&", r"\&")
        lines.append(header + " & " + " & ".join(cells) + r" \\")

    # ── Augment-legend footnote (matches PNG) ────────────────────────────
    # Each line is a separate \footnotesize row beneath the table; we use
    # text-mode superscripts via \textsuperscript so it renders without
    # math mode and matches the Unicode markers in the cell labels.
    def _tex_escape(s: str) -> str:
        return (s.replace("\\", r"\textbackslash{}").replace("&", r"\&")
                 .replace("%", r"\%").replace("_", r"\_")
                 .replace("≈", r"$\approx$").replace("×", r"$\times$"))
    sup_to_idx = {v: k for k, v in AUG_SUPERSCRIPT.items()}  # for ref only
    legend_tex_lines = [
        r"\multicolumn{" + str(n_tasks + 1) + r"}{l}{\footnotesize "
        + _tex_escape(line) + r"} \\"
        for line in AUG_LEGEND_LINES
    ]

    lines += [
        r"\bottomrule",
    ] + legend_tex_lines + [
        r"\end{tabular}",
        r"\end{table}",
    ]

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  TeX       : {path}")


def _print_stdout_summary(summ_df):
    """Quick per-task ranked stdout dump for sanity-checking the table."""
    print("\n" + "=" * 78)
    print("  Per-task ranking (sanity check)")
    print("=" * 78)
    for task in TASK_ORDER:
        t = summ_df[summ_df["task"] == task].sort_values(
            "balanced_acc_mean", ascending=False, na_position="last"
        )
        chance = CHANCE[task]
        print(f"\n  {task}   (chance bal_acc = {chance:.3f})")
        print(f"    {'model / variant / aug':40s} {'n':>2s}  "
              f"{'bal_acc':>13s}  {'auc':>13s}  flag")
        print("    " + "-" * 78)
        for _, r in t.iterrows():
            flag = (f"{int(r['n_degenerate'])} degenerate"
                    if r["n_degenerate"] else "")
            ba = (f"{r['balanced_acc_mean']:.3f}±{r['balanced_acc_std']:.3f}"
                  if pd.notna(r['balanced_acc_mean']) else "    -    ")
            au = (f"{r['auc_mean']:.3f}±{r['auc_std']:.3f}"
                  if pd.notna(r['auc_mean']) else "    -    ")
            label = _row_label(r['model'], r['variant'], r['augment'])
            print(f"    {label:40s} {int(r['n_seeds']):>2d}  "
                  f"{ba:>13s}  {au:>13s}  {flag}")


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("=" * 78)
    print("  06_render_cross_model_table — publication table")
    print(f"  Root: {args.root}")
    print(f"  Out : {args.out}")
    print("=" * 78)

    df = _collect_runs(args.root)
    if df.empty:
        print("\nNo runs found under any model tree — nothing to render.")
        sys.exit(1)

    df = _backfill_augment(df)
    df = _filter_cached_to_hp_winners(df)

    summ_df = _summarise(df)
    summ_df = summ_df[summ_df["task"].isin(TASK_ORDER)].reset_index(drop=True)

    # ── Relabel cached-head HP-tuned rows under their parent model name ──
    # All three (BrainDINO, BrainMVP, ViT-MAE) ran the same cached-head HP
    # sweep (04_head_finetune_from_embeddings.py); the result for each model
    # is one HP-winner row per task. Surface these as
    # "<Model> / frozen / none†" so a reader sees a single canonical row
    # per pretrained encoder and the † marker + legend explain that it's
    # the HP-winner from a sweep (not a single-HP run). Then drop the
    # single-HP "BrainMVP / frozen / none" rows (no T1c anyway; the
    # HP-tuned row supersedes them for the publication table).
    for cached_name, parent_name in [
        ("BrainDINO-cached", "BrainDINO"),
        ("BrainMVP-cached",  "BrainMVP"),
        ("ViT-MAE-cached",   "ViT-MAE75"),
    ]:
        mask = summ_df["model"] == cached_name
        summ_df.loc[mask, "model"]   = parent_name
        summ_df.loc[mask, "variant"] = "frozen"
        summ_df.loc[mask, "augment"] = "none" + HP_TUNED_MARKER

    # Drop the on-the-fly single-HP BrainMVP / frozen / none rows; the
    # HP-tuned cached variant displaces them in the publication table.
    summ_df = summ_df[~(
        (summ_df["model"] == "BrainMVP") &
        (summ_df["variant"] == "frozen") &
        (summ_df["augment"] == "none")
    )].reset_index(drop=True)

    # Merge the ViT-from-scratch family under one model name so the two sizes
    # appear as sibling sub-rows in the table (apples-to-apples size ablation):
    #   ViT-scratch / baseline / <aug>   <- existing ViT-B from-scratch sweep
    #   ViT-scratch / tiny     / <aug>   <- new ViT-Ti ablation (vit_tiny_baseline/)
    _vit_b_scratch = (summ_df["model"] == "ViT-scratch")
    summ_df.loc[_vit_b_scratch, "variant"] = "baseline"
    _vit_t_scratch = (summ_df["model"] == "ViT-Tiny")
    summ_df.loc[_vit_t_scratch, "model"]   = "ViT-scratch"
    summ_df.loc[_vit_t_scratch, "variant"] = "tiny"

    # ── Test metrics (legacy default) ─────────────────────────────────────
    fmt_test, num_test = _build_pivot(summ_df, metric="balanced_acc")
    _render_csv(summ_df, os.path.join(args.out, "cross_model_table.csv"))
    _render_png(fmt_test, num_test,
                os.path.join(args.out, "cross_model_table.png"),
                title_metric="test balanced accuracy")
    _render_tex(fmt_test, num_test,
                os.path.join(args.out, "cross_model_table.tex"),
                title_metric="test balanced accuracy")

    # ── Val metrics (early-stop / checkpoint-selection metric) ───────────
    fmt_val, num_val = _build_pivot(summ_df, metric="val_bacc")
    _render_png(fmt_val, num_val,
                os.path.join(args.out, "cross_model_table_val.png"),
                title_metric="val balanced accuracy (early-stop metric)")
    _render_tex(fmt_val, num_val,
                os.path.join(args.out, "cross_model_table_val.tex"),
                title_metric="val balanced accuracy (early-stop metric)")

    _print_stdout_summary(summ_df)

    print("\n" + "=" * 78)
    print(f"  {len(df)} runs across "
          f"{summ_df.groupby(['model','variant','augment']).ngroups} "
          f"model/variant/aug groups.")
    print("=" * 78)


if __name__ == "__main__":
    main()
