"""
07_render_prognostic_conversion_tables.py — T3/T4 + T1e cross-model tables
===========================================================================
Companion to `06_render_cross_model_table.py`. 06 renders the DIAGNOSTIC
table (T1/T1b/T1c/T1d/T2); this renders the conversion/prognosis tasks into
the SAME output dir (`mri_pipeline/outputs/cross_model/`) WITHOUT touching
06, reusing 06's machinery (run collection, augment back-fill, cached-head
HP-winner selection, mean±std summarise, degeneracy flagging, cached->parent
relabel, `_row_label`).

Layout mirrors the diagnostic table (rows = model/variant/augment), but each
TASK is a column GROUP with metric SUB-COLUMNS: test bACC | val bACC | AUC |
macro-F1. One PNG + one long-form CSV per group.
  • cross_model_T3T4 — T3a/T3b/T3c (conv to AD ≤3/5/7y, binary) + T4 (3-class
    horizon). T3d/10y excluded (degenerate). Frozen/cached rows first.
  • cross_model_T1e — pCN vs sCN (binary), separate smaller table.

Both-layouts collector: the BrainMVP T4 full_ft results live ONLY in the
double-nested `brainmvp_debug/brainmvp_debug/aug_*/BrainMVP_uniformer/
T4_conv_horizon/...`; 06's collector would miss them (its double-nested
fallback only fires when the flat glob found nothing, and flat BrainMVP runs
exist). We glob BOTH layouts per tree and union (dedup on path).

Metric mapping (via 06->05 `read_run`): `auc` = `auc_roc` (binary) /
`auc_roc_ovr` (T4); `f1` = `f1` (binary) / `macro_f1` (T4). Only val bACC is
stored (`config.best_val_balanced_acc`); there is no val AUC/F1.

Re-run after the HPC T4 cached sweep lands locally — the cached MODEL_TREES
globs use `*` in the task slot, so T4_conv_horizon is picked up automatically
and appears as `BrainMVP / frozen / none†` and `BrainDINO / frozen / none†`
rows (their T4 cells fill in) alongside the BrainMVP full_ft rows.

Usage
-----
  python mri_pipeline/07_render_prognostic_conversion_tables.py
  python mri_pipeline/07_render_prognostic_conversion_tables.py --root <derivatives> --out <dir>
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Import 06 as a module (filename starts with a digit). Importing runs 06's
# module-level setup (utf-8 stdout rebind + import of 05). We reuse 06's
# functions and read its HP_TUNED_MARKER / _row_label.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_06_PATH = os.path.join(_HERE, "06_render_cross_model_table.py")
_spec = importlib.util.spec_from_file_location("_06_cross_model", _06_PATH)
m06 = importlib.util.module_from_spec(_spec)
sys.modules["_06_cross_model"] = m06
_spec.loader.exec_module(m06)


# ---------------------------------------------------------------------------
# Task groups
# ---------------------------------------------------------------------------
T3T4_TASKS = ["T3a_conv3y", "T3b_conv5y", "T3c_conv7y", "T4_conv_horizon"]
T3T4_LABELS = {
    "T3a_conv3y":      "T3a: conv ≤3y",
    "T3b_conv5y":      "T3b: conv ≤5y",
    "T3c_conv7y":      "T3c: conv ≤7y",
    "T4_conv_horizon": "T4: horizon <3/3-7/≥7y",
}
T3T4_CHANCE = {"T3a_conv3y": 0.5, "T3b_conv5y": 0.5, "T3c_conv7y": 0.5,
               "T4_conv_horizon": 1 / 3}

T1E_TASKS = ["T1e_pcn_vs_scn"]
T1E_LABELS = {"T1e_pcn_vs_scn": "T1e: pCN vs sCN"}
T1E_CHANCE = {"T1e_pcn_vs_scn": 0.5}

# Metric sub-columns within each task group: (summary key, header).
METRIC_SUB = [
    ("balanced_acc", "test bACC"),
    ("auc",          "test AUC"),
    ("f1",           "test F1"),
    ("val_bacc",     "val bACC"),
]
# Validation-only sub-columns (early-stop metric + val AUC/F1 at the best epoch;
# read_run pulls val_auc/val_f1 from train_log.csv / a val_metrics block).
VAL_METRIC_SUB = [
    ("val_bacc", "val bACC"),
    ("val_auc",  "val AUC"),
    ("val_f1",   "val F1"),
]

# Row ordering: frozen (cached HP-winner) rows FIRST, then full_ft; within a
# variant by encoder, then augment.
ROW_VARIANT_ORDER = {"frozen": 0, "full_ft": 1, "lora": 2}
ROW_MODEL_ORDER = {"BrainDINO": 0, "BrainMVP": 1, "ViT-MAE75": 2}
ROW_AUG_ORDER = {"none": 0, "stochastic": 1, "plus_original": 2}

# Trimmed footnote — only the models/augments that appear in these tables.
AUG_LEGEND_TRIM = [
    "Augmentation key:",
    "⁰ none          = no train-time augmentation (frozen-encoder cached forward, or --augment=none).",
    "² stochastic    = BrainMVP trainer; RandFlip×3 + RandAffine + RandGaussianNoise + intensity jitter (always-on).",
    "³ plus_original = originals retained + K=1 deterministic (p=1.0) augmented copy (inner ops = BrainMVP stochastic).",
    "† HP-tuned      = cached-head sweep (frozen encoder + linear head on cached embeddings; trainer =",
    "                  04_head_finetune_from_embeddings.py). Grid: lr ∈ {1e-3,1e-4,1e-5} × drop_rate ∈ {0.1,0.2,0.3} ×",
    "                  label_smoothing ∈ {0.0,0.1} = 18 combos per (task, seed); epochs=50, wd=1e-5, patience=15.",
    "                  HP winner per task by mean val_bacc across 3 seeds.",
]
HP_KEY_TRIM = [
    "HP key (verified vs W&B run configs):",
    "  BrainMVP   / full_ft / *      lr=5e-5, wd=1e-4, drop=0.1, ls=0.0, epochs=200, patience=50, bs=4, AdamW.",
    "  {BrainMVP, BrainDINO, ViT-MAE75} / frozen / none†   cached-head HP sweep; HP-winner per task from the † grid above.",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=m06.DEFAULT_ROOT,
                   help="Directory containing the model output trees.")
    p.add_argument("--out", default=m06.DEFAULT_OUT,
                   help="Output dir (default: <repo>/mri_pipeline/outputs/cross_model).")
    return p.parse_args()


def _collect_runs_both_layouts(root):
    """06._collect_runs, but globs BOTH the flat and the double-nested layout
    per tree and unions (dedup on absolute path) — so the double-nested
    BrainMVP T4 full_ft runs are captured even though flat BrainMVP runs exist."""
    rows, seen = [], set()
    for model, rel in m06.MODEL_TREES:
        tree = rel.split("/")[0]
        files = set(glob.glob(os.path.join(root, rel)))                 # flat
        files |= set(glob.glob(os.path.join(root, tree, rel)))          # double-nested
        n_ok = 0
        for f in sorted(files):
            key = os.path.normcase(os.path.abspath(f))
            if key in seen:
                continue
            seen.add(key)
            try:
                rows.append(m06.read_run(f, model))
                n_ok += 1
            except Exception as exc:
                print(f"  [WARN] unreadable: {f}  ({exc})")
        if n_ok:
            print(f"  [ok]   {model:18s} — {n_ok} runs ({tree}/)")
    return pd.DataFrame(rows)


def _relabel_cached_to_parent(summ_df):
    """Surface cached-head HP-winner rows under their parent encoder name as
    `<Model> / frozen / none†` (same as 06.main)."""
    for cached_name, parent_name in [
        ("BrainDINO-cached", "BrainDINO"),
        ("BrainMVP-cached",  "BrainMVP"),
        ("ViT-MAE-cached",   "ViT-MAE75"),
    ]:
        mask = summ_df["model"] == cached_name
        summ_df.loc[mask, "model"]   = parent_name
        summ_df.loc[mask, "variant"] = "frozen"
        summ_df.loc[mask, "augment"] = "none" + m06.HP_TUNED_MARKER
    return summ_df


def _aug_clean(aug):
    if isinstance(aug, str) and aug.endswith(m06.HP_TUNED_MARKER):
        return aug[:-len(m06.HP_TUNED_MARKER)]
    return aug


def _row_sort_key(model, variant, augment):
    return (ROW_VARIANT_ORDER.get(variant, 9),
            ROW_MODEL_ORDER.get(model, 9),
            ROW_AUG_ORDER.get(_aug_clean(augment), 9),
            model, variant, augment)


def _cell(m, s):
    if m is None or (isinstance(m, float) and np.isnan(m)):
        return "—"
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return f"{m:.3f}"
    return f"{m:.3f}±{s:.3f}"


def _render_grouped_table(sub, tasks, labels, chance, out_png, title,
                          metric_sub=METRIC_SUB, val_only=False):
    """One table: rows = model/variant/augment (frozen first), columns = each
    task as a group of metric sub-columns (`metric_sub`) via a two-row header
    (task labels spanned over their sub-columns by hiding inner edges).
    The first sub-column is the one bolded (best per task)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed — skipping PNG.")
        return

    M, T = len(metric_sub), len(tasks)
    best_key = metric_sub[0][0]
    n_cols = 1 + M * T + 1  # model + (task × metric) + n

    groups = sorted({(r.model, r.variant, r.augment) for r in sub.itertuples()},
                    key=lambda t: _row_sort_key(*t))

    def _get(g, task, mk):
        cell = sub[(sub["model"] == g[0]) & (sub["variant"] == g[1])
                   & (sub["augment"] == g[2]) & (sub["task"] == task)]
        if cell.empty:
            return None, None, 0, 0
        r = cell.iloc[0]
        return (r.get(f"{mk}_mean"), r.get(f"{mk}_std"),
                int(r["n_seeds"]), int(r["n_degenerate"]))

    # ── header rows ──────────────────────────────────────────────────────
    hdrA = [""]                                  # task-span row
    for task in tasks:
        grp = [""] * M
        grp[(M - 1) // 2] = labels[task]
        hdrA += grp
    hdrA += [""]
    hdrB = ["Model / variant / augment"]         # metric row
    for _ in tasks:
        hdrB += [h for _, h in metric_sub]
    hdrB += ["n"]

    body = [hdrA, hdrB]
    deg_rows = []
    for g in groups:
        line = [m06._row_label(*g)]
        n_row, d_row = 0, 0
        for task in tasks:
            for mk, _ in metric_sub:
                m, s, n, d = _get(g, task, mk)
                line.append(_cell(m, s))
                n_row = max(n_row, n)
                d_row = max(d_row, d)
        line.append(f"{n_row}" + ("*" if d_row else ""))
        body.append(line)
        deg_rows.append(d_row > 0)

    # column widths from content
    col_chars = [max(len(str(row[j])) for row in body) for j in range(n_cols)]
    col_chars[0] += 6
    col_chars = [max(c, 5) for c in col_chars]
    total = sum(col_chars)

    fig, ax = plt.subplots(figsize=(0.105 * total, 0.55 + 0.4 * len(body)))
    ax.axis("off")
    tab = ax.table(cellText=body, loc="upper center", cellLoc="center",
                   colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1.0, 1.5)

    # header row 0 (task spans): hide inner vertical edges so each task label
    # spans its metric sub-columns (borders only, no shading).
    for ti in range(T):
        s = 1 + ti * M
        for k in range(M):
            c = s + k
            cell = tab[0, c]
            cell.set_text_props(weight="bold")
            edges = "TB" + ("L" if k == 0 else "") + ("R" if k == M - 1 else "")
            cell.visible_edges = edges
    # header row 1 (metric labels + model header)
    for j in range(n_cols):
        tab[1, j].set_text_props(weight="bold")
    tab[1, 0].set_text_props(weight="bold", ha="left")

    # body: left-align model col; bold best test-bACC per task (no shading)
    n_body = len(groups)
    for i in range(2, 2 + n_body):
        tab[i, 0].set_text_props(ha="left")
    for ti, task in enumerate(tasks):            # bold best (first metric) per task
        col = 1 + ti * M  # first metric sub-col
        best_i, best_v = None, -1.0
        for gi, g in enumerate(groups):
            m, *_ = _get(g, task, best_key)
            if m is not None and not (isinstance(m, float) and np.isnan(m)) and m > best_v:
                best_v, best_i = m, gi
        if best_i is not None:
            tab[2 + best_i, col].set_text_props(weight="bold")

    plt.title(title, pad=12, fontsize=12)

    chance_txt = ", ".join(f"{labels[t].split(':')[0]} {chance[t]:.2f}" for t in tasks)
    if val_only:
        intro = [
            "VALIDATION metrics only (model-selection / early-stop set, NOT the held-out test set). "
            "val bACC = balanced accuracy; val AUC = ROC (one-vs-rest for the 3-class T4); val F1 = macro-F1.",
            "Cells = mean ± std across seeds; bold = best val bACC per task; '—' = not recorded "
            "(trainer logged only val bACC; AUC/F1 fill after the --val_test recompute / cached-head re-run).",
        ]
    else:
        intro = [
            "val bACC = validation balanced accuracy (early-stop / checkpoint metric); test AUC = ROC (one-vs-rest for the 3-class T4); test F1 = macro-F1.",
            "Cells = mean ± std across seeds; bold = best test bACC per task; '—' = model did not run that task.",
        ]
    foot = intro + [
        f"Chance bACC: {chance_txt}.",
        "n = seeds at the winning HP combo; '*' = at least one degenerate seed (single-class prediction; see CSV n_degenerate).",
        "",
    ] + AUG_LEGEND_TRIM + [""] + HP_KEY_TRIM
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, family="monospace", linespacing=1.4)

    fig.savefig(out_png, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG       : {out_png}")


def _render_group(summ_df, out_dir, stem, tasks, labels, chance, title,
                  val_title=None):
    sub = summ_df[summ_df["task"].isin(tasks)].reset_index(drop=True)
    print("\n" + "=" * 78)
    print(f"  GROUP {stem}  ({len(sub)} model/task rows)")
    print("=" * 78)
    if sub.empty:
        print("  [skip] no runs for these tasks under --root.")
        return
    m06.TASK_ORDER = list(tasks)        # for _render_csv / _print_stdout_summary
    m06.TASK_LABELS = dict(labels)
    m06.CHANCE = dict(chance)

    # Combined test+val table (existing).
    m06._render_csv(sub, os.path.join(out_dir, f"{stem}.csv"))
    _render_grouped_table(sub, tasks, labels, chance,
                          os.path.join(out_dir, f"{stem}.png"), title)

    # Validation-only table (val bACC / AUC / F1) + matching CSV.
    _render_grouped_table(sub, tasks, labels, chance,
                          os.path.join(out_dir, f"{stem}_val.png"),
                          val_title or title.replace("TEST", "VALIDATION"),
                          metric_sub=VAL_METRIC_SUB, val_only=True)
    _val_cols = ["model", "variant", "augment", "task", "n_seeds", "n_degenerate",
                 "val_bacc_mean", "val_bacc_std", "val_auc_mean", "val_auc_std",
                 "val_f1_mean", "val_f1_std"]
    sub[[c for c in _val_cols if c in sub.columns]].to_csv(
        os.path.join(out_dir, f"{stem}_val.csv"), index=False, float_format="%.4f")
    print(f"  CSV (val) : {os.path.join(out_dir, f'{stem}_val.csv')}")

    m06._print_stdout_summary(sub)


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("=" * 78)
    print("  07_render_prognostic_conversion_tables — T3/T4 + T1e")
    print(f"  Root: {args.root}")
    print(f"  Out : {args.out}")
    print("=" * 78)

    df = _collect_runs_both_layouts(args.root)
    if df.empty:
        print("\nNo runs found under any model tree — nothing to render.")
        sys.exit(1)

    df = m06._backfill_augment(df)
    df = m06._filter_cached_to_hp_winners(df)
    summ_df = m06._summarise(df)
    summ_df = _relabel_cached_to_parent(summ_df)

    t4_ft = summ_df[(summ_df["task"] == "T4_conv_horizon")
                    & (summ_df["variant"] == "full_ft")]
    print(f"\n  T4 full_ft groups captured: {len(t4_ft)} "
          f"(augments: {sorted(t4_ft['augment'].unique().tolist())})")
    if t4_ft.empty:
        print("  [WARN] no T4 full_ft rows — check the double-nested glob / scp.")

    _render_group(summ_df, args.out, "cross_model_T3T4",
                  T3T4_TASKS, T3T4_LABELS, T3T4_CHANCE,
                  "MRI prognostic/conversion tasks — T3 (conv ≤3/5/7y) + T4 (horizon)\n"
                  "cross-model results — TEST bACC / val bACC / AUC / macro-F1 (mean ± std across seeds)",
                  val_title="MRI prognostic/conversion tasks — T3 (conv ≤3/5/7y) + T4 (horizon)\n"
                  "cross-model results — VALIDATION bACC / AUC / macro-F1 (mean ± std across seeds)")
    _render_group(summ_df, args.out, "cross_model_T1e",
                  T1E_TASKS, T1E_LABELS, T1E_CHANCE,
                  "MRI conversion task T1e (pCN vs sCN) — cross-model results\n"
                  "TEST bACC / val bACC / AUC / macro-F1 (mean ± std across seeds)",
                  val_title="MRI conversion task T1e (pCN vs sCN) — cross-model results\n"
                  "VALIDATION bACC / AUC / macro-F1 (mean ± std across seeds)")

    print("\n" + "=" * 78)
    print("  Done. NOTE: T4 full_ft test sets are tiny (n~12) and prone to overfit — read the")
    print("  'val' column (≈1.0) against TEST bACC, and AUC alongside bACC. T4 cached rows fill")
    print("  in once the HPC sweep lands locally (re-run this script).")
    print("=" * 78)


if __name__ == "__main__":
    main()
