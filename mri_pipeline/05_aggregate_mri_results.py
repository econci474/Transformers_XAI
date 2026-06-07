"""
05_aggregate_mri_results.py — unified cross-model MRI results comparison
=======================================================================
Scans the per-run `metrics.json` of every MRI model sweep and produces one
comparison table (+ CSVs + an optional figure) for the thesis.

Models aggregated (each a separate output tree under --root):
  ViT-MAE75    — MAE-pretrained ViT-B/3D  (vit_outputs_debug/)
  ViT-scratch  — from-scratch ViT-B/3D    (vit_baseline/)
  Spasov-CNN   — vanilla / separable 3D CNN (cnn3d_outputs/)
  BrainMVP     — UniFormer foundation model (brainmvp_debug/)

Runs incrementally — a model tree that doesn't exist yet is skipped with a
note, so this is useful now and fills in as each sweep lands.

The four pipelines write a near-common `metrics.json` schema (`config` /
`test_metrics` / `test_diagnostics`, plus `test_metrics_subject` for the ViT
and CNN). Metric keys differ binary-vs-multiclass (`auc_roc` / `f1` vs
`auc_roc_ovr` / `macro_f1`) — read defensively.

Outputs (to --out, default <root>/mri_results_summary/):
  mri_results_runs.csv      — long-form, one row per run
  mri_results_summary.csv   — grouped mean +/- std over seeds
  mri_results_comparison.png — grouped bar chart (with --plot)
  + a per-task comparison table printed to stdout

Usage
-----
  python mri_pipeline/05_aggregate_mri_results.py
  python mri_pipeline/05_aggregate_mri_results.py --root D:/.../derivatives --plot
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

DEFAULT_ROOT = r"D:/ADNI_BIDS_project/derivatives"

# (model label, metrics.json glob relative to --root)
MODEL_TREES = [
    # ViT-MAE75 split between two sweeps now:
    #   - frozen rows come from vit_outputs_debug (the original llrd=0.7 sweep;
    #     llrd is a no-op when only the head trains, so identical to llrd=1.0).
    #   - full_ft rows come from vit_outputs_hi_lr (llrd=1.0, originally called
    #     the "hi_lr" sweep but actually a LLRD-off comparison). LLRD-off beats
    #     LLRD-on by +0.017 to +0.037 test bacc on every full_ft task, so the
    #     llrd=1.0 numbers are the canonical ViT-MAE75 full_ft row going forward.
    #     vit_outputs_hi_lr also fills T1d which the original sweep lacked.
    # frozen @ llrd=0.7 for aug=random. Lives in vit_outputs_debug (the
    # original sweep). LLRD is a no-op when only the head trains, so
    # llrd=0.7 vs 1.0 doesn't matter mathematically -- this row is "the
    # canonical frozen/random row" regardless of what llrd_gamma it was
    # trained at.
    ("ViT-MAE75",        "vit_outputs_debug/ViT_B_mae75/*/seed_*/frozen/metrics.json"),
    # full_ft + (frozen for non-random augs) @ llrd=1.0 sweeps. Augmentation
    # variants land in sibling subtrees so they don't collide:
    #   vit_outputs_hi_lr/ViT_B_mae75/...                    <- aug=random (legacy unsuffixed default)
    #   vit_outputs_hi_lr/aug_none/ViT_B_mae75/...           <- aug=none
    #   vit_outputs_hi_lr/aug_plus_original/ViT_B_mae75/...  <- aug=plus_original
    # full_ft globs cover all augs; frozen glob covers only the aug_*/ subtree
    # (frozen/random comes from vit_outputs_debug above, so we don't re-glob
    # the unsuffixed hi_lr tree for frozen and accidentally duplicate).
    ("ViT-MAE75",        "vit_outputs_hi_lr/ViT_B_mae75/*/seed_*/full_ft/metrics.json"),
    ("ViT-MAE75",        "vit_outputs_hi_lr/aug_*/ViT_B_mae75/*/seed_*/full_ft/metrics.json"),
    ("ViT-MAE75",        "vit_outputs_hi_lr/aug_*/ViT_B_mae75/*/seed_*/frozen/metrics.json"),
    ("ViT-scratch",      "vit_baseline/ViT_B_scratch/*/seed_*/*/metrics.json"),
    # ViT-MAE conversion-task sweep (T3a/b/c + T4): full_ft {none,random,plus_orig}
    # + frozen {random,plus_orig}. Dedicated tree so it stays separate from the
    # diagnostic debug/hi_lr provenance split. Each cell saves val+test inline.
    ("ViT-MAE75",        "vit_outputs_conv/aug_*/ViT_B_mae75/*/seed_*/*/metrics.json"),
    # NOTE: ViT-Tiny ablation and AG-MS3D rescue2 (flips+strong) were started
    # but only got 1 W&B-uploaded run each before being abandoned. Their
    # MODEL_TREES entries are commented out so the cross-model table doesn't
    # surface the lone seed-0 cells. The on-disk shims (under
    # vit_tiny_baseline/ and agms3d_outputs_rescue2/) are left in place in
    # case the sweeps get revived later; revert this commit to re-enable.
    #   ("ViT-Tiny",         "vit_tiny_baseline/aug_*/ViT_T_scratch/*/seed_*/*/metrics.json"),
    #   ("ViT-Tiny",         "vit_tiny_baseline/ViT_T_scratch/*/seed_*/*/metrics.json"),
    #   ("AG-MS3D-r2",       "agms3d_outputs_rescue2/AGMS3DCNN_vanilla_slim/*/seed_*/metrics.json"),
    ("Spasov-CNN",       "cnn3d_outputs/Spasov3DCNN_*/*/seed_*/metrics.json"),
    # Legacy AG-MS3D run (pre-rescue separable backbone, collapsed on
    # 14/15 cells) — kept in the aggregator as evidence-of-failure.
    ("AG-MS3D-sep",      "agms3d_outputs/AGMS3DCNN/*/seed_*/metrics.json"),
    # Post-rescue vanilla-backbone run (--lr 1e-3, --label_smoothing 0.1).
    ("AG-MS3D-vanilla",  "agms3d_outputs/AGMS3DCNN_vanilla/*/seed_*/metrics.json"),
    ("BrainMVP",         "brainmvp_debug/aug_*/BrainMVP_uniformer/*/seed_*/*/metrics.json"),
    # BrainDINO supervised finetune sweeps. Three globs cover (a) the legacy
    # un-nested layout (residual frozen+head runs predating the strategy-
    # scoped refactor), (b) today's LoRA sweep under lora/, and (c) today's
    # full_ft sweep under ft/. The trainer encodes strategy + augment in the
    # path so read_run can still parse the variant from each match.
    ("BrainDINO",        "braindino_outputs/aug_*/BrainDINO_vitb16_*/*/seed_*/metrics.json"),
    ("BrainDINO",        "braindino_outputs/lora/aug_*/BrainDINO_vitb16_*/*/seed_*/metrics.json"),
    ("BrainDINO",        "braindino_outputs/ft/aug_*/BrainDINO_vitb16_*/*/seed_*/metrics.json"),
    ("BrainDINO",        "braindino_outputs/frozen/aug_*/BrainDINO_vitb16_*/*/seed_*/metrics.json"),
    # Cached-embedding head sweeps (one row per pretrained encoder, all
    # under aug_none since the encoder forward is deterministic). HP
    # leaves: <task>/seed_<n>/lr<>_d<>_ls<>/metrics.json.
    # Renamed dir from aug_none/ -> aug_none_hp_tuned/ to clarify that this
    # tree is the cached-head HP sweep (18 HP combos per task/seed), distinct
    # from the on-the-fly BrainDINO/frozen/aug_none cells under frozen/.
    # See sync_brainmvp_cached_from_wandb.py if you need to update the shim
    # writer's output path to match.
    ("BrainDINO-cached", "braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached/*/seed_*/*/metrics.json"),
    # Wide T2-only frozen-head HP sweep WINNER (lr1e-03_d0.3_ls0.1_hlinear_focal_std).
    # Lives in its own tree, so the cached HP-winner filter compares it against the
    # aug_none_hp_tuned T2 cells and picks whichever has the higher mean val-bACC
    # (the wide winner, 0.63 > 0.60) — this is the HP-tuned BrainDINO T2 row used by
    # the T2 late-fusion arm. T2-only, so it can't affect other tasks' selection.
    ("BrainDINO-cached", "braindino_outputs/aug_none_T2_wide_winner/BrainDINO_vitb16_frozen_cached/*/seed_*/*/metrics.json"),
    ("BrainMVP-cached",  "brainmvp_debug/aug_none/BrainMVP_uniformer_frozen_cached/*/seed_*/*/metrics.json"),
    ("ViT-MAE-cached",   "vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached/*/seed_*/*/metrics.json"),
]

TASK_ORDER = ["T1_binary", "T1b_binary", "T1c_binary", "T1d_binary", "T2_multiclass"]
METRICS = ["balanced_acc", "auc", "subj_balanced_acc", "accuracy", "f1"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="Directory containing the model output trees.")
    p.add_argument("--out", default=None,
                   help="Output dir for the CSVs/figure "
                        "(default: <root>/mri_results_summary).")
    p.add_argument("--plot", action="store_true",
                   help="Also render mri_results_comparison.png.")
    return p.parse_args()


def _is_degenerate(diag: dict) -> bool:
    """True if some class is never predicted (a confusion-matrix column is 0)."""
    cm = (diag or {}).get("confusion_matrix")
    if not cm:
        return False
    a = np.asarray(cm)
    if a.ndim != 2 or a.shape[0] < 2:
        return False
    return any(a[:, j].sum() == 0 for j in range(a.shape[0]))


def _read_val_metrics(d: dict, run_dir: str, best_epoch) -> dict:
    """Validation AUC/F1 at the best (early-stop) epoch.

    Two sources, in order of preference:
      1. An explicit `val_metrics` block in metrics.json — written by the
         `--val_test` checkpoint-recompute mode for trainers that never logged
         val AUC/F1 during training (BrainMVP / AG-MS3D / 3D-CNN).
      2. train_log.csv at `config.best_epoch` (1-based, matches the `epoch`
         column) — the per-epoch log the ViT / BrainDINO / cached-head trainers
         write (ViT/BrainDINO carry val_auc+val_f1; cached-head carries val_auc
         only, so val_f1 stays None there until Part B re-runs it).
    Returns {"val_auc": float|None, "val_f1": float|None}."""
    out = {"val_auc": None, "val_f1": None}
    log = os.path.join(run_dir, "train_log.csv")
    if os.path.exists(log):
        try:
            tl = pd.read_csv(log)
            if "epoch" in tl.columns and len(tl):
                row = (tl[tl["epoch"] == best_epoch]
                       if best_epoch is not None else tl.iloc[0:0])
                # Fallback when best_epoch is missing/not in the log: reproduce the
                # trainer's early-stop pick (max val_bacc, tie-break lowest val_loss)
                # so val_auc/f1 are read at the same epoch the checkpoint was saved.
                if row.empty and "val_bacc" in tl.columns and tl["val_bacc"].notna().any():
                    t = tl[tl["val_bacc"].notna()].copy()
                    by = ["val_bacc", "val_loss"] if "val_loss" in t.columns else ["val_bacc"]
                    asc = [False, True] if "val_loss" in t.columns else [False]
                    row = t.sort_values(by, ascending=asc).iloc[0:1]
                if not row.empty:
                    r0 = row.iloc[0]
                    for col in ("val_auc", "val_f1"):
                        if col in tl.columns and pd.notna(r0[col]):
                            out[col] = float(r0[col])
        except Exception:
            pass
    vm = d.get("val_metrics") or {}
    if vm:
        a = vm.get("auc", vm.get("auc_roc", vm.get("auc_roc_ovr")))
        f = vm.get("f1", vm.get("macro_f1"))
        if a is not None:
            out["val_auc"] = float(a)
        if f is not None:
            out["val_f1"] = float(f)
    return out


def read_run(path: str, model: str) -> dict:
    """Parse one metrics.json into a flat record (defensive — schemas vary)."""
    with open(path) as f:
        d = json.load(f)
    cfg = d.get("config", {})
    tm = d.get("test_metrics", {}) or {}
    tms = d.get("test_metrics_subject", {}) or {}
    val_xtra = _read_val_metrics(d, os.path.dirname(path), cfg.get("best_epoch"))
    return {
        "model":    model,
        # ViT/BrainMVP report `strategy`; the CNN reports `model_kind`.
        "variant":  cfg.get("strategy") or cfg.get("model_kind") or "-",
        "augment":  cfg.get("augment") or "-",
        "task":     cfg.get("task", "?"),
        "seed":     cfg.get("seed"),
        "best_epoch":   cfg.get("best_epoch"),
        # binary -> auc_roc/f1 ; multiclass -> auc_roc_ovr/macro_f1.
        "balanced_acc": tm.get("balanced_acc"),
        "accuracy":     tm.get("accuracy"),
        "auc":          tm.get("auc_roc", tm.get("auc_roc_ovr")),
        "f1":           tm.get("f1", tm.get("macro_f1")),
        "subj_balanced_acc": tms.get("balanced_acc"),
        "n_test":       cfg.get("n_test"),
        "degenerate":   _is_degenerate(d.get("test_diagnostics", {})),
        # Used by the cross-model aggregator to pick HP winners for cached
        # head sweeps (which produce many HP-leaves per task/seed). For
        # single-HP runs this is the same as best across training.
        "best_val_bacc": cfg.get("best_val_balanced_acc"),
        # Val AUC/F1 at the best epoch (train_log.csv or a val_metrics block).
        "val_auc":      val_xtra["val_auc"],
        "val_f1":       val_xtra["val_f1"],
        "path":         path,
    }


def group_label(row, aug_varies: dict) -> str:
    """Readable comparison-unit label, e.g. 'ViT-MAE75/full_ft',
    'Spasov-CNN/vanilla', 'BrainMVP/full_ft/aug_none'. The augment is appended
    only for a model that actually varied it (BrainMVP sweeps none / stochastic
    / plus_original) — for the ViT it is one fixed setting, not a comparison
    axis, so appending it would just be clutter."""
    lab = f"{row['model']}/{row['variant']}"
    if aug_varies.get(row["model"]) and row["augment"] not in ("-", None, ""):
        lab += f"/{row['augment']}"
    return lab


def _fmt(mean, std):
    if mean is None or (isinstance(mean, float) and np.isnan(mean)):
        return "   -   "
    return f"{mean:.3f}±{std:.3f}"


def main():
    args = parse_args()
    out_dir = args.out or os.path.join(args.root, "mri_results_summary")

    print("=" * 78)
    print("  05_aggregate_mri_results — cross-model MRI comparison")
    print(f"  Root: {args.root}")
    print("=" * 78)

    rows = []
    for model, rel in MODEL_TREES:
        files = sorted(glob.glob(os.path.join(args.root, rel)))
        tree = rel.split("/")[0]
        if not files:
            print(f"  [skip] {model:13s} — no metrics.json under {tree}/")
            continue
        n_ok = 0
        for f in files:
            try:
                rows.append(read_run(f, model))
                n_ok += 1
            except Exception as exc:
                print(f"  [WARN] unreadable: {f}  ({exc})")
        print(f"  [ok]   {model:13s} — {n_ok} runs ({tree}/)")

    if not rows:
        print("\nNo runs found under any model tree — nothing to aggregate.")
        sys.exit(1)

    df = pd.DataFrame(rows)
    # augment is a real comparison axis only where a model varied it.
    aug_varies = {m: df.loc[df["model"] == m, "augment"].nunique() > 1
                  for m in df["model"].unique()}
    df["group"] = df.apply(lambda r: group_label(r, aug_varies), axis=1)
    os.makedirs(out_dir, exist_ok=True)

    # ── Long-form per-run CSV ─────────────────────────────────────────────────
    runs_csv = os.path.join(out_dir, "mri_results_runs.csv")
    df.sort_values(["model", "variant", "augment", "task", "seed"]).to_csv(
        runs_csv, index=False)

    # ── Grouped summary (mean +/- std over seeds) ─────────────────────────────
    summ = []
    for keys, sub in df.groupby(["group", "model", "variant", "augment", "task"],
                                dropna=False):
        rec = dict(zip(["group", "model", "variant", "augment", "task"], keys))
        rec["n_seeds"] = len(sub)
        rec["n_degenerate"] = int(sub["degenerate"].sum())
        for m in METRICS:
            v = pd.to_numeric(sub[m], errors="coerce").dropna().to_numpy()
            rec[f"{m}_mean"] = round(float(v.mean()), 4) if len(v) else None
            rec[f"{m}_std"] = round(float(v.std()), 4) if len(v) else None
        summ.append(rec)
    summ_df = pd.DataFrame(summ)
    summ_csv = os.path.join(out_dir, "mri_results_summary.csv")
    summ_df.to_csv(summ_csv, index=False)

    # ── Printed per-task comparison table ─────────────────────────────────────
    tasks = [t for t in TASK_ORDER if t in set(df["task"])]
    tasks += sorted(set(df["task"]) - set(tasks))
    for task in tasks:
        t = summ_df[summ_df["task"] == task].sort_values(
            "balanced_acc_mean", ascending=False, na_position="last")
        chance = 0.333 if "multiclass" in task else 0.5
        print(f"\n  {task}   (chance balanced acc = {chance})")
        print(f"    {'model / variant':32s} {'n':>2s}  {'bal_acc':>13s}  "
              f"{'auc':>13s}  {'subj_bal_acc':>13s}  flag")
        print("    " + "-" * 86)
        for _, r in t.iterrows():
            flag = f"{r['n_degenerate']} degenerate" if r["n_degenerate"] else ""
            print(f"    {r['group']:32s} {int(r['n_seeds']):>2d}  "
                  f"{_fmt(r['balanced_acc_mean'], r['balanced_acc_std']):>13s}  "
                  f"{_fmt(r['auc_mean'], r['auc_std']):>13s}  "
                  f"{_fmt(r['subj_balanced_acc_mean'], r['subj_balanced_acc_std']):>13s}"
                  f"  {flag}")

    # ── Optional figure ───────────────────────────────────────────────────────
    if args.plot:
        _render_plot(summ_df, tasks, os.path.join(out_dir, "mri_results_comparison.png"))

    print("\n" + "=" * 78)
    print(f"  {len(df)} runs across {df['group'].nunique()} model-variant(s).")
    print(f"  Long-form : {runs_csv}")
    print(f"  Summary   : {summ_csv}")
    print("=" * 78)


def _render_plot(summ_df: pd.DataFrame, tasks: list, path: str):
    """Grouped bar chart — test balanced accuracy +/- std, x = task, bars = model."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib not installed — skipping --plot.")
        return

    groups = sorted(summ_df["group"].unique())
    n_g = len(groups)
    width = 0.8 / max(n_g, 1)
    cmap = plt.get_cmap("tab20" if n_g > 10 else "tab10")

    fig, ax = plt.subplots(figsize=(max(9, len(tasks) * n_g * 0.45), 6))
    for j, g in enumerate(groups):
        means, stds = [], []
        for task in tasks:
            row = summ_df[(summ_df["group"] == g) & (summ_df["task"] == task)]
            if len(row):
                means.append(row["balanced_acc_mean"].iloc[0] or np.nan)
                stds.append(row["balanced_acc_std"].iloc[0] or 0.0)
            else:
                means.append(np.nan)
                stds.append(0.0)
        x = np.arange(len(tasks)) + (j - (n_g - 1) / 2) * width
        ax.bar(x, means, width, yerr=stds, capsize=2, label=g, color=cmap(j % 20))

    ax.set_xticks(np.arange(len(tasks)))
    ax.set_xticklabels(tasks, rotation=20, ha="right")
    ax.set_ylabel("Test balanced accuracy")
    ax.set_title("MRI models — test balanced accuracy by task (mean ± std over seeds)")
    ax.axhline(0.5, ls="--", lw=0.8, color="grey")
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Figure    : {path}")


if __name__ == "__main__":
    main()
