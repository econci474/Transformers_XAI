"""
05_aggregate_vit.py
===================
Aggregates per-run metrics.json files from 04_supervised_finetuning_ViT.py
into summary tables (mean +/- std across seeds 0, 1, 2), split by the
EFFECTIVE cohort (which sessions actually fed the model).

Why split by cohort
-------------------
The 04 sweep was launched with --long 1 for all 4 tasks, but the legacy
session_policy='baseline_only' on T3a/T3b silently dropped non-bl scans.
So T1/T2 trained on ~500 (bl + m12) scans while T3a/T3b trained on ~22 bl
scans only. Putting them in one table is misleading; this aggregator emits:

  vit_summary_bl.csv         - bl-only runs (effective_cohort == 'bl')
  vit_summary_bl_m12.csv     - bl + m12 runs (effective_cohort == 'bl..m12')
  vit_summary_<other>.csv    - any other cohort variants encountered
  vit_summary_all.csv        - master table with effective_cohort column
  vit_results_per_run.csv    - raw per-run rows (also has effective_cohort)

Effective cohort is derived from each metrics.json's config block via the
session_policy / long_mode / max_months / session fields. The patcher
mri_pipeline/_patch_metrics_session.py backfills these for legacy runs that
predate the schema.

Usage:
  python mri_pipeline/05_aggregate_vit.py
  python mri_pipeline/05_aggregate_vit.py --out_dir /alt/outputs/root
"""

import argparse
import json
from pathlib import Path

import pandas as pd

THIS_DIR = Path(__file__).resolve().parent

# -- CLI ----------------------------------------------------------------------
p = argparse.ArgumentParser()
p.add_argument("--out_dir", type=str,
               default=str(THIS_DIR / "outputs" / "ViT_B_mae75"),
               help="Root containing {task}/seed_{seed}/{strategy}/metrics.json")
p.add_argument("--exclude", nargs="*", default=["smoke_test"],
               help="Top-level subdirs to skip (default: smoke_test)")
args = p.parse_args()

OUT_DIR = Path(args.out_dir)


def derive_effective_cohort(cfg: dict) -> str:
    """
    Read the human-readable cohort label from cfg['session'].
    Falls back to a derivation from session_policy/long_mode/max_months for
    legacy files that predate the schema (and the patcher hasn't been run on).
    """
    s = cfg.get("session")
    if s is not None:
        return str(s)
    sp = cfg.get("session_policy", "current")
    if sp == "baseline_only":
        return "bl"
    long_mode = cfg.get("long_mode")
    max_months = cfg.get("max_months")
    if long_mode is None:
        return "bl"
    if long_mode == "all":
        return "bl..mAll"
    if long_mode == "cutoff" and max_months is not None:
        if int(max_months) == 12:
            return "bl + m12"
        return f"bl..m{int(max_months)}"
    return "unknown"


# -- Collect all metrics.json files -------------------------------------------
json_files = [jf for jf in OUT_DIR.rglob("metrics.json")
              if not any(part in args.exclude for part in jf.parts)]
print(f"Found {len(json_files)} metrics.json files under {OUT_DIR}")
if args.exclude:
    print(f"  (excluded: {args.exclude})")

if len(json_files) == 0:
    print("No results found. Run 04_supervised_finetuning_ViT.py first.")
    raise SystemExit(0)

rows = []
for jf in json_files:
    try:
        with open(jf) as f:
            data = json.load(f)
        cfg = data["config"]
        metrics = data["test_metrics"]
        cohort = derive_effective_cohort(cfg)
        row = {
            "Model":            cfg.get("model_id", "ViT_B_mae75"),
            "Task":             cfg["task"],
            "Task_Desc":        cfg.get("task_description", cfg["task"]),
            "Seed":             cfg["seed"],
            "Strategy":         cfg["strategy"],
            "Cohort":           cohort,
            "Session":          cfg.get("session"),
            "Long_mode":        cfg.get("long_mode"),
            "Max_months":       cfg.get("max_months"),
            "Session_policy":   cfg.get("session_policy"),
            "Epochs":           cfg["epochs"],
            "Best_Epoch":       cfg.get("best_epoch", ""),
            "LR":               cfg["lr"],
            "N_Train":          cfg.get("n_train", ""),
            "N_Val":             cfg.get("n_val", ""),
            "N_Test":           cfg.get("n_test", ""),
            "Path":             str(jf.relative_to(OUT_DIR)),
        }
        row.update(metrics)
        rows.append(row)
    except Exception as e:
        print(f"  [ERROR] {jf}: {e}")

raw_df = pd.DataFrame(rows)
print(f"Loaded {len(raw_df)} runs.")
print(f"Cohorts present: {sorted(raw_df['Cohort'].unique())}\n")

per_run_path = OUT_DIR / "vit_results_per_run.csv"
raw_df.to_csv(per_run_path, index=False)
print(f"Saved per-run results -> {per_run_path}")

# -- Aggregate across seeds, grouped by cohort -------------------------------
meta_cols = {"Model", "Task", "Task_Desc", "Seed", "Strategy", "Cohort",
             "Session", "Long_mode", "Max_months", "Session_policy",
             "Epochs", "Best_Epoch", "LR", "N_Train", "N_Val", "N_Test", "Path"}
metric_cols = [c for c in raw_df.columns if c not in meta_cols]


def build_summary(sub_df: pd.DataFrame) -> pd.DataFrame:
    """Group by (Model, Task, Strategy) and compute mean +/- std strings."""
    if sub_df.empty:
        return pd.DataFrame()
    group_cols = ["Model", "Task", "Task_Desc", "Strategy"]
    agg = (sub_df
           .groupby(group_cols)[metric_cols]
           .agg(["mean", "std", "count"])
           .reset_index())
    agg.columns = [
        f"{m}_{stat}" if stat in ("mean", "std", "count") else m
        for m, stat in agg.columns
    ]
    out_rows = []
    for _, row in agg.iterrows():
        entry = {
            "Model":    row["Model"],
            "Task":     row["Task_Desc"],
            "Strategy": row["Strategy"],
            "N_seeds":  int(row[f"{metric_cols[0]}_count"]) if metric_cols else 0,
        }
        for m in metric_cols:
            mu = row.get(f"{m}_mean")
            sd = row.get(f"{m}_std")
            if pd.notna(mu):
                entry[m] = f"{mu:.4f} +/- {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}"
            else:
                entry[m] = ""
        out_rows.append(entry)
    return pd.DataFrame(out_rows)


# Master table (all cohorts, with Cohort column for filtering)
all_rows = []
for cohort, sub in raw_df.groupby("Cohort"):
    summary = build_summary(sub)
    summary.insert(2, "Cohort", cohort)
    all_rows.append(summary)
master_summary = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
master_path = OUT_DIR / "vit_summary_all.csv"
master_summary.to_csv(master_path, index=False)
print(f"Saved master summary  -> {master_path}")

# Per-cohort tables (e.g. vit_summary_bl.csv, vit_summary_bl_m12.csv)
def _safe_filename(s: str) -> str:
    out = []
    for ch in s:
        if ch.isalnum():
            out.append(ch)
        elif ch in (" ", "+"):
            out.append("_")
    safe = "".join(out).strip("_")
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe or "unknown"


for cohort, sub in raw_df.groupby("Cohort"):
    summary = build_summary(sub)
    safe = _safe_filename(cohort)
    out_path = OUT_DIR / f"vit_summary_{safe}.csv"
    summary.to_csv(out_path, index=False)
    print(f"Saved cohort='{cohort}' -> {out_path}  ({len(sub)} runs)")

# -- Print to console ---------------------------------------------------------
print(f"\n{'='*80}")
print("  ViT SUMMARY  (mean +/- std across seeds, split by effective cohort)")
print(f"{'='*80}")

key_metric_candidates = ["balanced_acc", "auc_roc", "auc_roc_ovr",
                         "f1", "macro_f1", "auc_pr"]

for cohort, sub in raw_df.groupby("Cohort"):
    summary = build_summary(sub)
    print(f"\n  --- COHORT: {cohort}  ({len(sub)} runs, "
          f"{sub['Task'].nunique()} tasks, "
          f"seeds={sorted(int(s) for s in sub['Seed'].unique())}) ---")
    for task_desc in summary["Task"].unique():
        s = summary[summary["Task"] == task_desc]
        key_metrics = [c for c in key_metric_candidates if c in s.columns]
        print(f"\n  [{task_desc}]")
        print(s[["Strategy", "N_seeds"] + key_metrics].to_string(index=False))
