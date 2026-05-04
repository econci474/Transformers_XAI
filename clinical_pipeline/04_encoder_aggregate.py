"""
04_encoder_aggregate.py
=======================
Aggregates per-run metrics.json files from 03_encoder_finetune.py into
a summary table (mean ± std across seeds 0, 1, 2).

Outputs:
  outputs/encoder_results_per_run.csv   — one row per run
  outputs/encoder_summary.csv           — mean ± std per (model, task, strategy)

Usage:
  python 04_encoder_aggregate.py
  python 04_encoder_aggregate.py --out_dir /path/to/outputs
  python 04_encoder_aggregate.py --include_baseline   # merges tabular baseline CSV too
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path

# ── CLI ───────────────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--out_dir", type=str,
               default=r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\clinical_pipeline\outputs")
p.add_argument("--include_baseline", action="store_true",
               help="Also merge tabular baseline results from baseline/baseline_model_comparison.csv")
args = p.parse_args()

OUT_DIR = Path(args.out_dir)

# ── Collect all metrics.json files ────────────────────────────────────────────
json_files = list(OUT_DIR.rglob("metrics.json"))
print(f"Found {len(json_files)} metrics.json files under {OUT_DIR}")

if len(json_files) == 0:
    print("No results found. Run 03_encoder_finetune.py first.")
    import sys; sys.exit(0)

rows = []
for jf in json_files:
    try:
        with open(jf) as f:
            data = json.load(f)
        cfg     = data["config"]
        metrics = data["test_metrics"]
        row = {
            "Model":    cfg["model_id"].split("/")[-1],
            "Model_ID": cfg["model_id"],
            "Task":     cfg["task"],
            "Task_Desc":cfg.get("task_description", cfg["task"]),
            "Seed":     cfg["seed"],
            "Strategy": cfg["strategy"],
            "Epochs":   cfg["epochs"],
            "Best_Epoch": cfg.get("best_epoch", ""),
            "LR":       cfg["lr"],
            "N_Train":  cfg.get("n_train", ""),
            "N_Test":   cfg.get("n_test", ""),
        }
        row.update(metrics)
        rows.append(row)
    except Exception as e:
        print(f"  [ERROR] {jf}: {e}")

raw_df = pd.DataFrame(rows)
print(f"Loaded {len(raw_df)} runs.\n")

# Save per-run results
per_run_path = OUT_DIR / "encoder_results_per_run.csv"
raw_df.to_csv(per_run_path, index=False)
print(f"Saved per-run results → {per_run_path}")

# ── Aggregate across seeds ────────────────────────────────────────────────────
meta_cols   = ["Model", "Model_ID", "Task", "Task_Desc", "Seed",
               "Strategy", "Epochs", "Best_Epoch", "LR", "N_Train", "N_Test"]
metric_cols = [c for c in raw_df.columns if c not in meta_cols]

agg = (raw_df
       .groupby(["Model", "Task", "Task_Desc", "Strategy"])[metric_cols]
       .agg(["mean", "std"])
       .reset_index())
agg.columns = [
    f"{m}_{stat}" if stat in ("mean", "std") else m
    for m, stat in agg.columns
]

# Build readable "mean ± std" strings
summary_rows = []
for _, row in agg.iterrows():
    entry = {
        "Model":    row["Model"],
        "Task":     row["Task_Desc"],
        "Strategy": row["Strategy"],
    }
    for m in metric_cols:
        mu = row.get(f"{m}_mean", float("nan"))
        sd = row.get(f"{m}_std",  float("nan"))
        if pd.notna(mu):
            entry[m] = f"{mu:.4f} ± {sd:.4f}" if pd.notna(sd) else f"{mu:.4f}"
        else:
            entry[m] = ""
    summary_rows.append(entry)

summary_df = pd.DataFrame(summary_rows)

# ── Optionally merge tabular baselines ────────────────────────────────────────
if args.include_baseline:
    baseline_path = OUT_DIR / "baseline" / "baseline_model_comparison.csv"
    if baseline_path.exists():
        bl_df = pd.read_csv(baseline_path)
        # Standardise columns to match
        bl_df["Strategy"] = "tabular"
        summary_df = pd.concat([bl_df, summary_df], ignore_index=True, sort=False)
        print(f"Merged tabular baselines from {baseline_path}")
    else:
        print(f"  [WARNING] Baseline file not found: {baseline_path}")

summary_path = OUT_DIR / "encoder_summary.csv"
summary_df.to_csv(summary_path, index=False)
print(f"Saved summary          → {summary_path}")

# ── Print to console ──────────────────────────────────────────────────────────
print(f"\n{'='*80}")
print("  ENCODER SUMMARY  (mean ± std across seeds 0, 1, 2)")
print(f"{'='*80}")

# Print task by task for readability
for task_desc in summary_df["Task"].unique():
    print(f"\n  [{task_desc}]")
    sub = summary_df[summary_df["Task"] == task_desc]
    # Pick key metric columns that exist
    key_metrics = [c for c in ["BalancedAcc", "AUC_ROC", "AUC_ROC_OvR", "F1", "MacroF1", "AUC_PR"] if c in sub.columns]
    print(sub[["Model", "Strategy"] + key_metrics].to_string(index=False))
