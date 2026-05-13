"""
16_aggregate_classification_results.py
=======================================
Walk a `classification_results/` tree of metrics.json files (written by
15_train_classification_head.py) and compile a comparison table:

  rows = (upstream_tag, aggregation)
  cols = test_auc, test_balanced_accuracy, test_f1, test_sensitivity,
         test_specificity   — each as "mean ± std" across seeds.

Also writes a flat long-form CSV with one row per (upstream, aggregation, seed)
for downstream plotting.

Usage
-----
  python 16_aggregate_classification_results.py
  python 16_aggregate_classification_results.py --root classification_results/
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = ["auc", "balanced_accuracy", "f1", "sensitivity", "specificity"]


def collect(root: Path) -> pd.DataFrame:
    rows = []
    for mfile in root.rglob("metrics.json"):
        with open(mfile) as f:
            m = json.load(f)
        row = {
            "label_mode": m.get("label_mode", "bl_multi"),
            "upstream_tag": m["upstream_tag"],
            "pool": m.get("pool", "cls"),  # legacy runs default to cls
            "aggregation": m["aggregation"],
            "seed": m["seed"],
        }
        for k in METRICS:
            row[f"val_{k}"] = m["val"].get(k)
            row[f"test_{k}"] = m["test"].get(k)
        row["best_val_auc"] = m.get("fit_info", {}).get("best_val_auc")
        row["n_epochs"] = m.get("fit_info", {}).get("n_epochs")
        rows.append(row)
    if not rows:
        raise SystemExit(f"No metrics.json found under {root}")
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(["label_mode", "upstream_tag", "pool", "aggregation"])
    out_rows = []
    for (label_mode, upstream, pool, agg), grp in g:
        row = {"label_mode": label_mode, "upstream_tag": upstream,
               "pool": pool, "aggregation": agg, "n_seeds": len(grp)}
        for k in METRICS:
            col = f"test_{k}"
            vals = grp[col].dropna().to_numpy()
            if len(vals) == 0:
                row[col] = ""
                continue
            mean = float(np.mean(vals))
            std = float(np.std(vals, ddof=0))
            row[col] = f"{mean:.3f} ± {std:.3f}"
            row[f"{col}_mean"] = mean
            row[f"{col}_std"] = std
        out_rows.append(row)
    return pd.DataFrame(out_rows).sort_values(
        ["label_mode", "upstream_tag", "pool", "aggregation"]
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=Path("classification_results"))
    p.add_argument("--long-csv", type=Path,
                   default=Path("classification_results/comparison_long.csv"))
    p.add_argument("--summary-csv", type=Path,
                   default=Path("classification_results/comparison_summary.csv"))
    args = p.parse_args()

    df_long = collect(args.root)
    df_long.to_csv(args.long_csv, index=False)
    print(f"Wrote {args.long_csv}  ({len(df_long)} rows)")

    df_sum = summarise(df_long)
    df_sum.to_csv(args.summary_csv, index=False)
    print(f"Wrote {args.summary_csv}  ({len(df_sum)} rows)")

    # Pretty-print the headline table.
    cols = ["label_mode", "upstream_tag", "pool", "aggregation", "n_seeds",
            "test_auc", "test_balanced_accuracy", "test_f1",
            "test_sensitivity", "test_specificity"]
    with pd.option_context("display.max_columns", None,
                           "display.width", 220,
                           "display.expand_frame_repr", False):
        print("\n" + "=" * 100)
        print(f"  Cross-run summary  →  {args.summary_csv}")
        print("=" * 100)
        print(df_sum[cols].to_string(index=False))


if __name__ == "__main__":
    main()
