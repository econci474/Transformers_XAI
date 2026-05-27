"""
06_fetch_wandb_tables.py -- pull per-project HP + val_metrics tables from W&B
============================================================================
For each named W&B project, dump every run's config (hyperparameters) and
summary (final metric values) into a single CSV. Runs in `state == "running"`
are skipped by default so the table reflects only finished sweeps.

The five projects this script targets for the thesis MRI comparison:
    vit_mae_frozen_cached
    brainmvp_frozen_cached
    braindino_frozen_cached
    brainmvp_debug
    vit_debugging

Output:
    <out-dir>/wandb_<project>.csv     one CSV per project

Caveat -- val_* metric availability:
    The cached-head trainer (04_head_finetune_from_embeddings.py) only logs
    val_bacc and val_loss to W&B. Other val_* metrics the user asked for
    (val_auc, val_f1, val_sens, val_spec, val_prec, val_recall, val_TP,
    val_FP, val_TN, val_FN) are computed on the TEST set in metrics.json
    but never logged to W&B. Those columns will be NaN in the per-project
    CSVs from the *_frozen_cached projects. The brainmvp_debug and
    vit_debugging projects come from earlier non-cached trainers and may
    log a richer set -- the script dumps whatever it finds.

Usage (Windows laptop):
    python mri_pipeline/cached_head_sweep/06_fetch_wandb_tables.py
    python mri_pipeline/cached_head_sweep/06_fetch_wandb_tables.py \
        --entity <wandb-entity> \
        --projects vit_mae_frozen_cached braindino_frozen_cached \
        --out-dir mri_pipeline/cached_head_sweep \
        --include-running

Requirements:
    pip install wandb pandas       # wandb auth via `wandb login` (one-time)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

DEFAULT_PROJECTS = [
    "vit_mae_frozen_cached",
    "brainmvp_frozen_cached",
    "braindino_frozen_cached",
    "brainmvp_debug",
    "vit_debugging",
]

# The val_* metrics the user explicitly wants in each project's CSV. Whatever
# is missing in W&B will land as NaN; the script also writes everything else
# the project's runs logged so nothing useful is dropped silently.
REQUESTED_VAL_COLS = [
    "val_bacc",
    "val_auc",
    "val_f1",
    "val_sens",
    "val_spec",
    "val_prec",
    "val_recall",
    "val_TP",
    "val_FP",
    "val_TN",
    "val_FN",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--entity", default=None,
                   help="W&B entity (user or team). Defaults to the logged-in "
                        "user's default entity.")
    p.add_argument("--projects", nargs="+", default=DEFAULT_PROJECTS,
                   help=f"Project names to pull. Default: {DEFAULT_PROJECTS}")
    p.add_argument("--out-dir", default=os.path.dirname(os.path.abspath(__file__)),
                   help="Directory to write per-project CSVs.")
    p.add_argument("--include-running", action="store_true",
                   help="Also include runs whose state is 'running'.")
    p.add_argument("--per-page", type=int, default=500,
                   help="W&B Api pagination chunk size.")
    return p.parse_args()


def _flatten_config(cfg: dict) -> dict:
    """W&B run.config is a Config object; coerce to a flat dict whose values
    are JSON-serialisable scalars or short strings. Nested dicts are dotted."""
    if not cfg:
        return {}
    out = {}
    for k, v in dict(cfg).items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                out[f"config/{k}.{kk}"] = vv
        else:
            out[f"config/{k}"] = v
    return out


def _flatten_summary(summary) -> dict:
    """run.summary is a SummarySubDict; pull its keys with a 'summary/'
    prefix so they don't collide with config keys of the same name."""
    if summary is None:
        return {}
    out = {}
    # SummarySubDict supports .items(); some versions don't, so be defensive.
    try:
        items = summary.items()
    except AttributeError:
        try:
            items = dict(summary).items()
        except Exception:
            return {}
    for k, v in items:
        # Strip internal/_wandb keys that aren't useful in a comparison table
        if str(k).startswith("_"):
            continue
        # If the summary value is a dict (e.g. a confusion matrix logged as
        # an artifact), serialise to a string so pandas can write it.
        if isinstance(v, dict):
            out[f"summary/{k}"] = str(v)
        else:
            out[f"summary/{k}"] = v
    return out


def _project_to_rows(api, entity: str, project: str, include_running: bool,
                     per_page: int):
    """Iterate all runs in <entity>/<project>, return list of flat row dicts."""
    project_path = f"{entity}/{project}" if entity else project
    try:
        runs = api.runs(project_path, per_page=per_page)
    except Exception as exc:
        print(f"  [SKIP] {project} -- {exc}")
        return []

    rows = []
    n_running = 0
    for run in runs:
        if run.state == "running" and not include_running:
            n_running += 1
            continue

        row = {
            "run_id":      run.id,
            "run_name":    run.name,
            "state":       run.state,
            "created_at":  str(run.created_at),
            "tags":        ", ".join(run.tags or []),
            "url":         run.url,
            "runtime_sec": run.summary.get("_runtime") if run.summary else None,
        }
        row.update(_flatten_config(run.config))
        row.update(_flatten_summary(run.summary))
        rows.append(row)

    print(f"  [ok]  {project:35s} -- {len(rows)} runs"
          + (f"  (+{n_running} still running, skipped)" if n_running else ""))
    return rows


def _ensure_requested_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Make sure each REQUESTED_VAL_COL appears in the output, with NaN if
    not logged by the trainer. The W&B keys live under 'summary/...'."""
    for col in REQUESTED_VAL_COLS:
        key = f"summary/{col}"
        if key not in df.columns:
            df[key] = pd.NA
    return df


def _reorder_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Put identity + requested val cols + config first, summary metrics last,
    so the user can spot the requested metrics without horizontal-scrolling."""
    id_cols = ["run_id", "run_name", "state", "created_at", "tags",
               "url", "runtime_sec"]
    requested = [f"summary/{c}" for c in REQUESTED_VAL_COLS]
    config_cols = sorted([c for c in df.columns if c.startswith("config/")])
    summary_cols = sorted([c for c in df.columns
                           if c.startswith("summary/") and c not in requested])
    other = [c for c in df.columns
             if c not in id_cols + requested + config_cols + summary_cols]
    return df[id_cols + requested + config_cols + summary_cols + other]


def _emit_missing_summary(df: pd.DataFrame, project: str):
    """Print which of the requested val_* columns came back empty."""
    missing = []
    populated = []
    for col in REQUESTED_VAL_COLS:
        key = f"summary/{col}"
        if key in df.columns and df[key].notna().any():
            populated.append(col)
        else:
            missing.append(col)
    if missing:
        print(f"        not in W&B for {project}: {', '.join(missing)}")
    if populated:
        print(f"        populated:              {', '.join(populated)}")


def main():
    args = parse_args()

    try:
        import wandb
    except ImportError:
        sys.exit("`wandb` is not installed. Run `pip install wandb` first.")

    api = wandb.Api()
    entity = args.entity or api.default_entity
    if not entity:
        sys.exit("Could not resolve W&B entity. Either pass --entity or run "
                 "`wandb login` first to set a default.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("  06_fetch_wandb_tables -- per-project HP + summary CSVs")
    print(f"  Entity   : {entity}")
    print(f"  Out      : {out_dir}")
    print(f"  Projects : {len(args.projects)}")
    print("=" * 78)

    for project in args.projects:
        rows = _project_to_rows(
            api, entity, project, args.include_running, args.per_page)
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df = _ensure_requested_cols(df)
        df = _reorder_cols(df)
        out_path = out_dir / f"wandb_{project}.csv"
        df.to_csv(out_path, index=False)
        print(f"        CSV: {out_path}  shape={df.shape}")
        _emit_missing_summary(df, project)

    print("=" * 78)
    print("  Done.")


if __name__ == "__main__":
    main()
