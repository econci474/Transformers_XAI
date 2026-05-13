#!/usr/bin/env python
"""
Extract BMFM training metrics from TensorBoard event files across all
regression runs under D:\\ADNI_SNP_Omni2.5M_20140220\\bmfm_outputs_regression\\.

Two modes:

1. Multi-run (default):
       python 10_unpack_bmfm_results.py
       python 10_unpack_bmfm_results.py --root <path/to/bmfm_outputs_regression>

   Walks `<root>/*/lightning_logs/version_*/`, writes per-run
   `bmfm_training_metrics_{full,summary}.csv` next to each `lightning_logs/`,
   and emits two top-level summaries:
       - <root>/all_runs_comparison.csv  (non-empty versions; hparams + final
         metrics + target z-score stats from the matching train.csv +
         r2_from_val_loss)
       - <root>/failed_runs.csv          (versions that only logged
         hp_metric=-1.0, with job_id / host / mtime for SLURM-log triage)

2. Single-run back-compat:
       python 10_unpack_bmfm_results.py --logdir <path/to/lightning_logs>

   Behaves as the original script: writes per-run CSVs next to the given
   lightning_logs dir and prints the per-split breakdown for its latest
   version. No top-level summaries are written.
"""
import argparse
import datetime as _dt
import re
from pathlib import Path

import pandas as pd
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

# ── Defaults ─────────────────────────────────────────────────────────────────
DEFAULT_ROOT = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\bmfm_outputs_regression")
INPUTS_ROOT = Path(r"D:\ADNI_SNP_Omni2.5M_20140220\bmfm_inputs")

# Map output-dataset directory name → input train.csv (relative to INPUTS_ROOT)
# Used to look up the target z-score distribution so MSE / R² are interpretable.
DATASET_TO_INPUT_TRAIN = {
    "bmfm_gwas_signed_regression_without_ukb": (
        Path("bmfm_gwas_signed_regression_without_ukb") / "train.csv",
    ),
    "bmfm_gwas_signed_regression_without_ukb_augmented_lr1e6": (
        Path("bmfm_gwas_signed_regression_without_ukb_augmented") / "train.csv",
    ),
    "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom": (
        Path("bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom")
        / "forward_and_reverse" / "train.csv",
        Path("bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom")
        / "forward" / "train.csv",
    ),
    # CSD3 combos runs (chromosome-context 8 kb windows, n≈13.8k rows, σ≈31.87)
    "bmfm_gwas_combos_lr1e6_2ep": (
        Path("bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos")
        / "forward_and_reverse" / "train.csv",
    ),
    "bmfm_gwas_combos_ref_lr1e6_2ep": (
        Path("bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos")
        / "forward_and_reverse" / "train.csv",
    ),
}

# Metric tags surfaced in the comparison table.
COMPARISON_METRICS = ["epoch", "train/loss_epoch", "validation/loss", "test/loss"]

# Hparams keys to scrape from hparams.yaml.
# Only `learning_rate` (from trainer_config) is reliably present and varied.
# `max_length` in hparams.yaml refers to the model's text-generation output
# length (default 20), NOT the dataloader's sequence length — Hydra overrides
# for the latter (and for batch_size / max_epochs) do not appear in the
# Lightning hparams snapshot, so we don't try to scrape them.
HPARAM_KEYS = ["learning_rate"]


# ── hparams.yaml regex parser ────────────────────────────────────────────────
# Lightning writes !!python/object tags, so yaml.safe_load rejects the file and
# yaml.full_load would execute arbitrary code. Regex is sufficient for the
# scalar keys we care about and is safe.
def parse_hparams_text(text: str) -> dict:
    out = {}
    for key in HPARAM_KEYS:
        m = re.search(rf"^\s*{re.escape(key)}:\s*(\S+)\s*$", text, re.MULTILINE)
        if not m:
            out[key] = None
            continue
        val = m.group(1)
        # Strip wrapping quotes if any.
        if (val.startswith('"') and val.endswith('"')) or (
            val.startswith("'") and val.endswith("'")
        ):
            val = val[1:-1]
        # Try numeric coercion; fall back to string.
        try:
            out[key] = float(val) if "." in val or "e" in val.lower() else int(val)
        except ValueError:
            out[key] = val
    return out


def read_hparams(version_dir: Path) -> dict:
    f = version_dir / "hparams.yaml"
    if not f.exists():
        return {k: None for k in HPARAM_KEYS}
    return parse_hparams_text(f.read_text(encoding="utf-8", errors="replace"))


# ── Event file metadata (used for failure triage) ────────────────────────────
EVENT_FILENAME_RE = re.compile(r"^events\.out\.tfevents\.\d+\.(.+)$")


def event_file_info(version_dir: Path) -> dict:
    events = sorted(version_dir.glob("events.out.tfevents.*"))
    if not events:
        return {"event_file_bytes": None, "event_file_host": None, "event_file_mtime": None}
    # Use the largest event file as the canonical one (resumes / multi-write).
    canonical = max(events, key=lambda p: p.stat().st_size)
    stat = canonical.stat()
    host = None
    m = EVENT_FILENAME_RE.match(canonical.name)
    if m:
        host = m.group(1)
    return {
        "event_file_bytes": stat.st_size,
        "event_file_host": host,
        "event_file_mtime": _dt.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


# ── Target z-score stats per dataset (cached) ────────────────────────────────
_target_cache: dict[str, dict] = {}


def target_z_stats(dataset_name: str) -> dict:
    if dataset_name in _target_cache:
        return _target_cache[dataset_name]
    candidates = DATASET_TO_INPUT_TRAIN.get(dataset_name)
    stats = {"target_train_n": None, "target_train_mean": None, "target_train_std": None}
    if not candidates:
        _target_cache[dataset_name] = stats
        return stats
    for rel in candidates:
        p = INPUTS_ROOT / rel
        if p.exists():
            try:
                z = pd.read_csv(p, usecols=["z_score"])["z_score"]
                stats = {
                    "target_train_n": int(len(z)),
                    "target_train_mean": float(z.mean()),
                    "target_train_std": float(z.std(ddof=0)),
                }
            except Exception as e:
                print(f"  [warn] failed to read {p}: {e}")
            break
    else:
        print(f"  [warn] no input train.csv found for dataset {dataset_name}")
    _target_cache[dataset_name] = stats
    return stats


# ── Core extraction (one lightning_logs dir → per-version DataFrames) ────────
def extract_lightning_logs(logdir: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Walk every version_* under `logdir`, collect scalar events.

    Returns (df_full, df_summary, version_meta) where version_meta is a list of
    one dict per version with hparams + event-file metadata + scalar tag set.
    """
    version_dirs = sorted([d for d in logdir.glob("version_*") if d.is_dir()])
    if not version_dirs:
        version_dirs = [logdir]

    all_records: list[dict] = []
    version_meta: list[dict] = []

    for vdir in version_dirs:
        version_name = vdir.name if vdir != logdir else "root"
        ea = EventAccumulator(str(vdir))
        ea.Reload()
        tags = ea.Tags().get("scalars", [])

        meta = {
            "version": version_name,
            "scalar_tags": list(tags),
            **read_hparams(vdir),
            **event_file_info(vdir),
        }
        version_meta.append(meta)

        if not tags:
            continue
        for tag in sorted(tags):
            for e in ea.Scalars(tag):
                all_records.append({
                    "version": version_name,
                    "metric": tag,
                    "step": e.step,
                    "value": e.value,
                    "wall_time": e.wall_time,
                })

    if all_records:
        df_full = pd.DataFrame(all_records)
        df_summary = (
            df_full.sort_values(["version", "metric", "step"])
                   .groupby(["version", "metric"])
                   .last()
                   .reset_index()[["version", "metric", "step", "value"]]
        )
    else:
        df_full = pd.DataFrame(columns=["version", "metric", "step", "value", "wall_time"])
        df_summary = pd.DataFrame(columns=["version", "metric", "step", "value"])

    return df_full, df_summary, version_meta


def write_per_run_csvs(logdir: Path, df_full: pd.DataFrame, df_summary: pd.DataFrame) -> None:
    out_dir = logdir.parent
    df_full.to_csv(out_dir / "bmfm_training_metrics_full.csv", index=False)
    df_summary.to_csv(out_dir / "bmfm_training_metrics_summary.csv", index=False)


# ── Classify each version as empty (only hp_metric=-1.0) or real ─────────────
def is_empty_version(meta: dict, df_summary: pd.DataFrame) -> bool:
    tags = meta["scalar_tags"]
    if tags == ["hp_metric"]:
        row = df_summary[(df_summary["version"] == meta["version"]) & (df_summary["metric"] == "hp_metric")]
        if not row.empty and float(row["value"].iloc[0]) == -1.0:
            return True
    return False


# ── Job ID extraction from version name ──────────────────────────────────────
def job_id_from_version(version: str) -> str | None:
    if not version.startswith("version_"):
        return None
    rest = version[len("version_"):]
    # SLURM job IDs are large integers (≥ 10⁷); local-Lightning versions are
    # small (0, 1, 2, …). Treat only the former as job IDs.
    if rest.isdigit() and int(rest) >= 10_000_000:
        return rest
    return None


# ── Build the cross-run comparison + failed-runs tables ──────────────────────
def summary_value(df_summary: pd.DataFrame, version: str, metric: str):
    row = df_summary[(df_summary["version"] == version) & (df_summary["metric"] == metric)]
    if row.empty:
        return None
    return float(row["value"].iloc[0])


def build_comparison_and_failed(
    per_dataset: list[tuple[str, pd.DataFrame, list[dict]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    comp_rows: list[dict] = []
    fail_rows: list[dict] = []
    for dataset_name, df_summary, version_meta in per_dataset:
        z = target_z_stats(dataset_name)
        for meta in version_meta:
            v = meta["version"]
            if is_empty_version(meta, df_summary):
                fail_rows.append({
                    "dataset": dataset_name,
                    "version": v,
                    "job_id": job_id_from_version(v),
                    "event_file_bytes": meta["event_file_bytes"],
                    "event_file_host": meta["event_file_host"],
                    "event_file_mtime": meta["event_file_mtime"],
                    "target_train_std": z["target_train_std"],
                })
                continue
            # Non-empty: collect comparison row.
            val_loss = summary_value(df_summary, v, "validation/loss")
            std = z["target_train_std"]
            r2 = None
            if val_loss is not None and std is not None and std > 0:
                r2 = 1.0 - val_loss / (std ** 2)
            comp_rows.append({
                "dataset": dataset_name,
                "version": v,
                "learning_rate": meta["learning_rate"],
                "target_train_n": z["target_train_n"],
                "target_train_mean": z["target_train_mean"],
                "target_train_std": z["target_train_std"],
                "last_epoch": summary_value(df_summary, v, "epoch"),
                "train_loss_epoch": summary_value(df_summary, v, "train/loss_epoch"),
                "validation_loss": val_loss,
                "test_loss": summary_value(df_summary, v, "test/loss"),
                "r2_from_val_loss": r2,
            })
    return pd.DataFrame(comp_rows), pd.DataFrame(fail_rows)


# ── Single-run mode: replicate original per-split breakdown print ────────────
def print_per_split_breakdown(df_summary: pd.DataFrame) -> None:
    if df_summary.empty:
        return
    latest_version = df_summary["version"].iloc[-1]
    latest = df_summary[df_summary["version"] == latest_version]
    print("\n" + "=" * 72)
    print(f"  Per-Split Summary  (latest version: {latest_version})")
    print("=" * 72)
    for split in ["train", "validation", "test"]:
        split_metrics = latest[latest["metric"].str.startswith(split + "/")]
        if split_metrics.empty:
            continue
        print(f"\n  {split.upper()}:")
        for _, row in split_metrics.iterrows():
            print(f"    {row['metric']:40s}  {row['value']:.6f}  (step {int(row['step'])})")
    print()


# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(description="Unpack BMFM TensorBoard logs across all regression runs.")
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help="Root containing one subdirectory per dataset (each with lightning_logs/).")
    ap.add_argument("--logdir", type=Path, default=None,
                    help="Back-compat: process a single lightning_logs/ directory and exit.")
    args = ap.parse_args()

    # --- Single-run back-compat mode ---
    if args.logdir is not None:
        if not args.logdir.exists():
            raise FileNotFoundError(f"--logdir not found: {args.logdir}")
        df_full, df_summary, version_meta = extract_lightning_logs(args.logdir)
        if df_full.empty:
            print("No scalar metrics found.")
            raise SystemExit(1)
        write_per_run_csvs(args.logdir, df_full, df_summary)
        print(f"Wrote {args.logdir.parent / 'bmfm_training_metrics_full.csv'}")
        print(f"Wrote {args.logdir.parent / 'bmfm_training_metrics_summary.csv'}")
        print_per_split_breakdown(df_summary)
        return

    # --- Multi-run mode ---
    root: Path = args.root
    if not root.exists():
        raise FileNotFoundError(f"--root not found: {root}")

    dataset_logdirs = sorted(root.glob("*/lightning_logs"))
    if not dataset_logdirs:
        print(f"No */lightning_logs found under {root}")
        raise SystemExit(1)

    print(f"Found {len(dataset_logdirs)} dataset run(s) under {root}\n")

    per_dataset: list[tuple[str, pd.DataFrame, list[dict]]] = []
    for logdir in dataset_logdirs:
        dataset_name = logdir.parent.name
        print(f"── {dataset_name} ──")
        df_full, df_summary, version_meta = extract_lightning_logs(logdir)
        n_real = sum(not is_empty_version(m, df_summary) for m in version_meta)
        n_empty = len(version_meta) - n_real
        print(f"  versions: {len(version_meta)}  (real={n_real}, empty={n_empty})")
        write_per_run_csvs(logdir, df_full, df_summary)
        per_dataset.append((dataset_name, df_summary, version_meta))

    df_comp, df_fail = build_comparison_and_failed(per_dataset)

    comp_path = root / "all_runs_comparison.csv"
    fail_path = root / "failed_runs.csv"
    df_comp.to_csv(comp_path, index=False)
    df_fail.to_csv(fail_path, index=False)

    # --- Console output ---
    print("\n" + "=" * 96)
    print(f"  Cross-run comparison  →  {comp_path}")
    print("=" * 96)
    if df_comp.empty:
        print("  (no non-empty versions found)")
    else:
        with pd.option_context("display.max_columns", None, "display.width", 220, "display.float_format", "{:.4g}".format):
            print(df_comp.to_string(index=False))

    print("\n" + "=" * 96)
    print(f"  Failed runs to debug  →  {fail_path}")
    print("=" * 96)
    if df_fail.empty:
        print("  (none — all versions logged real metrics)")
    else:
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(df_fail.to_string(index=False))
        job_ids = [j for j in df_fail["job_id"].tolist() if isinstance(j, str) and j]
        if job_ids:
            joined = ",".join(job_ids)
            print(
                "\n  → Pull SLURM logs from CSD3:\n"
                f"      mkdir -p debug_logs && \\\n"
                f"      scp 'ec474@login.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_SNP/logs/*{{{joined}}}*' debug_logs/"
            )
    print()


if __name__ == "__main__":
    main()
