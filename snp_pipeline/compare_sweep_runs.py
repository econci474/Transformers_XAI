"""
compare_sweep_runs.py
======================
After all SLURM array tasks from 09_bmfm_gwas_regression_finetuning_sweep.sh
have completed, this script reads each run_config.json and finds:

  1. The best val_loss checkpoint across all runs
  2. A comparison table of all configs

Usage (on CSD3 or locally after downloading the sweep output):
  python snp_pipeline/compare_sweep_runs.py
  python snp_pipeline/compare_sweep_runs.py --sweep-dir ~/ADNI_SNP/bmfm_gwas_sweep
"""
import argparse
import json
import pathlib

_p = argparse.ArgumentParser()
_p.add_argument("--sweep-dir", default=None,
                help="Root directory containing sweep run subdirectories. "
                     "Default: ~/ADNI_SNP/bmfm_gwas_sweep")
args = _p.parse_args()

SWEEP_DIR = pathlib.Path(args.sweep_dir) if args.sweep_dir \
            else pathlib.Path.home() / "ADNI_SNP" / "bmfm_gwas_sweep"

print(f"\nScanning sweep directory: {SWEEP_DIR}")

configs = []
for cfg_path in sorted(SWEEP_DIR.rglob("run_config.json")):
    try:
        cfg = json.loads(cfg_path.read_text())
        # Extract val_loss from best checkpoint filename
        ckpt = cfg.get("best_checkpoint", "")
        val_loss = None
        if ckpt:
            for part in pathlib.Path(ckpt).name.split("-"):
                if part.startswith("val_loss="):
                    try:
                        val_loss = float(part.split("=")[1].replace(".ckpt", ""))
                    except ValueError:
                        pass
        cfg["_val_loss"] = val_loss
        configs.append(cfg)
    except Exception as e:
        print(f"  [WARN] Could not read {cfg_path}: {e}")

if not configs:
    print("No run_config.json files found. Ensure sweep has completed.")
    raise SystemExit(1)

# Sort by val_loss ascending
configs.sort(key=lambda c: c.get("_val_loss") or float("inf"))

# Print table
print(f"\n{'='*80}")
print(f"  {'Epochs':>7}  {'Loss':>5}  {'LR':>8}  {'val_loss':>10}  Best checkpoint")
print(f"  {'-'*7}  {'-'*5}  {'-'*8}  {'-'*10}  {'-'*40}")
for cfg in configs:
    vl = cfg.get("_val_loss")
    vl_str = f"{vl:.4f}" if vl is not None else "N/A"
    ckpt_name = pathlib.Path(cfg.get("best_checkpoint", "N/A")).name
    marker = " ← BEST" if cfg == configs[0] else ""
    print(f"  {cfg.get('max_epochs', '?'):>7}  "
          f"{cfg.get('loss', '?'):>5}  "
          f"{cfg.get('learning_rate', '?'):>8}  "
          f"{vl_str:>10}  {ckpt_name}{marker}")

print(f"{'='*80}")

best = configs[0]
print(f"\n  Best run summary:")
print(f"    max_epochs       : {best.get('max_epochs')}")
print(f"    loss             : {best.get('loss')}")
print(f"    learning_rate    : {best.get('learning_rate')}")
print(f"    val_loss         : {best.get('_val_loss')}")
print(f"    best_checkpoint  : {best.get('best_checkpoint')}")
print(f"    output_directory : {best.get('output_directory')}")

# Write summary CSV
import csv
out_csv = SWEEP_DIR / "sweep_comparison.csv"
fieldnames = ["max_epochs", "loss", "learning_rate", "batch_size",
              "accumulate_grad_batches", "_val_loss", "best_checkpoint",
              "exit_code", "finished_at"]
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    w.writeheader()
    w.writerows(configs)
print(f"\n  Comparison table saved → {out_csv}")
