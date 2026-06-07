"""
09_collect_val_test_metrics.py — recompute val + test metrics from checkpoints
==============================================================================
Fills the val AUC / F1 gaps in the cross-model table for the supervised MRI
models that only ever logged val *balanced accuracy* (BrainMVP, AG-MS3D,
Spasov-CNN / "3D-CNN"). For each existing run it loads the saved `best_model.pt`
and recomputes val + test metrics via each trainer's `--val_test` mode, which
patches a `val_metrics` block into the run's metrics.json **without retraining**
(so the published test numbers are preserved — `--val_test` verifies the
recomputed test bACC matches the stored one before overwriting it).

`05_aggregate_mri_results.py:read_run` then reads `val_metrics` (or train_log at
best_epoch), so re-running `06` → `06c` surfaces the new val AUC/F1 columns.

This is designed to run on **Colab A100** (the HPC checkpoints were deleted to
free space; they survive on `D:` and get staged to Drive) but works anywhere the
checkpoints + input volumes + splits are reachable. It is GPU-optional — eval of
the small val/test sets is fine on CPU, just slower for the 3D models.

Per model it discovers runs with the same globs as the aggregator
(`05.MODEL_TREES`), reads each metrics.json `config`, and invokes the matching
trainer with `--val_test --out_dir <the run dir>` plus the architecture args
needed to rebuild a model that matches the checkpoint.

Usage (Colab, after staging Drive paths)::

    python mri_pipeline/09_collect_val_test_metrics.py \
        --derivs-root   /content/drive/MyDrive/ADNI/derivatives \
        --brainmvp-inputs /content/drive/MyDrive/ADNI/brainmvp_inputs \
        --cnn-inputs      /content/drive/MyDrive/ADNI/cnn_inputs \
        --data-dir        /content/drive/MyDrive/ADNI/no_cdr_stratified_post_exclusion/tabular/baseline \
        --matched-labels  /content/drive/MyDrive/ADNI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv \
        --brainmvp-ckpt   /content/drive/MyDrive/ADNI/ViT_pretrained/BrainMVP_uniformer.pt \
        --models brainmvp cnn3d agms3d
    # add --dry-run to print the commands without running
    # add --only-missing to skip runs whose metrics.json already has val_metrics
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent

# Reuse the aggregator's run-discovery globs + config reader (single source of
# truth for where each model's runs live on disk).
_spec = importlib.util.spec_from_file_location(
    "_05_aggregate", str(_HERE / "05_aggregate_mri_results.py"))
_05 = importlib.util.module_from_spec(_spec)
sys.modules["_05_aggregate"] = _05
_spec.loader.exec_module(_05)


# Trainer scripts, keyed by the --models token.
_VIT = _HERE / "04_supervised_finetuning_ViT.py"
TRAINERS = {
    "brainmvp":    _HERE / "brain_mvp" / "04_supervised_finetuning_BrainMVP.py",
    "agms3d":      _HERE / "3d_cnn_vit" / "train_agms3d.py",
    "cnn3d":       _HERE / "3d_conv_net" / "train_3dcnn.py",
    "vit_mae":     _VIT,
    "vit_scratch": _VIT,
}

# Which MODEL_TREES labels (from 05) map to each trainer token.
MODELS_BY_TOKEN = {
    "brainmvp":    ["BrainMVP"],
    "agms3d":      ["AG-MS3D-sep", "AG-MS3D-vanilla"],
    "cnn3d":       ["Spasov-CNN"],
    "vit_mae":     ["ViT-MAE75"],
    "vit_scratch": ["ViT-scratch"],
}


def _long_args(cfg):
    """Reconstruct the split-selection flags BrainMVP was trained with."""
    lm = cfg.get("long_mode")
    if lm == "all":
        return ["--long", "all"]
    if lm == "cutoff" and cfg.get("max_months"):
        return ["--long", str(int(cfg["max_months"]) // 12)]
    # default: baseline session
    return ["--session", "bl"]


def _build_cmd(token, run_dir, cfg, args):
    """Build the trainer CLI for a --val_test recompute of one run."""
    py = sys.executable
    script = str(TRAINERS[token])
    common = ["--task", cfg["task"], "--seed", str(cfg["seed"]),
              "--out_dir", str(run_dir), "--val_test",
              "--num_workers", str(args.num_workers),
              "--matched_labels_csv", args.matched_labels,
              "--data_dir", args.data_dir]

    if token == "brainmvp":
        if not args.brainmvp_ckpt:
            return None, "needs --brainmvp-ckpt"
        cmd = [py, script, *common,
               "--strategy", cfg.get("strategy", "full_ft"),
               "--augment", cfg.get("augment", "stochastic"),
               "--pretrained_ckpt", args.brainmvp_ckpt,
               "--brainmvp_inputs_dir", args.brainmvp_inputs,
               *_long_args(cfg)]
        return cmd, None

    if token == "cnn3d":
        # model_kind is "vanilla"/"separable"; --model takes that verbatim.
        model_kind = cfg.get("model_kind") or "vanilla"
        cmd = [py, script, *common,
               "--model", model_kind,
               "--cnn_inputs_dir", args.cnn_inputs]
        return cmd, None

    if token == "agms3d":
        # Rescue runs store backbone/head/base_filters; the legacy "AGMS3DCNN"
        # tree predates them -> infer (separable backbone, large head, base 32).
        backbone = cfg.get("backbone") or "separable"
        head = cfg.get("head") or "large"
        base = cfg.get("base_filters") or 32
        cmd = [py, script, *common,
               "--backbone", backbone, "--head", head,
               "--base_filters", str(base),
               "--cnn_inputs_dir", args.cnn_inputs]
        return cmd, None

    if token in ("vit_mae", "vit_scratch"):
        # ViT rebuilds the model (loading the MAE pretrained ckpt for full_ft/
        # frozen) BEFORE overwriting with best_model.pt, so non-scratch needs
        # --vit-ckpt. vit_size/strategy/augment must match the run; drop/attn/ls
        # don't change tensor shapes but are passed to stay faithful.
        strategy = cfg.get("strategy", "full_ft")
        cmd = [py, script, *common,
               "--strategy", strategy,
               "--vit_size", cfg.get("vit_size", "base"),
               "--augment", cfg.get("augment", "random"),
               "--vit_inputs_dir", args.vit_inputs,
               *_long_args(cfg)]
        if strategy != "scratch":
            if not args.vit_ckpt:
                return None, "needs --vit-ckpt (MAE pretrained)"
            cmd += ["--pretrained_ckpt", args.vit_ckpt]
        for k_cfg, flag in [("drop_path_rate", "--drop_path_rate"),
                            ("attn_dropout", "--attn_dropout"),
                            ("label_smoothing", "--label_smoothing")]:
            if cfg.get(k_cfg) is not None:
                cmd += [flag, str(cfg[k_cfg])]
        return cmd, None

    return None, f"unknown token {token}"


def _already_has_val_metrics(metrics_path):
    try:
        with open(metrics_path) as f:
            return bool((json.load(f) or {}).get("val_metrics"))
    except Exception:
        return False


def _discover(token, derivs_root):
    """Yield (run_dir, metrics_path, cfg) for every run of the model(s) mapped
    to `token`, using 05's globs (with the double-nested rsync fallback)."""
    out = []
    for label in MODELS_BY_TOKEN[token]:
        for model, rel in _05.MODEL_TREES:
            if model != label:
                continue
            tree = rel.split("/")[0]
            files = sorted(glob.glob(os.path.join(derivs_root, rel)))
            if not files:
                files = sorted(glob.glob(os.path.join(derivs_root, tree, rel)))
            for f in files:
                try:
                    cfg = json.load(open(f)).get("config", {})
                except Exception as exc:
                    print(f"  [WARN] unreadable {f}: {exc}")
                    continue
                if "task" in cfg and cfg.get("seed") is not None:
                    out.append((os.path.dirname(f), f, cfg))
    return out


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--derivs-root", required=True,
                   help="Root holding the model run trees (brainmvp_debug/, "
                        "cnn3d_outputs/, agms3d_outputs/, ...).")
    p.add_argument("--brainmvp-inputs", default="")
    p.add_argument("--cnn-inputs", default="")
    p.add_argument("--vit-inputs", default="",
                   help="vit_inputs dir (128^3 volumes) for ViT recompute.")
    p.add_argument("--vit-ckpt", default="",
                   help="ViT-B MAE pretrained .pth.tar (needed to build the model "
                        "before loading best_model.pt; not needed for scratch).")
    p.add_argument("--data-dir", required=True,
                   help="Post-exclusion splits root (seed_*/{train,val,test}.csv).")
    p.add_argument("--matched-labels", required=True)
    p.add_argument("--brainmvp-ckpt", default="",
                   help="BrainMVP_uniformer.pt (needed to build the model "
                        "before loading best_model.pt).")
    p.add_argument("--models", nargs="+", default=["brainmvp", "cnn3d", "agms3d"],
                   choices=list(TRAINERS.keys()))
    p.add_argument("--num-workers", type=int, default=0)
    p.add_argument("--only-missing", action="store_true",
                   help="Skip runs whose metrics.json already has a val_metrics block.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the commands without running them.")
    return p.parse_args()


def main():
    args = parse_args()
    print("=" * 78)
    print("  09_collect_val_test_metrics — recompute val+test from checkpoints")
    print(f"  derivs-root: {args.derivs_root}")
    print(f"  models     : {args.models}")
    print("=" * 78)

    n_ok = n_skip = n_fail = n_nockpt = 0
    for token in args.models:
        runs = _discover(token, args.derivs_root)
        print(f"\n[{token}] {len(runs)} run(s) discovered.")
        for run_dir, metrics_path, cfg in runs:
            tag = f"{token}:{cfg.get('task')}/seed_{cfg.get('seed')}:{Path(run_dir).name}"
            if args.only_missing and _already_has_val_metrics(metrics_path):
                print(f"  [skip] {tag} (already has val_metrics)"); n_skip += 1
                continue
            if not (Path(run_dir) / "best_model.pt").exists():
                print(f"  [no-ckpt] {tag} (no best_model.pt)"); n_nockpt += 1
                continue
            cmd, err = _build_cmd(token, run_dir, cfg, args)
            if cmd is None:
                print(f"  [skip] {tag}: {err}"); n_skip += 1
                continue
            if args.dry_run:
                print("  [dry-run] " + " ".join(cmd)); n_ok += 1
                continue
            print(f"  [run ] {tag}")
            rc = subprocess.run(cmd, cwd=str(_REPO)).returncode
            if rc == 0:
                n_ok += 1
            else:
                print(f"  [FAIL] {tag} (rc={rc})"); n_fail += 1

    print("\n" + "=" * 78)
    print(f"  done: {n_ok} ok, {n_skip} skipped, {n_nockpt} missing-ckpt, "
          f"{n_fail} failed.")
    print("  Re-run 06_render_cross_model_table.py then 06c to surface val AUC/F1.")
    print("=" * 78)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
