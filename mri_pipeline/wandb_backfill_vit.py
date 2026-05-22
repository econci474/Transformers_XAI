"""
wandb_backfill_vit.py
=====================
Retroactively upload finished ViT augmentation-sweep runs to Weights & Biases.

WHY THIS EXISTS
---------------
The ViT runs under
    D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_none/
    D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_plus_original/
    D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_plus_original_1-2/
    D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_random/
were trained on the HPC *without* wandb being initialised. All the
information wandb would have logged live is still on disk, so we can
replay each finished run into a fresh wandb run after the fact.

THE AUGMENTATION CONDITIONS
---------------------------
Each output tree is one augmentation condition. Crucially, the condition
CANNOT always be read from a run's config: the two plus-original sweeps
both carry augment="plus_original" in metrics.json and differ only in
aug_copies (1 vs 2). So each directory is given an explicit label here:

    none            no augmentation
    plus_original   originals + augmented copies, 1:1 ratio (aug_copies=1)
    plus_original2  originals + augmented copies, 1:2 ratio (aug_copies=2)
    random          random augmentation

The label is what is used for the wandb run name and the augmentation tag,
so the 1:1 and 1:2 sweeps stay distinct in the UI.

WHAT EACH FINISHED RUN DIRECTORY CONTAINS
-----------------------------------------
    <root>/ViT_B_mae75/<task>/seed_<N>/<strategy>/
        metrics.json        - config + final test_metrics + test_diagnostics
        train_log.csv        - per-epoch: epoch,lr,train_loss,train_acc,
                                          val_loss,val_acc
        test_predictions.csv - y_true,y_pred,prob_0,prob_1[,prob_2...]
        best_model.pt / last_checkpoint.pt   <- NOT uploaded (too large)

A directory that has only dataset_manifest.csv is an unfinished run and is
skipped automatically.

WHAT THIS SCRIPT DOES PER FINISHED RUN
--------------------------------------
  1. wandb.init() a fresh run in project "vit_orig_vs_aug", with the full
     metrics.json["config"] as the wandb config.
  2. Replays train_log.csv row-by-row with wandb.log(..., step=epoch) so the
     train/val curves appear exactly as if logged live.
  3. Writes the final test_metrics into run.summary (test/accuracy, ...).
  4. Logs a confusion matrix and the raw test-prediction table.
  5. wandb.finish().

The .pt checkpoints are deliberately NOT uploaded (354 MB - 1 GB each;
tens of GB across the sweep). Flip UPLOAD_BEST_CKPT below if you ever want
the best_model.pt as a wandb artifact.

USAGE  (run inside the `mri` conda env, where wandb is installed)
----------------------------------------------------------------
    wandb login                                  # once, if not already
    conda run -n mri python mri_pipeline/wandb_backfill_vit.py --dry-run
    conda run -n mri python mri_pipeline/wandb_backfill_vit.py
    # upload only one (or some) condition(s):
    conda run -n mri python mri_pipeline/wandb_backfill_vit.py --only plus_original2

NOTE ON RE-RUNNING
------------------
Each invocation creates fresh wandb runs. Running the script twice for the
same condition uploads duplicates - use --only to upload just the new
condition(s). If you need to redo a condition, delete its runs in the
wandb UI first.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import wandb

# ─────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────

# The wandb project all runs are uploaded into.
WANDB_PROJECT = "vit_orig_vs_aug"

# Augmentation condition  ->  output tree on disk.
# The KEY is the label used for wandb run names / tags; it is set
# explicitly (not derived from each run's config) so the 1:1 and 1:2
# plus-original sweeps - which share augment="plus_original" in their
# metrics.json - get distinct, non-colliding names.
OUTPUT_ROOTS = {
    "none":          Path(r"D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_none"),
    "plus_original":  Path(r"D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_plus_original"),       # 1:1 orig:aug
    "plus_original2": Path(r"D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_plus_original_1-2"),   # 1:2 orig:aug
    "random":        Path(r"D:/ADNI_BIDS_project/derivatives/vit_outputs_aug_random"),
}

# Set True to additionally upload best_model.pt as a wandb artifact.
# Off by default - these files are 350 MB - 1 GB each.
UPLOAD_BEST_CKPT = False


# ─────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────

def find_finished_runs(root: Path) -> list[Path]:
    """
    Return the run directories under `root` that look complete.

    A run is "finished" if it has BOTH metrics.json (final results) and
    train_log.csv (per-epoch history). Directories with only
    dataset_manifest.csv are unfinished jobs and are excluded.
    """
    runs = []
    for metrics_path in sorted(root.rglob("metrics.json")):
        run_dir = metrics_path.parent
        if (run_dir / "train_log.csv").is_file():
            runs.append(run_dir)
    return runs


def run_identity(cfg: dict, run_dir: Path, label: str) -> dict:
    """
    Derive the wandb run name / grouping fields for one run.

    `label` is the augmentation-condition label of the directory the run
    came from (e.g. "plus_original2"); it is used verbatim for the run
    name and the augmentation tag, INSTEAD of cfg["augment"], because
    cfg["augment"] cannot distinguish the 1:1 and 1:2 plus-original sweeps.

    task / strategy / seed / model_id come from the config, falling back
    to parsing the directory path for any field an older run may not have.
    Layout: <root>/<model_id>/<task>/seed_<N>/<strategy>/
    """
    parts = run_dir.parts
    task = cfg.get("task") or parts[-3]
    strategy = cfg.get("strategy") or parts[-1]
    seed = cfg.get("seed")
    if seed is None:
        # parts[-2] is "seed_<N>"
        seed = parts[-2].replace("seed_", "")
    model_id = cfg.get("model_id") or parts[-4]

    # Human-readable, unique-per-run name. A non-empty `label` is used as a
    # prefix (e.g. "plus_original2-T1_binary-full_ft-seed0"); when there is
    # no label the name is just "<task>-<strategy>-seed<N>".
    name = f"{task}-{strategy}-seed{seed}"
    if label:
        name = f"{label}-{name}"
    return {
        "name": name,
        "label": label or "",
        "task": str(task),
        "strategy": str(strategy),
        "seed": seed,
        "model_id": str(model_id),
    }


def log_training_history(train_log_csv: Path) -> int:
    """
    Replay train_log.csv into the active wandb run.

    Each CSV row is one epoch; it is logged with step=epoch so the wandb
    x-axis matches the real training epoch. Returns the number of epochs
    logged.
    """
    df = pd.read_csv(train_log_csv)
    for _, row in df.iterrows():
        epoch = int(row["epoch"])
        # Group metrics under train/ and val/ namespaces so wandb shows two
        # tidy sections. "lr" and "epoch" are logged flat for easy x-axes.
        wandb.log(
            {
                "epoch": epoch,
                "lr": float(row["lr"]),
                "train/loss": float(row["train_loss"]),
                "train/acc": float(row["train_acc"]),
                "val/loss": float(row["val_loss"]),
                "val/acc": float(row["val_acc"]),
            },
            step=epoch,
        )
    return len(df)


def log_test_results(run, metrics: dict, run_dir: Path) -> None:
    """
    Write final test metrics into run.summary and log the confusion matrix
    and raw prediction table from test_predictions.csv (if present).
    """
    # Final scalar test metrics -> run summary (test/accuracy, test/auc_roc...)
    for key, value in metrics.get("test_metrics", {}).items():
        run.summary[f"test/{key}"] = value

    # Record which epoch was selected as best, if the config has it.
    best_epoch = metrics.get("config", {}).get("best_epoch")
    if best_epoch is not None:
        run.summary["best_epoch"] = best_epoch

    # Per-class test metrics -> summary (test/<class>/f1, ...).
    per_class = metrics.get("test_diagnostics", {}).get("per_class", {})
    for cls_idx, cls in per_class.items():
        cname = cls.get("name", cls_idx)
        for m in ("precision", "recall", "f1", "support"):
            if m in cls:
                run.summary[f"test/{cname}/{m}"] = cls[m]

    # Confusion matrix + raw predictions from test_predictions.csv.
    pred_csv = run_dir / "test_predictions.csv"
    if pred_csv.is_file():
        preds = pd.read_csv(pred_csv)

        # Class names, in label order, from test_diagnostics.labels
        # (e.g. {"0": "CN", "1": "MCI+AD"}); fall back to integer labels.
        labels_map = metrics.get("test_diagnostics", {}).get("labels", {})
        if labels_map:
            class_names = [labels_map[str(i)] for i in
                           range(len(labels_map))]
        else:
            class_names = [str(c) for c in
                           sorted(preds["y_true"].unique())]

        wandb.log({
            "test/confusion_matrix": wandb.plot.confusion_matrix(
                y_true=preds["y_true"].tolist(),
                preds=preds["y_pred"].tolist(),
                class_names=class_names,
            )
        })

        # Raw per-sample predictions as a browsable wandb table.
        wandb.log({"test/predictions": wandb.Table(dataframe=preds)})


# ─────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="List the runs that would be uploaded, then exit "
                         "without contacting wandb.")
    ap.add_argument("--project", default=WANDB_PROJECT,
                    help=f"wandb project name (default: {WANDB_PROJECT}).")
    ap.add_argument("--only", nargs="+", metavar="LABEL",
                    help="Upload only the named augmentation condition(s), "
                         f"e.g. --only plus_original2. "
                         f"Choices: {list(OUTPUT_ROOTS)}. Default: all.")
    ap.add_argument("--root", metavar="DIR",
                    help="Upload a single arbitrary output tree instead of "
                         "the built-in augmentation sweep. Use with "
                         "--project (and optionally --label) for one-off "
                         "uploads, e.g. the viscode2 longitudinal sweep.")
    ap.add_argument("--label", default="", metavar="LABEL",
                    help="Run-name prefix / tag for runs from --root. If "
                         "omitted, runs are named <task>-<strategy>-seed<N> "
                         "with no prefix.")
    args = ap.parse_args()

    # ── Resolve which output tree(s) to upload ───────────────────────────
    # `selected` maps label -> root directory.
    if args.root:
        # Single-directory mode: one ad-hoc tree; ignores OUTPUT_ROOTS/--only.
        root = Path(args.root)
        if not root.is_dir():
            sys.exit(f"[error] --root directory not found: {root}")
        selected = {args.label: root}
    else:
        selected = OUTPUT_ROOTS
        if args.only:
            unknown = [lbl for lbl in args.only if lbl not in OUTPUT_ROOTS]
            if unknown:
                sys.exit(f"[error] unknown --only label(s): {unknown}. "
                         f"Known labels: {list(OUTPUT_ROOTS)}")
            selected = {lbl: OUTPUT_ROOTS[lbl] for lbl in args.only}

    # ── Discover every finished run across the selected output trees ─────
    # all_runs holds (label, run_dir) pairs so each run keeps the label of
    # the directory it came from.
    all_runs: list[tuple[str, Path]] = []
    for label, root in selected.items():
        if not root.is_dir():
            print(f"[warn] output root not found, skipping: {root}")
            continue
        found = find_finished_runs(root)
        print(f"[scan] {label or root.name}: {len(found)} finished run(s)")
        all_runs.extend((label, rd) for rd in found)

    if not all_runs:
        print("[done] no finished runs found - nothing to upload.")
        return

    print(f"\n[total] {len(all_runs)} finished run(s) to upload "
          f"-> wandb project '{args.project}'\n")

    # ── Dry run: just show what would happen ─────────────────────────────
    if args.dry_run:
        for label, run_dir in all_runs:
            metrics = json.loads((run_dir / "metrics.json").read_text())
            ident = run_identity(metrics.get("config", {}), run_dir, label)
            print(f"  would upload: {ident['name']}")
        print("\n[dry-run] no data sent. Re-run without --dry-run to upload.")
        return

    # ── Upload each run ──────────────────────────────────────────────────
    n_ok = 0
    for label, run_dir in all_runs:
        metrics = json.loads((run_dir / "metrics.json").read_text())
        cfg = metrics.get("config", {})
        ident = run_identity(cfg, run_dir, label)

        # One fresh wandb run. group=<task> and the tags let you slice the
        # original-vs-augmented comparison in the UI any way you like.
        run = wandb.init(
            project=args.project,
            name=ident["name"],
            group=ident["task"],
            job_type=ident["strategy"],
            tags=[t for t in (ident["label"], ident["strategy"],
                              ident["model_id"]) if t],
            config=cfg,
            reinit=True,
        )

        try:
            n_epochs = log_training_history(run_dir / "train_log.csv")
            log_test_results(run, metrics, run_dir)

            if UPLOAD_BEST_CKPT:
                ckpt = run_dir / "best_model.pt"
                if ckpt.is_file():
                    art = wandb.Artifact(f"{ident['name']}-best", type="model")
                    art.add_file(str(ckpt))
                    run.log_artifact(art)

            print(f"  [ok] {ident['name']}  ({n_epochs} epochs)")
            n_ok += 1
        finally:
            run.finish()

    print(f"\n[done] uploaded {n_ok}/{len(all_runs)} run(s) "
          f"to wandb project '{args.project}'.")


if __name__ == "__main__":
    sys.exit(main())
