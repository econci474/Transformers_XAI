"""
03e_build_embeddings_manifest.py
================================
Catalogue + integrity check for the per-patient clinical-encoder embeddings that
03_encoder_finetune.py writes inline (one `embeddings/` dir per run).

Primary role — **manifest builder**:
  - scans `<out_dir>/**/embeddings/embeddings_meta.json`,
  - runs the LEAKAGE provenance check (exported `split` labels match the source
    seed CSV partition, folds disjoint), and
  - writes `<out_dir>/embeddings_manifest.csv` — the catalogue the multimodal
    fusion step reads to discover encoders and join on (Patient_ID, seed, split).

Optional `--harvest`:
  - if any run still has a `best_checkpoint/` (e.g. one trained with --save_model),
    (re)export its embeddings from the saved weights. Most runs keep no weights
    (embeddings are produced inline), so this is usually a no-op.

Imports TASK_CONFIG from 03_encoder_finetune.py via importlib (digit-prefixed
module name), mirroring snp_pipeline/49 importing 46.

Env: clinical (GPU optional — ~616 short forward passes per run are cheap on CPU).

Usage:
  python clinical_pipeline/03e_build_embeddings_manifest.py \
      --out_dir  /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion \
      --data_dir /home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline
  # add --harvest to (re)export from any best_checkpoint/ dirs that exist.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from _encoder_embed_lib import (  # noqa: E402
    build_manifest, build_seed_frame, load_encoder, make_provenance,
    pooled_predict, verify_split_integrity, write_run_embeddings,
)

# ── Import TASK_CONFIG from 03 (module name starts with a digit) ──────────────
_spec = importlib.util.spec_from_file_location(
    "enc_finetune", HERE / "03_encoder_finetune.py")
_ft = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ft)
TASK_CONFIG = _ft.TASK_CONFIG


def _infer_run(run_dir: Path, out_root: Path):
    """run_dir = out_root/<model_slug>/<task>/seed_<N>/<strategy>."""
    rel = run_dir.relative_to(out_root).parts
    model_slug, task, seed_s, strategy = rel[0], rel[1], rel[2], rel[3]
    return model_slug, task, int(seed_s.replace("seed_", "")), strategy


def _passes_filters(model_slug, task, seed, strategy, args):
    return not (
        (args.model and model_slug != args.model) or
        (args.task and task != args.task) or
        (args.seed is not None and seed != args.seed) or
        (args.strategy and strategy != args.strategy)
    )


def harvest(out_root: Path, args):
    """(Re)export embeddings from any best_checkpoint/ dirs (optional path)."""
    import torch
    ckpts = sorted(out_root.rglob("best_checkpoint"))
    print(f"[harvest] found {len(ckpts)} best_checkpoint/ dir(s)")
    for ckpt in ckpts:
        run_dir = ckpt.parent
        try:
            model_slug, task, seed, strategy = _infer_run(run_dir, out_root)
        except Exception:
            print(f"  [skip] cannot parse run path: {run_dir}"); continue
        if not _passes_filters(model_slug, task, seed, strategy, args):
            continue
        if task not in TASK_CONFIG:
            print(f"  [skip] unknown task {task}: {run_dir}"); continue
        emb_npz = run_dir / "embeddings" / "embeddings.npz"
        if emb_npz.exists() and not args.overwrite:
            print(f"  [skip] embeddings exist: {run_dir}"); continue

        task_cfg = TASK_CONFIG[task]
        print(f"  [harvest] {model_slug}/{task}/seed_{seed}/{strategy}")
        model, tok = load_encoder(ckpt, args.device)
        frame = build_seed_frame(args.data_dir, seed, task_cfg)
        emb_pool, emb_head, probs, preds = pooled_predict(
            model, tok, frame["text"].tolist(),
            max_length=args.max_length, batch_size=32, device=args.device)
        prov = make_provenance(model, model_slug, model.config._name_or_path, task,
                               strategy, seed, args.max_length, task_cfg.get("label_map"),
                               "Generated_Text", args.data_dir, source="harvest",
                               checkpoint_path=ckpt)
        write_run_embeddings(run_dir / "embeddings", frame, emb_pool, emb_head,
                             probs, preds, prov)
        del model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", required=True, help="Root of the encoder run outputs.")
    ap.add_argument("--data_dir", default=None,
                    help="Verbose baseline dir with seed_*/{train,val,test}.csv. Optional for "
                         "manifest/verify (each run's own data_dir from meta is used); REQUIRED "
                         "for --harvest.")
    ap.add_argument("--device", default=None)
    ap.add_argument("--max_length", type=int, default=1024)
    ap.add_argument("--harvest", action="store_true",
                    help="(Re)export embeddings from any best_checkpoint/ dirs.")
    ap.add_argument("--overwrite", action="store_true",
                    help="With --harvest, overwrite existing embeddings.")
    ap.add_argument("--model"); ap.add_argument("--task")
    ap.add_argument("--seed", type=int); ap.add_argument("--strategy")
    args = ap.parse_args()

    out_root = Path(args.out_dir)

    if args.harvest:
        if not args.data_dir:
            ap.error("--harvest requires --data_dir")
        harvest(out_root, args)

    # Leakage / provenance check per run. Each run records its OWN data_dir (the split
    # tree it was trained on) — tasks like T4 use a different tree (verbose/baseline_T4),
    # so verify against meta["data_dir"], falling back to the CLI --data_dir.
    issues = []
    for meta_path in sorted(out_root.rglob("embeddings/embeddings_meta.json")):
        run_embed_dir = meta_path.parent
        try:
            meta = json.loads(meta_path.read_text())
        except Exception as e:
            issues.append((str(run_embed_dir), f"unreadable meta: {e}")); continue
        run_data_dir = meta.get("data_dir") or args.data_dir
        if not run_data_dir:
            issues.append((str(run_embed_dir), "no data_dir in meta and no --data_dir given"))
            continue
        ok, msg = verify_split_integrity(run_embed_dir, run_data_dir, meta.get("seed"))
        if not ok:
            issues.append((str(run_embed_dir), msg))

    out_csv, man = build_manifest(out_root)
    print(f"\nManifest → {out_csv}  ({len(man)} run(s))")
    if issues:
        print(f"[WARN] {len(issues)} run(s) FAILED the split-integrity check:")
        for r, msg in issues:
            print(f"  - {r}: {msg}")
    else:
        print("Split-integrity check passed for all runs (no leakage of fold assignment).")


if __name__ == "__main__":
    main()
