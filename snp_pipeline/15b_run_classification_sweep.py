"""
15b_run_classification_sweep.py
================================
Drive the full frozen-embedding classification sweep over

    6 upstreams x 2 label_modes x 2 pools x 6 aggregations x 3 seeds = 432 runs

calling 15_train_classification_head.py once per cell. Resumable: skips any
run whose metrics.json already exists under --output-root.

Run from the Transformers_XAI directory:

    python snp_pipeline/15b_run_classification_sweep.py \
        --output-root classification_results_wBCE --loss weighted_bce
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

# upstream_tag -> embeddings .npz filename (Colab export dir)
EMB_DIR = Path("D:/ADNI_SNP_Omni2.5M_20140220/bmfm_outputs_regression/"
                "BMFM-DNA-REF_Colab/embeddings")
UPSTREAMS = {
    "base_ref":                  "base_ref.npz",
    "without_ukb_ref":           "without_ukb_ref.npz",
    "without_ukb_augmented_ref": "without_ukb_augmented.npz",
    "by_chrom_fr_ref":           "by_chrom_fr_ref.npz",
    "combos_fr_ref":             "combos_fr_ref.npz",
    "combos_singletons_fr_ref":  "combos_singletons_fr_ref.npz",
}
WINDOWS = Path("D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/"
               "patient_window_sequences/windows.tsv")
SPLITS = {
    "bl_multi":     Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                         "no_cdr_stratified/tabular/baseline"),
    "ever_convert": Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                         "no_cdr_stratified_ever_convert/tabular/baseline"),
}
POOLS = ["cls", "mean"]
AGGREGATIONS = [
    "mean_all", "chrom_mean_then_mean", "attn_all",
    "chrom_mean_then_attn", "chrom_attn_then_attn", "concat_chrom_means",
]
SEEDS = [0, 1, 2]
SCRIPT = Path(__file__).with_name("15_train_classification_head.py")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--loss", choices=["bce", "weighted_bce"], default="weighted_bce")
    ap.add_argument("--device", default=None,
                    help="Pass-through to script 15 (default: its auto-detect).")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for tag, fname in UPSTREAMS.items():
        p = EMB_DIR / fname
        if not p.exists():
            sys.exit(f"[fatal] missing embeddings for {tag}: {p}")
    if not WINDOWS.exists():
        sys.exit(f"[fatal] missing windows.tsv: {WINDOWS}")

    grid = [
        (tag, lm, pool, agg, seed)
        for tag in UPSTREAMS
        for lm in SPLITS
        for pool in POOLS
        for agg in AGGREGATIONS
        for seed in SEEDS
    ]
    total = len(grid)
    print(f"Sweep: {total} runs  loss={args.loss}  -> {args.output_root}")

    t0 = time.time()
    done = skipped = failed = 0
    for i, (tag, lm, pool, agg, seed) in enumerate(grid, 1):
        out_dir = args.output_root / lm / tag / pool / agg / f"seed_{seed}"
        metrics = out_dir / "metrics.json"
        tag_str = f"{tag:<26} {lm:<12} {pool:<4} {agg:<22} seed={seed}"
        if metrics.exists():
            skipped += 1
            print(f"[{i:3d}/{total}] {tag_str}  SKIP (exists)")
            continue

        cmd = [
            sys.executable, str(SCRIPT),
            "--embeddings", str(EMB_DIR / UPSTREAMS[tag]),
            "--windows", str(WINDOWS),
            "--splits-root", str(SPLITS[lm]),
            "--upstream-tag", tag,
            "--label-mode", lm,
            "--pool", pool,
            "--aggregation", agg,
            "--seed", str(seed),
            "--loss", args.loss,
            "--output-root", str(args.output_root),
        ]
        if args.device:
            cmd += ["--device", args.device]

        if args.dry_run:
            print(f"[{i:3d}/{total}] {tag_str}  DRY: {' '.join(cmd)}")
            continue

        rt0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        dt = time.time() - rt0
        if r.returncode == 0 and metrics.exists():
            done += 1
            # pull test_auc out of the child's stdout tail for a live readout
            auc = ""
            for line in reversed(r.stdout.splitlines()):
                if "'auc'" in line and "TEST" in line:
                    auc = line.strip()[:90]
                    break
            print(f"[{i:3d}/{total}] {tag_str}  OK  ({dt:4.0f}s)  {auc}")
        else:
            failed += 1
            print(f"[{i:3d}/{total}] {tag_str}  FAIL rc={r.returncode} ({dt:.0f}s)")
            print("  stderr tail:", r.stderr.strip()[-600:])

    el = time.time() - t0
    print(f"\nDone. ok={done} skipped={skipped} failed={failed}  "
          f"total_time={el/60:.1f} min")


if __name__ == "__main__":
    main()
