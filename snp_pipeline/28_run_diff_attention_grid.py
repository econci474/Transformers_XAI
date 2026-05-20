r"""
28_run_diff_attention_grid.py   (LOCAL, env: snp)
=================================================
Driver for 28_train_diff_attention_classifier.py over the experiment grid:
  {snp_set} × {model} × {aggregation} × {mlp_depth} × {delta} × {seed}
Default = 2 × 4 × 2 × 2 × 2 × 3 = 192 runs (each seconds on CPU — the
expensive part is the upstream 26→27 extraction). Each cell is an isolated
subprocess (a crash in one does not kill the grid); skip-if-done unless
--force. Mirrors the 19_compare_frozen_models grid-loop pattern.

Usage:
  conda run -n snp python snp_pipeline/28_run_diff_attention_grid.py \
      --output-root snp_pipeline/outputs/diff_attn
  ... --dry-run            # print the plan, run nothing
  ... --snp-sets 9W27B --models bmfm_ref --smoke   # quick subset
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from itertools import product
from pathlib import Path

CELL = Path(__file__).parent / "28_train_diff_attention_classifier.py"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snp-sets", nargs="+", default=["9W27B", "9W27B5D"])
    ap.add_argument("--models", nargs="+",
                    default=["bmfm_ref", "bmfm_snp", "ntv2",
                             "caduceus_ph", "caduceus_ps"])
    ap.add_argument("--aggregations", nargs="+",
                    default=["global_attn", "chrom_hier"])
    ap.add_argument("--mlp-depths", nargs="+", default=["2L", "3L"])
    ap.add_argument("--deltas", nargs="+", default=["off", "on"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    ap.add_argument("--snp-encoder", action="store_true")
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="passthrough args to the cell script "
                         "(e.g. --extra --diffs-subdir fm_embeddings --wandb)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    grid = list(product(args.snp_sets, args.models, args.aggregations,
                        args.mlp_depths, args.deltas, args.seeds))
    print(f"[grid] {len(grid)} cells "
          f"({len(args.snp_sets)}×{len(args.models)}×"
          f"{len(args.aggregations)}×{len(args.mlp_depths)}×"
          f"{len(args.deltas)}×{len(args.seeds)})")
    if args.dry_run:
        for g in grid:
            print("  ", g)
        return

    done = fail = skip = 0
    t0 = time.time()
    for i, (ss, mdl, agg, dep, dl, sd) in enumerate(grid, 1):
        rdir = (args.output_root / ss / mdl / agg / dep /
                f"delta_{dl}" / f"seed_{sd}")
        if (rdir / "metrics.json").exists() and not args.force:
            skip += 1
            continue
        cmd = [sys.executable, str(CELL),
               "--snp-set", ss, "--model", mdl, "--aggregation", agg,
               "--mlp-depth", dep, "--delta", dl, "--seed", str(sd),
               "--output-root", str(args.output_root)]
        if args.snp_encoder:
            cmd.append("--snp-encoder")
        if args.force:
            cmd.append("--force")
        if args.smoke:
            cmd.append("--smoke")
        cmd += args.extra
        print(f"\n[{i}/{len(grid)}] {ss}/{mdl}/{agg}/{dep}/delta_{dl}/s{sd}",
              flush=True)
        r = subprocess.run(cmd)
        if r.returncode == 0 and (rdir / "metrics.json").exists():
            done += 1
        else:
            fail += 1
            print(f"  [FAIL] rc={r.returncode}")
    dt = time.time() - t0
    print(f"\n[done] {done} ok, {skip} skipped, {fail} failed "
          f"({dt:.0f}s) → {args.output_root}")


if __name__ == "__main__":
    main()
