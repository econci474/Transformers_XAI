r"""
28_run_diff_attention_grid.py   (LOCAL, env: snp)
=================================================
**One-shot driver** for the full diff-attention experiment grid. Runs
every {snp_set × model × aggregation × mlp_depth × delta × seed} cell as
an isolated subprocess, then auto-renders the 28b summary table.

Defaults (the spec the user requested):
    snp_sets      9W27B  9W27B5D  APOE_cluster  APOE_causal     (4)
    models        bmfm_ref bmfm_snp ntv2 caduceus_ph caduceus_ps (5)
    aggregations  global_attn  chrom_hier                       (2)
    mlp_depths    2L  3L                                        (2)
    deltas        off  on                                       (2)
    seeds         0  1  2                                       (3)
    → 4 × 5 × 2 × 2 × 2 × 3 = 480 runs (skip-if-`metrics.json`-exists
    unless --force; each is seconds on CPU).

One-command usage:
    conda run -n snp python snp_pipeline/28_run_diff_attention_grid.py
        # defaults: --output-root snp_pipeline/outputs/diff_attn
        # auto-renders the 28b summary at the end

Useful flags:
    --dry-run               print the plan, run nothing
    --smoke                 fast budget per cell (epochs=3) for end-to-end check
    --no-render             skip the auto-summary at the end
    --force                 re-run cells whose metrics.json already exists
    --snp-sets / --models / --aggregations / --mlp-depths / --deltas / --seeds
                            restrict the grid (any axis)
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import json
import subprocess
import sys
import time
from itertools import product
from pathlib import Path

import pandas as pd

CELL = Path(__file__).parent / "28_train_diff_attention_classifier.py"
RENDER = Path(__file__).parent / "28b_render_diff_attention_summary.py"
DEFAULT_OUT = Path(__file__).parent / "outputs" / "diff_attn"


def _read_metrics(rdir: Path) -> dict | None:
    p = rdir / "metrics.json"
    if not p.exists():
        return None
    try:
        return json.load(open(p))
    except Exception:
        return None


def _fmt_metric(m: dict | None, sp: str, k: str) -> str:
    if not m:
        return "?"
    v = (m.get(sp) or {}).get(k)
    return f"{v:.3f}" if isinstance(v, (int, float)) else "?"


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snp-sets", nargs="+",
                    default=["9W27B", "9W27B5D",
                             "APOE_cluster", "APOE_causal"])
    ap.add_argument("--models", nargs="+",
                    default=["bmfm_ref", "bmfm_snp", "ntv2",
                             "caduceus_ph", "caduceus_ps"])
    ap.add_argument("--aggregations", nargs="+",
                    default=["global_attn", "chrom_hier"])
    ap.add_argument("--mlp-depths", nargs="+", default=["2L", "3L"])
    ap.add_argument("--deltas", nargs="+", default=["off", "on"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    ap.add_argument("--snp-encoder", action="store_true")
    ap.add_argument("--output-root", type=Path, default=DEFAULT_OUT,
                    help=f"default: {DEFAULT_OUT}")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="passthrough args to the cell script "
                         "(e.g. --extra --wandb)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--no-render", action="store_true",
                    help="skip the 28b summary render at the end")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    grid = list(product(args.snp_sets, args.models, args.aggregations,
                        args.mlp_depths, args.deltas, args.seeds))
    N = len(grid)
    print(f"[grid] {N} cells "
          f"({len(args.snp_sets)} sets × {len(args.models)} models × "
          f"{len(args.aggregations)} agg × {len(args.mlp_depths)} depth × "
          f"{len(args.deltas)} δ × {len(args.seeds)} seeds)  "
          f"→ {args.output_root}", flush=True)
    if args.dry_run:
        for g in grid:
            print("  ", g)
        return

    args.output_root.mkdir(parents=True, exist_ok=True)
    done = fail = skip = 0
    headline = []                  # rows for top-N preview at the end
    t0 = time.time()
    for i, (ss, mdl, agg, dep, dl, sd) in enumerate(grid, 1):
        rdir = (args.output_root / ss / mdl / agg / dep /
                f"delta_{dl}" / f"seed_{sd}")
        tag = f"{ss}/{mdl}/{agg}/{dep}/δ{dl}/s{sd}"
        if (rdir / "metrics.json").exists() and not args.force:
            skip += 1
            m = _read_metrics(rdir)
            print(f"[{i:>3}/{N}] SKIP {tag:54s}  "
                  f"test bAcc={_fmt_metric(m,'test','balanced_accuracy')} "
                  f"AUC={_fmt_metric(m,'test','auc')}", flush=True)
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
        t = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        dt = time.time() - t
        m = _read_metrics(rdir)
        ok = r.returncode == 0 and (rdir / "metrics.json").exists()
        if ok:
            done += 1
            print(f"[{i:>3}/{N}] ok   {tag:54s}  "
                  f"test bAcc={_fmt_metric(m,'test','balanced_accuracy')} "
                  f"AUC={_fmt_metric(m,'test','auc')} "
                  f"val bAcc={_fmt_metric(m,'val','balanced_accuracy')} "
                  f"({dt:.0f}s)", flush=True)
        else:
            fail += 1
            print(f"[{i:>3}/{N}] FAIL {tag:54s}  rc={r.returncode} "
                  f"({dt:.0f}s)", flush=True)
            # dump the last few lines of stderr for debugging
            err = (r.stderr or "").strip().splitlines()
            for line in err[-8:]:
                print(f"        {line}")
    dt = time.time() - t0
    print(f"\n[grid] {done} ok, {skip} skipped, {fail} failed in {dt:.0f}s "
          f"({dt/max(1,done):.1f}s/cell ok-avg)", flush=True)

    if args.no_render or done + skip == 0:
        print("[render] skipped (--no-render or no usable cells)")
        return

    # ── auto-render 28b ──────────────────────────────────────────────────
    print(f"\n[render] running 28b summary on {args.output_root} …",
          flush=True)
    rr = subprocess.run(
        [sys.executable, str(RENDER), "--root", str(args.output_root)],
        capture_output=True, text=True)
    print(rr.stdout)
    if rr.returncode:
        print("[render] FAIL:", rr.stderr[-1500:])
        return

    # ── headline: top rows by test balanced_accuracy ─────────────────────
    flat = args.output_root / "diff_attention_summary.csv"
    if flat.exists():
        df = pd.read_csv(flat)
        df = df.sort_values("balanced_accuracy", ascending=False).head(12)
        with pd.option_context("display.width", 220,
                                "display.max_columns", None):
            print("\n=== TOP-12 configurations by TEST balanced_accuracy ===")
            print(df[["snp_set", "model_name", "aggregation_mode",
                      "mlp_depth", "delta", "seed", "balanced_accuracy",
                      "AUROC", "F1", "sensitivity", "specificity"]
                    ].to_string(index=False))
    print(f"\nAll outputs under: {args.output_root}")
    print(f"  per-run trees   : <snp_set>/<model>/<agg>/<depth>/delta_*/seed_*/")
    print(f"  summary CSVs    : summary_long.csv, summary_mean_std.csv,")
    print(f"                    diff_attention_summary.csv")
    print(f"  rendered table  : diff_attention_table.png, .md")


if __name__ == "__main__":
    main()
