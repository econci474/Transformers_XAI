r"""
28c_render_per_set_tables.py
=============================
For each SNP set under `snp_pipeline/outputs/diff_attn/`, save a per-set
cropped version of the master diff-attention table (the one rendered by
28b at the root). Each per-set PNG/MD contains the same 40 rows for
that set (5 models × 2 agg × 2 depth × 2 δ) in the same house style.

Output per set:
  snp_pipeline/outputs/diff_attn/<set>/diff_attention_table.{png,md}

Reuses `collect`, `summarise`, `render` from 28b verbatim (importlib).

Usage:
  conda run -n snp python snp_pipeline/28c_render_per_set_tables.py \
      --root snp_pipeline/outputs/diff_attn
"""
from __future__ import annotations

import argparse
import importlib.util as _il
from pathlib import Path

import pandas as pd

_S28B = Path(__file__).parent / "28b_render_diff_attention_summary.py"
_spec = _il.spec_from_file_location("s28b", _S28B)
_s28b = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s28b)
collect   = _s28b.collect
summarise = _s28b.summarise
render    = _s28b.render


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()

    long_df = collect(args.root)
    print(f"[load] {len(long_df)} rows from metrics.json tree")
    summ = summarise(long_df)
    sets = sorted(summ["snp_set"].dropna().unique())
    print(f"[sets] {len(sets)}: {sets}")

    for ss in sets:
        sub = summ[summ["snp_set"] == ss].reset_index(drop=True)
        if sub.empty:
            print(f"  [skip] {ss}: no rows")
            continue
        outdir = args.root / ss
        outdir.mkdir(parents=True, exist_ok=True)
        out_png = outdir / "diff_attention_table.png"
        out_md  = outdir / "diff_attention_table.md"
        render(sub, out_png, out_md)
        print(f"  [out] {ss}: {len(sub)} rows -> {out_png.name}")

    print("\nDone.")


if __name__ == "__main__":
    main()
