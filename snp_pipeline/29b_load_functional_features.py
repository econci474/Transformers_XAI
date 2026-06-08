r"""
29b_load_functional_features.py   (importable + CLI; env: snp / colab)
=====================================================================
Concatenate the 6 AlphaGenome npz outputs into a single (128, 25) feature
matrix aligned to the recover_all_pool canonical rsID order, ready for
per-SNP concat in the diff-attention trainer (script 30v3).

Order (matches the per-modality npz column order; STRICT FILTER ALLOW-LIST
post-2026-06-07 — see 28_extract_alphagenome_scores.py docstring):
  [0..3]   eqtl                — 4 RNA-seq tracks (3 DLPFC + 1 Hippocampus)
  [4..6]   accessibility       — 3 DNase tracks (DLPFC + Hippo + Hippo-astro)
  [7..8]   TSS                 — 2 Hippocampus CAGE tracks
  [9]      splicing_merged     — 1 merged scalar
  [10..12] splicing_components — 3 max-abs (sites, usage, junctions)
  [13..23] chip_histone        — 11 CHIP_HISTONE tracks
                                  (DLPFC × {H3K27ac, H3K27me3, H3K36me3,
                                            H3K4me1, H3K4me3, H3K9ac}
                                   + Layer-of-hippo × {H3K27ac, H3K36me3,
                                                       H3K4me1, H3K4me3,
                                                       H3K9ac})
  [24]     chip_tf             — 1 CHIP_TF track (DLPFC × CTCF)
Total: 4 + 3 + 2 + 1 + 3 + 11 + 1 = 25 features per SNP.

Default column counts come from `effect_scores.shape[1]` of each npz; if
your run produces different counts, the helper auto-pads with zeros to a
fixed width per modality so downstream concat stays deterministic. The
fixed widths are recorded in `FUNCTIONAL_FEATURE_LAYOUT`.

Usage as a library:
  from snp_pipeline import importlib
  fn = importlib.util.spec_from_file_location("ff", "snp_pipeline/29b_load_functional_features.py")
  mod = importlib.util.module_from_spec(fn); fn.loader.exec_module(mod)
  feat, names, dim = mod.load_functional_features("D:/.../fm_embeddings_short_seq_1kb",
                                                    target_rsids)
  # feat: (128, 25) float32

Usage as CLI:
  python 29b_load_functional_features.py --base D:/ADNI_SNP_Omni2.5M_20140220 \\
      --snps-tsv .../recover_all_pool_snps.tsv --out functional_features.npz
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FUNCTIONAL_FEATURE_LAYOUT = [
    ("eqtl",                4),  # 3 DLPFC RNA_SEQ + 1 Hippocampus RNA_SEQ
    ("accessibility",       3),  # DLPFC DNase + Hippo DNase + Hippo-astro DNase
    ("TSS",                 2),  # Hippocampus CAGE ×2
    ("splicing_merged",     1),
    ("splicing_components", 3),
    ("chip_histone",       11),  # DLPFC × 6 marks + Layer-of-hippo × 5 marks
    ("chip_tf",             1),  # DLPFC × CTCF
]
TOTAL_FEATURES = sum(n for _, n in FUNCTIONAL_FEATURE_LAYOUT)  # 25


def _load_modality_matrix(parent: Path, modality: str, target_rsids: list[str],
                            width: int) -> tuple[np.ndarray, list[str]]:
    """Load alphagenome_<modality>/recover_all_pool_snp_embeddings_<modality>.npz;
    return (128, width) matrix aligned to target_rsids (zero-pad if narrower).
    For modality='splicing_components', take the `components` key instead of
    `effect_scores`."""
    npz_path = (parent / f"alphagenome_{modality if modality != 'splicing_components' else 'splicing_merged'}"
                / f"recover_all_pool_snp_embeddings_{modality if modality != 'splicing_components' else 'splicing_merged'}.npz")
    if not npz_path.exists():
        print(f"  [WARN] missing {npz_path} — using zeros")
        return np.zeros((len(target_rsids), width), dtype=np.float32), [f"{modality}_{i}" for i in range(width)]
    z = np.load(npz_path, allow_pickle=True)
    rsids = [str(x) for x in z["rsids"]]
    if modality == "splicing_components":
        if "components" not in z:
            print(f"  [WARN] no 'components' key in {npz_path} — using zeros")
            return np.zeros((len(target_rsids), width), dtype=np.float32), [f"splicing_components_{i}" for i in range(width)]
        X = z["components"].astype(np.float32)
        col_names_z = (list(z["component_names"]) if "component_names" in z
                        else [f"comp_{i}" for i in range(X.shape[1])])
    else:
        X = z["effect_scores"].astype(np.float32)
        col_names_z = (list(z["track_biosamples"])
                        if "track_biosamples" in z and len(z["track_biosamples"]) == X.shape[1]
                        else [f"{modality}_{i}" for i in range(X.shape[1])])

    # Align to target_rsids
    idx_map = {rs: i for i, rs in enumerate(rsids)}
    out = np.zeros((len(target_rsids), width), dtype=np.float32)
    actual_cols = min(X.shape[1], width)
    for i, rs in enumerate(target_rsids):
        if rs in idx_map:
            out[i, :actual_cols] = X[idx_map[rs], :actual_cols]
    # If actual count < width, log
    if X.shape[1] != width:
        print(f"  [info] {modality}: npz has {X.shape[1]} cols, layout expects {width} — zero-pad/truncate")
    col_names = list(col_names_z[:actual_cols]) + [""] * max(0, width - actual_cols)
    return out, col_names


def load_functional_features(parent: Path, target_rsids: list[str]) -> tuple[np.ndarray, list[str]]:
    """Returns ((128, 23) float32, list of 23 column names)."""
    blocks: list[np.ndarray] = []
    col_names: list[str] = []
    for modality, width in FUNCTIONAL_FEATURE_LAYOUT:
        X, cn = _load_modality_matrix(parent, modality, target_rsids, width)
        blocks.append(X)
        col_names.extend([f"{modality}::{n}" for n in cn])
    feat = np.concatenate(blocks, axis=1)
    assert feat.shape == (len(target_rsids), TOTAL_FEATURES), \
        f"expected ({len(target_rsids)}, {TOTAL_FEATURES}), got {feat.shape}"
    return feat, col_names


# ───────────────────────── Modality summaries (for v3 trainer) ─────────────

# Map modality → list of AlphaGenome output_type values that belong to it.
# 2026-06-07 amendment: added chip_histone + chip_tf (CHIP_HISTONE / CHIP_TF
# scorers added to script 28 alongside the strict-filter allow-list). 6
# modalities now drive `load_modality_summaries()` → 12 (abs+signed) keys.
MODALITY_OUTPUT_TYPES = {
    "eqtl":         ["RNA_SEQ"],
    "acc":          ["DNASE", "ATAC"],
    "tss":          ["CAGE", "PROCAP"],
    "splice":       ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"],
    "chip_histone": ["CHIP_HISTONE"],
    "chip_tf":      ["CHIP_TF"],
}


def load_modality_summaries(parent: Path, target_rsids: list[str]) -> dict[str, np.ndarray]:
    """For each of 4 modalities, derive (signed, abs) per-SNP scalar from the
    tidy_long.tsv.gz long-format AlphaGenome output (already brain-filtered).

    signed = max_signed across all tracks within the modality (i.e. value with
    largest |raw_score|, sign preserved).
    abs    = |signed|.

    Returns dict with keys eqtl_signed, eqtl_abs, acc_signed, acc_abs,
    tss_signed, tss_abs, splice_signed, splice_abs — each (N_SNPs,) float32
    aligned to target_rsids order, zero for SNPs absent from tidy_long.
    """
    p = parent / "alphagenome_eqtl" / "recover_all_pool_tidy_long.tsv.gz"
    if not p.exists():
        raise FileNotFoundError(
            f"tidy_long not found at {p} — required for modality summaries. "
            f"Re-run AlphaGenome extraction or copy from local D: drive.")
    df = pd.read_csv(p, sep="\t", compression="gzip",
                      usecols=["rsID", "output_type", "raw_score"])
    df = df[df["output_type"].notna() & df["raw_score"].notna()]
    df["abs_score"] = df["raw_score"].abs()

    rs_idx = {rs: i for i, rs in enumerate(target_rsids)}
    N = len(target_rsids)
    out: dict[str, np.ndarray] = {}
    for modality, ot_list in MODALITY_OUTPUT_TYPES.items():
        sub = df[df["output_type"].isin(ot_list)]
        if sub.empty:
            signed = np.zeros(N, dtype=np.float32)
        else:
            # For each rsID, pick the row with the largest |raw_score|
            idx_max = sub.groupby("rsID")["abs_score"].idxmax()
            best = sub.loc[idx_max, ["rsID", "raw_score"]]
            signed = np.zeros(N, dtype=np.float32)
            for _, row in best.iterrows():
                if row["rsID"] in rs_idx:
                    signed[rs_idx[row["rsID"]]] = np.float32(row["raw_score"])
        out[f"{modality}_signed"] = signed
        out[f"{modality}_abs"] = np.abs(signed).astype(np.float32)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True,
                    help="Root containing fm_embeddings_short_seq_1kb/ (where AlphaGenome npzs live)")
    ap.add_argument("--snps-tsv", type=Path, required=True,
                    help="recover_all_pool_snps.tsv (canonical rsID order)")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output npz path")
    args = ap.parse_args()

    target_rsids = pd.read_csv(args.snps_tsv, sep="\t", dtype=str)["rsID"].tolist()
    print(f"[targets] {len(target_rsids)} rsIDs")

    parent = args.base / "fm_embeddings_short_seq_1kb"
    feat, col_names = load_functional_features(parent, target_rsids)
    print(f"[feat] shape={feat.shape}, cols={len(col_names)}")

    np.savez_compressed(
        args.out,
        rsids=np.array(target_rsids, dtype=object),
        functional=feat,
        column_names=np.array(col_names, dtype=object),
        layout=np.array(FUNCTIONAL_FEATURE_LAYOUT, dtype=object),
    )
    print(f"[out] {args.out}  size={args.out.stat().st_size:,} B")


if __name__ == "__main__":
    main()
