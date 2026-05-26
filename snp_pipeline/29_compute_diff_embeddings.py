r"""
29_compute_diff_embeddings.py   (LOCAL on Windows OR Colab, env: snp/colab)
=========================================================================
For each (model, sequence-length) pair, derive `diff = alt_emb − ref_emb`
for the 128 recover_all_pool rsIDs in canonical order. Saves a unified
`recover_all_pool_snp_diff.npz` per (model × seq_length) ready for the
diff-attention sweep trainer (script 30).

Smart-delta strategy for 1001bp:
  - 79 NEW rsIDs come from `recover_all_pool_delta_snp_embeddings.npz`
    (extracted on Colab 2026-05-25)
  - 49 OVERLAP rsIDs come from any legacy npz under the model's dir
    (5W26B, 9W27B, sourceB, sourceW, sourceD, sourceK, APOE_cluster,
     APOE_causal, etc.) — first hit wins per rsID.

For 101bp:
  - All 128 rsIDs from `recover_all_pool_101_snp_embeddings.npz`
    (extracted on Colab 2026-05-26).

Output schema:
  rsids       : (128,)  object  — canonical order from recover_all_pool_snps.tsv
  snp_idx     : (128,)  int64
  diff_emb    : (128, H) float32 — ALT − REF
  model_id    : str
  hidden_dim  : int
  seq_length  : str    — '1001bp' or '101bp'
  provenance  : (128,) object — per-rsID source ('delta' | '<legacy_set>' | '101bp_full')

Usage:
  conda run -n snp python snp_pipeline/29_compute_diff_embeddings.py \\
      --base D:/ADNI_SNP_Omni2.5M_20140220 \\
      --seq-length 1001bp
  # or
  python snp_pipeline/29_compute_diff_embeddings.py \\
      --base /content/drive/MyDrive/ADNI_SNP \\
      --seq-length 101bp
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODELS = ["bmfm_ref", "bmfm_snp", "ntv2", "caduceus_ph", "caduceus_ps"]

# Legacy dir aliases — on Windows D:\ the caduceus dirs have a "_d256" suffix;
# on Drive they were saved without it. Try both.
LEGACY_DIR_ALIASES = {
    "caduceus_ph": ["caduceus_ph", "caduceus_ph_d256"],
    "caduceus_ps": ["caduceus_ps", "caduceus_ps_d256"],
}


def _resolve_model_dir(parent: Path, model: str) -> Path | None:
    for alias in [model] + LEGACY_DIR_ALIASES.get(model, []):
        d = parent / alias
        if d.exists():
            return d
    return None


def _load_npz_index(p: Path) -> dict[str, tuple[np.ndarray, np.ndarray, dict]]:
    """Open an npz; return {rsID: (ref_emb_row, alt_emb_row, meta)}."""
    z = np.load(p, allow_pickle=True)
    if "rsids" not in z.files:
        return {}
    rsids = [str(x) for x in z["rsids"]]
    ref = z["ref_emb"]
    alt = z["alt_emb"]
    meta = {k: (z[k].item() if z[k].dtype.kind in "OU" and z[k].shape == ()
                else str(z[k]))
            for k in z.files
            if k not in {"rsids", "snp_idx", "ref_emb", "alt_emb"}}
    out = {}
    for i, rs in enumerate(rsids):
        out[rs] = (ref[i], alt[i], meta)
    return out


def compute_diff_1001bp(model: str, delta_parent: Path, legacy_parent: Path,
                          target_rsids: list[str]) -> tuple[np.ndarray, list[str], dict]:
    """Smart-delta over delta npz (delta_parent) + all legacy npzs (legacy_parent)."""
    delta_dir = _resolve_model_dir(delta_parent, model)
    legacy_dir = _resolve_model_dir(legacy_parent, model)
    if delta_dir is None and legacy_dir is None:
        sys.exit(f"[ERROR] no model dir for {model} under either "
                  f"{delta_parent} or {legacy_parent}")

    rs2re_alt: dict[str, tuple[np.ndarray, np.ndarray, str]] = {}  # rsID → (ref, alt, source)
    meta_first = None

    # First-priority: the delta npz (from delta_parent)
    if delta_dir is not None:
        delta_path = delta_dir / "recover_all_pool_delta_snp_embeddings.npz"
        if delta_path.exists():
            for rs, (r, a, meta) in _load_npz_index(delta_path).items():
                rs2re_alt[rs] = (r, a, "delta")
                if meta_first is None:
                    meta_first = meta
            print(f"  [delta] {delta_path}: {len(rs2re_alt)} rsIDs")
        else:
            print(f"  [warn] no delta npz at {delta_path}")

    # Second-priority: scan all other legacy npzs (from legacy_parent)
    if legacy_dir is not None:
        for legacy_npz in sorted(legacy_dir.glob("*_snp_embeddings.npz")):
            if legacy_npz.name == "recover_all_pool_delta_snp_embeddings.npz":
                continue
            set_name = legacy_npz.name.replace("_snp_embeddings.npz", "")
            added = 0
            for rs, (r, a, meta) in _load_npz_index(legacy_npz).items():
                if rs in rs2re_alt:
                    continue
                if rs in target_rsids:
                    rs2re_alt[rs] = (r, a, set_name)
                    added += 1
                    if meta_first is None:
                        meta_first = meta
            if added:
                print(f"  [legacy] {set_name}: +{added} new rsIDs")

    # Assemble diff in target order
    H = None
    rows = []
    sources = []
    missing = []
    for rs in target_rsids:
        if rs not in rs2re_alt:
            missing.append(rs)
            continue
        r, a, src = rs2re_alt[rs]
        rows.append(a - r)
        sources.append(src)
        if H is None:
            H = len(r)
    if missing:
        print(f"  [WARN] {len(missing)} rsIDs missing across delta + legacy: "
              f"{missing[:8]}{'…' if len(missing) > 8 else ''}")
    if not rows:
        sys.exit("[ERROR] no rsIDs resolved")
    diff = np.stack(rows).astype(np.float32)
    return diff, sources, (meta_first or {})


def compute_diff_101bp(model: str, parent: Path, target_rsids: list[str]) -> tuple[np.ndarray, list[str], dict]:
    """All 128 from the single recover_all_pool_101 npz."""
    mdir = _resolve_model_dir(parent, model)
    if mdir is None:
        sys.exit(f"[ERROR] no model dir for {model} under {parent}")
    npz_path = mdir / "recover_all_pool_101_snp_embeddings.npz"
    if not npz_path.exists():
        sys.exit(f"[ERROR] missing {npz_path}")
    idx = _load_npz_index(npz_path)
    print(f"  [101bp] {npz_path.name}: {len(idx)} rsIDs")
    rows = []
    sources = []
    missing = []
    meta_first = None
    for rs in target_rsids:
        if rs not in idx:
            missing.append(rs)
            continue
        r, a, meta = idx[rs]
        rows.append(a - r)
        sources.append("101bp_full")
        if meta_first is None:
            meta_first = meta
    if missing:
        print(f"  [WARN] {len(missing)} rsIDs missing: "
              f"{missing[:8]}{'…' if len(missing) > 8 else ''}")
    diff = np.stack(rows).astype(np.float32)
    return diff, sources, (meta_first or {})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True,
                    help="Root containing fm_embeddings_short_seq_1kb/ (legacy) and "
                         "fm_embeddings_short_seq_101bp/ + (optionally) "
                         "fm_embeddings_short_seq_1kb_with_alphagenome/ (delta + AG).")
    ap.add_argument("--delta-1kb-dir", default="fm_embeddings_short_seq_1kb_with_alphagenome",
                    help="Subdir under --base holding delta + AG npzs (default: "
                         "fm_embeddings_short_seq_1kb_with_alphagenome). Used as both "
                         "delta source AND 1kb output target.")
    ap.add_argument("--legacy-1kb-dir", default="fm_embeddings_short_seq_1kb",
                    help="Subdir under --base holding legacy 49-overlap npzs.")
    ap.add_argument("--seq-length", choices=["1001bp", "101bp", "both"],
                    default="both")
    ap.add_argument("--models", nargs="+", default=MODELS)
    ap.add_argument("--snps-tsv", type=Path, default=None,
                    help="recover_all_pool_snps.tsv; default: "
                         "<base>/diff_attn_drive_upload/recover_all_pool/"
                         "recover_all_pool_snps.tsv")
    args = ap.parse_args()

    snps_tsv = args.snps_tsv or (args.base / "diff_attn_drive_upload"
                                  / "recover_all_pool" / "recover_all_pool_snps.tsv")
    snps = pd.read_csv(snps_tsv, sep="\t", dtype=str)
    target_rsids = snps["rsID"].tolist()
    print(f"[targets] {len(target_rsids)} rsIDs from {snps_tsv}")

    seq_lengths = ["1001bp", "101bp"] if args.seq_length == "both" else [args.seq_length]
    for sl in seq_lengths:
        if sl == "1001bp":
            delta_parent = args.base / args.delta_1kb_dir
            legacy_parent = args.base / args.legacy_1kb_dir
            out_parent = delta_parent
        else:
            delta_parent = args.base / "fm_embeddings_short_seq_101bp"
            legacy_parent = None
            out_parent = delta_parent

        if sl == "1001bp":
            print(f"\n=== {sl} (delta: {delta_parent}, legacy: {legacy_parent}) ===")
        else:
            print(f"\n=== {sl} (parent: {delta_parent}) ===")
        for model in args.models:
            print(f"\n[{model}]")
            if sl == "1001bp":
                diff, sources, meta = compute_diff_1001bp(model, delta_parent,
                                                            legacy_parent, target_rsids)
            else:
                diff, sources, meta = compute_diff_101bp(model, delta_parent, target_rsids)
            mdir = _resolve_model_dir(out_parent, model)
            out_path = mdir / "recover_all_pool_snp_diff.npz"
            np.savez_compressed(
                out_path,
                rsids=np.array(target_rsids, dtype=object),
                snp_idx=np.arange(len(target_rsids), dtype=np.int64),
                diff_emb=diff,
                model_id=meta.get("model_id", model),
                hidden_dim=diff.shape[1],
                seq_length=sl,
                provenance=np.array(sources, dtype=object),
            )
            print(f"  [out] {out_path}  shape=(128, {diff.shape[1]})  "
                  f"size={out_path.stat().st_size:,} B")
            # Provenance breakdown
            from collections import Counter
            src_counts = Counter(sources)
            print(f"  sources: {dict(src_counts)}")


if __name__ == "__main__":
    main()
