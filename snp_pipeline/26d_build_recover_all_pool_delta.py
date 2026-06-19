r"""
26d_build_recover_all_pool_delta.py   (LOCAL, env: snp)
======================================================
Smart-delta sequence bundle for Colab FM extraction: filter the full
128-SNP `recover_all_pool` down to only the rsIDs NOT already in any
existing per-model `.npz` under `fm_embeddings_short_seq_1kb/`.

Per [[feedback-smart-delta-fm-extraction]]: embed only NEW rsIDs on
Colab; assemble per-set npz locally from delta + legacy. This saves
~38% Colab wall-clock when extending the embedding pool.

Pre-flight delta (verified 2026-05-25 night):
  recover_all_pool : 128 target rsIDs
  legacy npz union : 49 already embedded across one or more of
                     bmfm_ref / bmfm_snp / ntv2 / caduceus_ph_d256 /
                     caduceus_ps_d256
  DELTA            : 79 NEW rsIDs (need Colab embedding)
  → 79 × {REF, ALT} × {fwd, rc} = **316 sequence rows**

Pipeline:
  1. Load `recover_all_pool_snps.tsv` (128 rows).
  2. Scan all `fm_embeddings_short_seq_1kb/*/*_snp_embeddings.npz` and
     union the `rsids` arrays.
  3. Filter the 128-row tsv to the 79 NEW rsIDs (set difference).
  4. Write `recover_all_pool_delta_snps.tsv` to
     `diff_attn_drive_upload/recover_all_pool_delta/`.
  5. Subprocess-invoke `26_build_fm_snp_sequences.py --set
     recover_all_pool_delta --gwas-root diff_attn_drive_upload --emit-rc`.
  6. Bundle ready to upload to Drive.

Usage:
  conda run -n snp python snp_pipeline/26d_build_recover_all_pool_delta.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
UPLOAD_ROOT = BASE / "diff_attn_drive_upload"
EXISTING_EMBED_ROOT = BASE / "fm_embeddings_short_seq_1kb"
SCRIPT_26 = Path(__file__).with_name("26_build_fm_snp_sequences.py")


def _scan_existing_rsids() -> set[str]:
    """Union of `rsids` arrays across every per-model per-set npz under
    fm_embeddings_short_seq_1kb/."""
    out: set[str] = set()
    for npz in EXISTING_EMBED_ROOT.rglob("*_snp_embeddings.npz"):
        try:
            data = np.load(npz, allow_pickle=True)
            if "rsids" in data:
                for rs in data["rsids"]:
                    out.add(str(rs))
        except Exception as e:
            print(f"  [warn] {npz.name}: {e}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-set", default="recover_all_pool",
                    help="Full set name to derive delta from")
    ap.add_argument("--delta-set", default="recover_all_pool_delta")
    ap.add_argument("--upload-root", type=Path, default=UPLOAD_ROOT)
    args = ap.parse_args()

    # ── Stage 1: read full 128-SNP source set ───────────────────────────
    src_tsv = args.upload_root / args.source_set / f"{args.source_set}_snps.tsv"
    if not src_tsv.exists():
        sys.exit(f"[ERROR] source set TSV not found: {src_tsv}")
    full = pd.read_csv(src_tsv, sep="\t", dtype=str)
    print(f"[full] {src_tsv.name}: {len(full)} SNPs")

    # ── Stage 2: scan legacy embeddings for already-embedded rsIDs ──────
    print(f"[scan] enumerating existing embeddings under {EXISTING_EMBED_ROOT} …")
    legacy = _scan_existing_rsids()
    print(f"[scan] union of rsIDs across all legacy npz files: {len(legacy)}")

    # ── Stage 3: compute delta = full target − legacy ───────────────────
    target = set(full["rsID"])
    delta = sorted(target - legacy)
    overlap = sorted(target & legacy)
    print(f"\n[delta] full target:      {len(target)} rsIDs")
    print(f"[delta] legacy overlap:   {len(overlap)} rsIDs (skip)")
    print(f"[delta] NEW (to embed):   {len(delta)} rsIDs")

    if not delta:
        sys.exit("[done] No new rsIDs to embed — full set already covered.")

    # ── Stage 4: filter + write delta snps.tsv ──────────────────────────
    delta_df = full[full["rsID"].isin(delta)].reset_index(drop=True)
    out_dir = args.upload_root / args.delta_set
    out_dir.mkdir(parents=True, exist_ok=True)
    delta_tsv = out_dir / f"{args.delta_set}_snps.tsv"
    delta_df.to_csv(delta_tsv, sep="\t", index=False)
    print(f"\n[delta-tsv] {delta_tsv} ({len(delta_df)} SNPs)")
    print("First 10 NEW rsIDs by source:")
    print(delta_df[["rsID", "gene", "source", "origin", "CHR",
                     "BP_GRCh38", "effect_allele", "other_allele"]].head(10)
                                                            .to_string(index=False))
    # Also write a no-op patient_dosage placeholder so script 26 finds it
    dos = out_dir / f"{args.delta_set}_patient_dosage.tsv"
    if not dos.exists():
        pd.DataFrame(columns=["PTID"] + list(delta_df["rsID"])).to_csv(
            dos, sep="\t", index=False)
        print(f"[delta-dos] {dos}  (header-only placeholder)")

    # ── Stage 5: invoke script 26 to build sequences ────────────────────
    print(f"\n[26] building 1001-bp REF/ALT fwd+rc sequences for {len(delta_df)} SNPs …")
    cmd = [sys.executable, str(SCRIPT_26),
           "--set", args.delta_set,
           "--gwas-root", str(args.upload_root),
           "--emit-rc"]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=False, text=True)
    if r.returncode != 0:
        sys.exit(f"[ERROR] script 26 failed (rc={r.returncode})")

    # ── Stage 6: report bundle contents ─────────────────────────────────
    print(f"\n[bundle] {out_dir}")
    for p in sorted(out_dir.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(out_dir)}  ({p.stat().st_size:,} B)")

    seq_tsv = out_dir / "fm_sequences" / f"{args.delta_set}_snp_sequences.tsv"
    if seq_tsv.exists():
        sdf = pd.read_csv(seq_tsv, sep="\t", dtype=str)
        n_snps = sdf["rsID"].nunique()
        n_rows = len(sdf)
        seqlens = sdf["sequence"].str.len().value_counts() if "sequence" in sdf.columns else {}
        print(f"\n[verify] {n_rows} sequence rows / {n_snps} unique SNPs")
        if len(seqlens):
            print(f"  sequence length: {dict(seqlens.head())}")

    print(f"\n[ready] upload `{out_dir.name}/` to Drive for Colab run")


if __name__ == "__main__":
    main()
