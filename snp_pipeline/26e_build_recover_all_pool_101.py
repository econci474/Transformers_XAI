r"""
26e_build_recover_all_pool_101.py   (LOCAL, env: snp)
====================================================
Build 101 bp SNP-centred REF/ALT (+ reverse-complement) sequences for
the full 128-SNP `recover_all_pool` set — the shorter context window
variant. Existing legacy embeddings are 1001 bp so NONE can be reused;
all 128 SNPs need fresh extraction at 101 bp.

Pipeline:
  1. Copy `recover_all_pool/recover_all_pool_snps.tsv` to
     `recover_all_pool_101/recover_all_pool_101_snps.tsv`.
  2. Subprocess-invoke `26_build_fm_snp_sequences.py
     --set recover_all_pool_101 --gwas-root diff_attn_drive_upload
     --flank 50 --emit-rc` → 128 × 4 = **512 sequence rows × 101 bp**.

Output → `diff_attn_drive_upload/recover_all_pool_101/`:
  recover_all_pool_101_snps.tsv (128 rows, same as 1001bp)
  fm_sequences/recover_all_pool_101_snp_sequences.tsv (512 rows, 101 bp)
  fm_sequences/recover_all_pool_101_snp_manifest.tsv (128 rows + QC)
  fm_sequences/recover_all_pool_101_seqbuild_{qc.json,summary.txt}

Bundle ready to upload to Drive for Colab 5-model extraction.

Usage:
  conda run -n snp python snp_pipeline/26e_build_recover_all_pool_101.py
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
UPLOAD_ROOT = BASE / "diff_attn_drive_upload"
SCRIPT_26 = Path(__file__).with_name("26_build_fm_snp_sequences.py")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-set", default="recover_all_pool")
    ap.add_argument("--target-set", default="recover_all_pool_101")
    ap.add_argument("--upload-root", type=Path, default=UPLOAD_ROOT)
    ap.add_argument("--flank", type=int, default=50,
                    help="bp each side; window = 1 + 2*flank (default 50 → 101bp)")
    args = ap.parse_args()

    # ── Stage 1: copy source snps.tsv → target dir under target name ────
    src_tsv = args.upload_root / args.source_set / f"{args.source_set}_snps.tsv"
    if not src_tsv.exists():
        sys.exit(f"[ERROR] source set TSV not found: {src_tsv}")

    out_dir = args.upload_root / args.target_set
    out_dir.mkdir(parents=True, exist_ok=True)
    tgt_tsv = out_dir / f"{args.target_set}_snps.tsv"
    shutil.copy2(src_tsv, tgt_tsv)
    nrows = len(pd.read_csv(tgt_tsv, sep="\t"))
    print(f"[copy] {src_tsv.name}  →  {tgt_tsv.name}  ({nrows} SNPs)")

    # Header-only patient_dosage placeholder (script 26 won't drop on absence,
    # but consistent with prior bundles)
    dos = out_dir / f"{args.target_set}_patient_dosage.tsv"
    if not dos.exists():
        snps_df = pd.read_csv(tgt_tsv, sep="\t")
        pd.DataFrame(columns=["PTID"] + list(snps_df["rsID"])).to_csv(
            dos, sep="\t", index=False)
        print(f"[dos]  {dos.name}  (header-only placeholder)")

    # ── Stage 2: invoke script 26 with --flank 50 ────────────────────────
    expected_bp = 1 + 2 * args.flank
    print(f"\n[26] building {expected_bp}-bp REF/ALT fwd+rc sequences for "
          f"{nrows} SNPs …")
    cmd = [sys.executable, str(SCRIPT_26),
           "--set", args.target_set,
           "--gwas-root", str(args.upload_root),
           "--flank", str(args.flank),
           "--emit-rc"]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=False, text=True)
    if r.returncode != 0:
        sys.exit(f"[ERROR] script 26 failed (rc={r.returncode})")

    # ── Stage 3: report bundle contents ─────────────────────────────────
    print(f"\n[bundle] {out_dir}")
    for p in sorted(out_dir.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(out_dir)}  ({p.stat().st_size:,} B)")

    seq_tsv = out_dir / "fm_sequences" / f"{args.target_set}_snp_sequences.tsv"
    if seq_tsv.exists():
        sdf = pd.read_csv(seq_tsv, sep="\t", dtype=str)
        n_snps = sdf["rsID"].nunique()
        n_rows = len(sdf)
        seqlens = (sdf["sequence"].str.len().value_counts()
                    if "sequence" in sdf.columns else {})
        print(f"\n[verify] {n_rows} sequence rows / {n_snps} unique SNPs")
        if len(seqlens):
            print(f"  sequence length: {dict(seqlens.head())}")

    print(f"\n[ready] upload `{out_dir.name}/` to Drive for Colab extraction")


if __name__ == "__main__":
    main()
