r"""
26b_build_delta_sequences.py   (LOCAL, env: snp)
================================================
Build the delta-pool sequence bundle for the 4 newly-rescued SNPs in
the missingness-resolved (--geno 0.05) PRS bundle. These are the only
rsIDs not already in the legacy fm_embeddings_short_seq_1kb/ npz files;
extracting just these on Colab is ~10x cheaper than re-running all 51.

Delta rsIDs:
  rs7401792    Bellenguez chr14   (--geno 0.05 rescue)
  rs117618017  Bellenguez+Wightman chr15   (--geno 0.05 rescue)
  rs4266886    Desikan PGS chr1   (published-PRS missingness-bypass exception)
  rs115186657  Wightman chr2   (published-PRS MAF-bypass exception)

For rsIDs in multiple sources we keep Bellenguez's β (B > D > W > K
priority), matching the 25e harmonisation; the sequence itself doesn't
depend on β — but the snps.tsv emitted here is consumed by script 26
which expects the standard 12-column schema.

Pipeline:
  1. Read the 4 rsIDs from the per-source TSVs in
     `GWAS_filtered_missingness_resolved_geno005/`.
  2. Pick the priority-winning row per rsID (B > D > W > K).
  3. Write `delta_pool_geno005_snps.tsv` to
     `diff_attn_drive_upload/delta_pool_geno005/`.
  4. Subprocess-invoke `26_build_fm_snp_sequences.py
     --set delta_pool_geno005 --gwas-root diff_attn_drive_upload
     --emit-rc` to produce `fm_sequences/delta_pool_geno005_snp_sequences.tsv`
     (16 rows: 4 SNPs × {REF fwd, ALT fwd, REF rc, ALT rc}) and
     `…_snp_manifest.tsv`.
  5. Bundle ready to upload to Drive.

Usage:
  conda run -n snp python snp_pipeline/26b_build_delta_sequences.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
GENO005_ROOT = BASE / "GWAS_filtered_missingness_resolved_geno005"
UPLOAD_ROOT = BASE / "diff_attn_drive_upload"
SCRIPT_26 = Path(__file__).with_name("26_build_fm_snp_sequences.py")

DELTA_RSIDS = ("rs7401792", "rs117618017", "rs4266886", "rs115186657")

# Priority order — B > D > W > K (matches 25d.PRIORITY)
SOURCE_LOOKUP_ORDER = (("Bellenguez", "sourceB"),
                       ("Desikan", "sourceD"),
                       ("Wightman", "sourceW"),
                       ("Kunkle", "sourceK"))

SNP_COLS = ["rsID", "gene", "lead_rsID", "source", "origin", "CHR",
            "BP_GRCh38", "effect_allele", "other_allele", "REF", "ALT",
            "beta_A1"]


def _row_from_source(rsid: str, src_dir: str) -> dict | None:
    p = GENO005_ROOT / src_dir / f"{src_dir}_snps.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", dtype=str)
    hit = df[df["rsID"] == rsid]
    if len(hit) == 0:
        return None
    r = hit.iloc[0]
    return {c: r[c] for c in SNP_COLS if c in df.columns}


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set-name", default="delta_pool_geno005",
                    help="Name of the synthetic set folder under --upload-root")
    ap.add_argument("--upload-root", type=Path, default=UPLOAD_ROOT)
    args = ap.parse_args()

    out_dir = args.upload_root / args.set_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Stage 1: build the synthetic snps.tsv ────────────────────────────
    rows = []
    for rsid in DELTA_RSIDS:
        row = None
        for src, src_dir in SOURCE_LOOKUP_ORDER:
            row = _row_from_source(rsid, src_dir)
            if row is not None:
                break
        if row is None:
            raise SystemExit(f"[ERROR] {rsid} not found in any source TSV "
                             f"under {GENO005_ROOT}")
        rows.append(row)

    snps_df = pd.DataFrame(rows, columns=SNP_COLS)
    snps_path = out_dir / f"{args.set_name}_snps.tsv"
    snps_df.to_csv(snps_path, sep="\t", index=False)
    print(f"[snps] {snps_path}  ({len(snps_df)} SNPs)")
    print(snps_df[["rsID", "gene", "source", "CHR", "BP_GRCh38",
                    "effect_allele", "other_allele", "beta_A1"]].to_string(
        index=False))

    # ── Stage 2: invoke 26 to build sequences ────────────────────────────
    print(f"\n[26] building 1001-bp REF/ALT fwd+rc sequences …")
    cmd = [sys.executable, str(SCRIPT_26),
           "--set", args.set_name,
           "--gwas-root", str(args.upload_root),
           "--emit-rc"]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=False, text=True)
    if r.returncode != 0:
        raise SystemExit(f"script 26 failed (rc={r.returncode})")

    # ── Stage 3: report bundle contents ──────────────────────────────────
    seq_dir = out_dir / "fm_sequences"
    print(f"\n[bundle] {out_dir}")
    for p in sorted(out_dir.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(out_dir)}  ({p.stat().st_size:,} B)")


if __name__ == "__main__":
    main()
