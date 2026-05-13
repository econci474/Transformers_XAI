"""
12_extract_patient_genotypes.py
================================
Extract per-patient additive dosage (0/1/2 copies of A1) for the 430 GW-sig
SNPs from the QC'd patient PLINK fileset.

Input
-----
  D:/ADNI_SNP_Omni2.5M_20140220/SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr.{bed,bim,fam}
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_labels.tsv      (label==1 rows = GW-sig)

Output
------
  D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/patient_genotypes/
    label1_rsids.txt              one rsID per line (input to plink --extract)
    patient_dosage.raw            raw plink output (--recode A)
    patient_dosage.tsv            cleaned: PTID × rsID dosage matrix (rows=patients,
                                                                    cols=SNPs,
                                                                    values=0/1/2 or NA)
    patient_dosage.bim            allele table (CHR rsID 0 BP A1 A2) for the 430 SNPs,
                                  copied straight from plink's extract output so
                                  downstream scripts know which allele = A1 (the dosage
                                  is the count of A1 copies).

Conventions
-----------
  - Patient BIM A2 is REF (file has the `_refcorr` suffix; this was set with
    `plink2 --ref-from-fa` per 01_notes.txt).
  - Patient BIM A1 is the alt allele the patient may carry.
  - --recode A produces one column per SNP, named "<rsID>_<A1>", with value
    = number of A1 copies (0, 1, 2, or NA).

Usage
-----
  conda activate snp
  python snp_pipeline/12_extract_patient_genotypes.py
  python snp_pipeline/12_extract_patient_genotypes.py --check   # dry-run counts only
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = Path(r"D:/ADNI_SNP_Omni2.5M_20140220")
BFILE_PREFIX = BASE_DIR / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr"
LABELS_TSV = BASE_DIR / "bmfm_inputs" / "external_gwas_labels.tsv"
OUT_DIR = BASE_DIR / "bmfm_inputs" / "patient_genotypes"

PLINK_BIN_CANDIDATES = ["plink", "plink2", "plink.exe"]


def find_plink() -> str:
    for cmd in PLINK_BIN_CANDIDATES:
        if shutil.which(cmd):
            return cmd
    sys.exit(
        "[ERROR] could not find plink/plink2/plink.exe on PATH.\n"
        f"  Candidates checked: {PLINK_BIN_CANDIDATES}"
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bfile", type=Path, default=BFILE_PREFIX,
                   help="PLINK .bed/.bim/.fam prefix (no extension).")
    p.add_argument("--labels", type=Path, default=LABELS_TSV,
                   help="external_gwas_labels.tsv (must have SNP, label columns).")
    p.add_argument("--outdir", type=Path, default=OUT_DIR,
                   help="Output directory.")
    p.add_argument("--check", action="store_true",
                   help="Print row/col counts and exit without running plink.")
    args = p.parse_args()

    for f, name in [
        (Path(str(args.bfile) + ".bed"), "BED"),
        (Path(str(args.bfile) + ".bim"), "BIM"),
        (Path(str(args.bfile) + ".fam"), "FAM"),
        (args.labels, "labels TSV"),
    ]:
        if not f.exists():
            sys.exit(f"[ERROR] {name} not found: {f}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # 1. Build the rsID list for label==1.
    labels = pd.read_csv(args.labels, sep="\t", low_memory=False)
    sig = labels.loc[labels["label"] == 1, "SNP"].drop_duplicates()
    n_sig = len(sig)
    print(f"label=1 SNPs in {args.labels.name}: {n_sig}")
    if n_sig == 0:
        sys.exit("[ERROR] no label==1 SNPs found.")

    # Sanity-check intersection with patient BIM.
    bim = pd.read_csv(Path(str(args.bfile) + ".bim"), sep=r"\s+",
                      header=None, usecols=[0, 1, 3, 4, 5],
                      names=["CHR", "rsID", "BP", "A1", "A2"], dtype=str,
                      engine="python")
    in_bim = sig.isin(bim["rsID"])
    n_in_bim = int(in_bim.sum())
    print(f"  of which present in patient BIM: {n_in_bim} / {n_sig}")
    if n_in_bim < n_sig:
        missing = sig[~in_bim].tolist()[:5]
        print(f"  [warn] {n_sig - n_in_bim} not in BIM, e.g. {missing}")

    rsids_path = args.outdir / "label1_rsids.txt"
    sig.to_csv(rsids_path, header=False, index=False)
    print(f"Wrote {rsids_path}  ({n_sig} rsIDs)")

    if args.check:
        n_patients = sum(1 for _ in open(Path(str(args.bfile) + ".fam"), encoding="utf-8"))
        print(f"[check] expected output shape: {n_patients} patients × {n_in_bim} SNPs.")
        return

    # 2. Run plink --extract --recode A.
    plink = find_plink()
    raw_prefix = args.outdir / "patient_dosage"
    cmd = [
        plink,
        "--bfile", str(args.bfile),
        "--extract", str(rsids_path),
        "--recode", "A",
        "--keep-allele-order",  # don't let plink flip A1/A2; we trust the refcorr BIM
        "--out", str(raw_prefix),
    ]
    print(f"\nRunning: {' '.join(cmd)}\n")
    r = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        sys.exit(f"[ERROR] plink exited with code {r.returncode}")

    raw_path = raw_prefix.with_suffix(".raw")
    if not raw_path.exists():
        sys.exit(f"[ERROR] expected plink output not found: {raw_path}")

    # 3. Reshape .raw → patient_dosage.tsv (rows=PTID, cols=rsID, values=0/1/2).
    print(f"\nReshaping {raw_path.name} → patient_dosage.tsv")
    df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    # Plink names columns: FID IID PAT MAT SEX PHENOTYPE rs1_A rs2_T ...
    fixed = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    snp_cols = [c for c in df.columns if c not in fixed]
    # Strip the "_<A1>" suffix to leave only the rsID.
    rsid_map = {c: c.rsplit("_", 1)[0] for c in snp_cols}
    out = df[["IID"] + snp_cols].copy().rename(columns={"IID": "PTID", **rsid_map})
    out_path = args.outdir / "patient_dosage.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    # Also save a slim BIM containing just the 430 rows so downstream scripts know
    # which A1/A2 each dosage column maps to.
    bim_slim = bim.merge(sig.rename("rsID"), on="rsID", how="inner")
    bim_slim_path = args.outdir / "patient_dosage.bim"
    bim_slim.to_csv(bim_slim_path, sep="\t", index=False)

    n_pat, n_snp = out.shape[0], out.shape[1] - 1
    miss = out.iloc[:, 1:].isna().sum().sum()
    print(f"\nWrote {out_path}: {n_pat} patients × {n_snp} SNPs  ({miss} missing cells)")
    print(f"Wrote {bim_slim_path}: {len(bim_slim)} SNP allele rows")


if __name__ == "__main__":
    main()
