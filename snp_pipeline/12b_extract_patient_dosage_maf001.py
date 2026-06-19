r"""
12b_extract_patient_dosage_maf001.py   (LOCAL, env: snp)
========================================================
Build a CANONICAL MAF>0.01-filtered dosage matrix covering the union of
all PRS-eligible rsIDs (Bellenguez 27 + Wightman 9 + Kunkle 3 + Desikan
9 + the 430 GW-sig SNPs from external_gwas_labels) via:
  plink --bfile <maf001> --extract <rsids.txt> --recode A --keep-allele-order

Output → D:\ADNI_SNP_Omni2.5M_20140220\patient_genotype_dosage_maf001\:
  patient_dosage.tsv     PTID + one column per rsID (0/1/2 ALT dosage)
  patient_dosage.bim     plink bim subset (rsID CHR BP A1 A2)
  patient_dosage.raw     plink raw output (kept for provenance)
  patient_dosage.log     plink log
  rsid_extract.txt       the extract list used
  extract_summary.txt    per-source counts + missing list

This is the SINGLE canonical dosage source for the supplementary
PRS/LD analyses S2/S3/S4 — replaces the patchwork of
patient_genotypes/patient_dosage.tsv + wfilt_dosage/patient_dosage.tsv +
ad-hoc Kunkle extracts.

Usage:
  conda run -n snp python snp_pipeline/12b_extract_patient_dosage_maf001.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
MAF_BED = BASE / "SNP_filtered_with_mri_rsid_clean_current_GRCh38_refcorr_maf001"
PLINK = BASE / "plink.exe"
LABELS = BASE / "bmfm_inputs" / "external_gwas_labels.tsv"
GBF = BASE / "genetic_baselines_filtered"
GB = BASE / "genetic_baselines"
OUT = BASE / "patient_genotype_dosage_maf001"


def gather_prs_rsids() -> dict[str, set[str]]:
    """Union of all PRS-eligible rsIDs (with source attribution)."""
    sources: dict[str, set[str]] = {}

    # Bellenguez 27 + Wightman in-430 + 430 GW-sig set
    lab = pd.read_csv(LABELS, sep="\t", low_memory=False)
    lab = lab[lab["label"] == 1]
    bell = set(lab.loc[lab["source"].astype(str).str.contains("Bellenguez"),
                        "SNP"].astype(str))
    sources["Bellenguez_27"] = bell
    sources["430_GWsig"] = set(lab["SNP"].astype(str))

    # Wightman in-430 + 3 re-extracted
    allow = pd.read_csv(GBF / "wightman_resolved_allowlist.tsv", sep="\t")
    wight_inarr = set(allow["pipeline_rsID"].astype(str))
    re_ex = pd.read_csv(GBF / "wfilt_enriched_weights.tsv", sep="\t")
    wight_re = set(re_ex.loc[
        (re_ex["novel"] == True) & re_ex["source"].astype(str)
        .str.contains("Wightman"), "rsID"].astype(str))
    sources["Wightman_9"] = wight_inarr | wight_re

    # Desikan 4 in_430 + 5 novel
    des = pd.read_csv(GB / "desikan_pgs_resolved.tsv", sep="\t")
    des_in430 = set(des.loc[des["beta_A1"].notna()
                              & (des["in_430"] == True), "rsID"].astype(str))
    de_novel = pd.read_csv(GBF / "wfilt_desikan_enriched_weights.tsv",
                            sep="\t")
    des_novel = set(de_novel.loc[
        (de_novel["novel"] == True) & de_novel["source"].astype(str)
        .str.contains("Desikan"), "rsID"].astype(str))
    sources["Desikan_9"] = des_in430 | des_novel

    # Kunkle 3 genotyped
    sources["Kunkle_3"] = {"rs7920721", "rs593742", "rs2830500"}
    return sources


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bfile", type=Path, default=MAF_BED)
    ap.add_argument("--plink", type=Path, default=PLINK)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print("[gather] building PRS-source rsID union …")
    sources = gather_prs_rsids()
    for s, st in sources.items():
        print(f"  {s:20s} {len(st):5d} rsIDs")
    union = set().union(*sources.values())
    print(f"  {'UNION':20s} {len(union):5d} unique rsIDs")

    # cross-reference against the MAF BIM
    bim = pd.read_csv(MAF_BED.with_suffix(".bim"), sep="\t", header=None,
                       names=["CHR", "rsID", "cM", "BP", "A1", "A2"],
                       dtype=str)
    in_bim = set(bim["rsID"]) & union
    missing = sorted(union - in_bim)
    print(f"\n[bim] {len(in_bim)}/{len(union)} rsIDs in MAF BIM "
          f"({len(missing)} missing)")
    if missing:
        per_src_miss = {s: sorted(st - in_bim) for s, st in sources.items()}
        for s, m in per_src_miss.items():
            print(f"  missing in {s} ({len(m)}): {m[:8]}{'…' if len(m)>8 else ''}")

    extract_path = args.out / "rsid_extract.txt"
    extract_path.write_text("\n".join(sorted(in_bim)) + "\n", "utf-8")

    print(f"\n[plink] extracting dosage for {len(in_bim)} SNPs …")
    prefix = args.out / "patient_dosage"
    cmd = [str(args.plink),
           "--bfile", str(args.bfile),
           "--extract", str(extract_path),
           "--recode", "A",
           "--keep-allele-order",
           "--out", str(prefix)]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("[plink] STDERR:", r.stderr[-2000:])
        sys.exit(r.returncode)
    print("[plink] OK")

    # Convert .raw → .tsv keyed by PTID
    raw_path = prefix.with_suffix(".raw")
    raw = pd.read_csv(raw_path, sep=r"\s+")
    # plink --recode A columns: FID IID PAT MAT SEX PHENOTYPE rsID_A1...
    # we drop FID/PAT/MAT/SEX/PHENOTYPE, keep IID→PTID
    pheno_cols = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    snp_cols = [c for c in raw.columns if c not in pheno_cols]
    # plink suffixes "_A1": "rs12345_T" — strip to rsID only
    df = raw[["IID"] + snp_cols].rename(columns={"IID": "PTID"})
    df.columns = ["PTID"] + [c.rsplit("_", 1)[0] for c in snp_cols]
    df["PTID"] = df["PTID"].astype(str)
    tsv_path = prefix.with_suffix(".tsv")
    df.to_csv(tsv_path, sep="\t", index=False, na_rep="")
    print(f"[out] {tsv_path}  ({len(df)} patients × {len(snp_cols)} SNPs)")

    # bim subset
    bim_sub = bim[bim["rsID"].isin(in_bim)][
        ["rsID", "CHR", "BP", "A1", "A2"]]
    bim_path = args.out / "patient_dosage.bim"
    bim_sub.to_csv(bim_path, sep="\t", index=False)
    print(f"[out] {bim_path}  ({len(bim_sub)} rows)")

    # summary
    lines = [
        "Canonical MAF>0.01-filtered patient dosage matrix",
        f"  BED source     : {args.bfile.name}",
        f"  Patients       : {len(df)}",
        f"  SNPs extracted : {len(snp_cols)}  (union of all 4 PRS sources "
        "+ the 430 GW-sig set)",
        "",
        "Per-source counts:",
    ]
    for s, st in sources.items():
        in_set = st & in_bim
        miss_set = st - in_bim
        lines.append(f"  {s:20s}  in_bim={len(in_set):5d}  "
                     f"missing={len(miss_set):5d}")
        if miss_set:
            ms = sorted(miss_set)
            lines.append(f"    missing rsIDs: {ms[:12]}"
                         f"{'…' if len(ms)>12 else ''}")
    (args.out / "extract_summary.txt").write_text("\n".join(lines),
                                                    "utf-8")
    print(f"[out] {args.out/'extract_summary.txt'}")


if __name__ == "__main__":
    main()
