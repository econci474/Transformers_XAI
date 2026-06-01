"""
08_extract_adni_freesurfer.py
=============================
Extract per-(RID, VISCODE2) FreeSurfer ROI tables from the ADNI UCSF
FreeSurfer derivative bundle (`D:/ADNI_Freesurfer_metrics.zip`) into two
flat CSVs that downstream `09_join_freesurfer_to_mri_master.py` can consume:

    derivatives/freesurfer_unified/fs_xs_per_scan.csv     # cross-sectional
    derivatives/freesurfer_unified/fs_long_per_scan.csv   # longitudinal-pipeline

Per the project plan (see `.tool/plans/snug-honking-aho.md`):

  Cross-sectional priority cascade (FS6 > FS5.1 > FS5.1 ADNI1 3T):
    UCSFFSX6_07_06_23_31May2026.csv             (FS 6.0)
    UCSFFSX51_08_01_16_31May2026.csv            (FS 5.1)
    UCSFFSX51_ADNI1_3T_05_20_15_31May2026.csv   (FS 5.1, ADNI1 3T)
  Longitudinal priority cascade (FSL51ALL > FSL51Y1):
    UCSFFSL51ALL_05_04_16_31May2026.csv         (template = all timepoints)
    UCSFFSL51Y1_05_04_16_31May2026.csv          (template = bl + Y1)

ROI whitelist (AD-signature + precuneus + posterior cingulate + frontal):
    16 cortical ROIs x 2 hemispheres x TA suffix  +  hippocampus L+R SV
    + ST10CV reference. See ROI_TABLE below.

QC handling:
    OVERALLQC == 'Fail'   -> drop
    OVERALLQC == 'Pass'   -> preferred
    OVERALLQC == 'Partial' or '' -> fallback
    Regional QC (TEMPQC, FRONTQC, HIPPOQC) carried through.

Dedup rule per (RID, VISCODE2):
    1. Drop FAIL rows.
    2. Within each pipeline-source, prefer Pass > Partial, then latest RUNDATE.
    3. Across pipeline-sources, prefer FS6 > FS51 > FS51_ADNI1_3T (XS) and
       FSL51ALL > FSL51Y1 (Long); the higher-priority pipeline wins even if
       a lower-priority row has Pass and the higher has Partial.

Inputs:  D:/ADNI_Freesurfer_metrics.zip  (read directly, no extraction)
Outputs: D:/ADNI_BIDS_project/derivatives/freesurfer_unified/{fs_xs,fs_long}_per_scan.csv

Usage:
    python clinical_pipeline/08_extract_adni_freesurfer.py
    python clinical_pipeline/08_extract_adni_freesurfer.py --dry-run
    python clinical_pipeline/08_extract_adni_freesurfer.py --zip D:/...zip --out-dir D:/...
"""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
DEFAULT_ZIP = Path(r"D:/ADNI_Freesurfer_metrics.zip")
DEFAULT_OUT = Path(r"D:/ADNI_BIDS_project/derivatives/freesurfer_unified")

# Latest snapshot per family (from header peek of the ZIP, 2026-06-01).
XS_SOURCES = [
    # (pipeline-source-label, entry name inside the ZIP)
    ("fs6",            "UCSFFSX6_07_06_23_31May2026.csv"),
    ("fs51",           "UCSFFSX51_08_01_16_31May2026.csv"),
    ("fs51_adni1_3t",  "UCSFFSX51_ADNI1_3T_05_20_15_31May2026.csv"),
]
LONG_SOURCES = [
    ("fs51all", "UCSFFSL51ALL_05_04_16_31May2026.csv"),
    ("fs51y1",  "UCSFFSL51Y1_05_04_16_31May2026.csv"),
]

# Priority order (smaller index = higher priority).
XS_PRIORITY   = {label: i for i, (label, _) in enumerate(XS_SOURCES)}
LONG_PRIORITY = {label: i for i, (label, _) in enumerate(LONG_SOURCES)}

# QC preference: 'Pass' (0) > 'Partial' (1) > '' / NaN (2). 'Fail' -> dropped.
QC_RANK = {"Pass": 0, "Partial": 1, "": 2}

# ROI whitelist: short_name -> (L_st_code, R_st_code, suffix)
# Suffix: TA = cortical thickness average (mm); SV = subcortical volume (mm3).
ROI_TABLE = [
    # short_name,              L,        R,        suffix
    ("entorhinal",            "ST24",   "ST83",   "TA"),
    ("inf_temporal",          "ST32",   "ST91",   "TA"),
    ("middle_temporal",       "ST40",   "ST99",   "TA"),
    ("fusiform",              "ST26",   "ST85",   "TA"),
    ("parahippocampal",       "ST44",   "ST103",  "TA"),
    ("precuneus",             "ST52",   "ST111",  "TA"),
    ("post_cingulate",        "ST50",   "ST109",  "TA"),
    ("sup_frontal",           "ST56",   "ST115",  "TA"),
    ("rost_mid_frontal",      "ST55",   "ST114",  "TA"),
    ("caud_mid_frontal",      "ST15",   "ST74",   "TA"),
    ("lat_orbitofrontal",     "ST36",   "ST95",   "TA"),
    ("med_orbitofrontal",     "ST39",   "ST98",   "TA"),
    ("frontal_pole",          "ST25",   "ST84",   "TA"),
    ("pars_opercularis",      "ST45",   "ST104",  "TA"),
    ("pars_triangularis",     "ST47",   "ST106",  "TA"),
    ("pars_orbitalis",        "ST46",   "ST105",  "TA"),
    ("hippocampus",           "ST29",   "ST88",   "SV"),
]
WHOLE_HEAD_REF = "ST10CV"   # lh aparc total cortical volume; reference

# Identifier + QC columns we keep from every source.
ID_COLS_KEEP = ["RID", "VISCODE2", "EXAMDATE", "IMAGEUID", "RUNDATE"]
QC_COLS_KEEP = ["OVERALLQC", "TEMPQC", "FRONTQC", "PARQC", "HIPPOQC",
                "LHIPQC", "RHIPQC"]  # FS51 has L/R, FS6 has HIPPOQC; both kept.

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def st_columns_to_keep() -> list[str]:
    """All ST* codes we extract (with hemisphere + suffix)."""
    cols = [WHOLE_HEAD_REF]
    for _short, l, r, suf in ROI_TABLE:
        cols.append(f"{l}{suf}")
        cols.append(f"{r}{suf}")
    return cols


def st_col_to_friendly(st_code_with_suf: str) -> str | None:
    """Map e.g. 'ST24TA' -> 'entorhinal_L_thickness'. None if not in whitelist."""
    if st_code_with_suf == WHOLE_HEAD_REF:
        return "st10cv"
    for short, l, r, suf in ROI_TABLE:
        if st_code_with_suf == f"{l}{suf}":
            return f"{short}_L_thickness" if suf == "TA" else f"{short}_L_volume"
        if st_code_with_suf == f"{r}{suf}":
            return f"{short}_R_thickness" if suf == "TA" else f"{short}_R_volume"
    return None


def read_csv_from_zip(zf: zipfile.ZipFile, entry: str) -> pd.DataFrame:
    with zf.open(entry, "r") as fh:
        raw = fh.read()
    return pd.read_csv(io.BytesIO(raw), low_memory=False)


def filter_and_label(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """Drop FAIL rows, keep only whitelisted cols + ID/QC, tag pipeline."""
    if "OVERALLQC" in df.columns:
        df = df[df["OVERALLQC"].fillna("").str.strip().str.lower() != "fail"].copy()

    keep_st = st_columns_to_keep()
    have = set(df.columns)
    keep_id = [c for c in ID_COLS_KEEP if c in have]
    keep_qc = [c for c in QC_COLS_KEEP if c in have]
    keep_st = [c for c in keep_st if c in have]
    missing_st = [c for c in st_columns_to_keep() if c not in have]
    if missing_st:
        print(f"  [warn] {label}: missing {len(missing_st)} ST cols "
              f"(e.g. {missing_st[:3]}); they will be NaN.", file=sys.stderr)
    kept = df[keep_id + keep_qc + keep_st].copy()
    kept["fs_pipeline"] = label
    for c in ID_COLS_KEEP:
        if c not in kept.columns:
            kept[c] = pd.NA
    for c in QC_COLS_KEEP:
        if c not in kept.columns:
            kept[c] = pd.NA
    for c in st_columns_to_keep():
        if c not in kept.columns:
            kept[c] = pd.NA
    # Normalize types.
    kept["RID"]      = pd.to_numeric(kept["RID"], errors="coerce").astype("Int64")
    kept["IMAGEUID"] = pd.to_numeric(kept["IMAGEUID"], errors="coerce").astype("Int64")
    kept["EXAMDATE"] = pd.to_datetime(kept["EXAMDATE"], errors="coerce")
    kept["RUNDATE"]  = pd.to_datetime(kept["RUNDATE"], errors="coerce")
    kept["VISCODE2"] = kept["VISCODE2"].astype("string").str.strip()
    return kept


def dedupe_priority(df: pd.DataFrame, priority: dict[str, int]) -> pd.DataFrame:
    """One row per (RID, VISCODE2). Priority: pipeline > QC > RUNDATE."""
    df = df.copy()
    df["_pri_pipeline"] = df["fs_pipeline"].map(priority).fillna(99).astype(int)
    df["_pri_qc"]       = df["OVERALLQC"].fillna("").map(QC_RANK).fillna(3).astype(int)
    # RUNDATE: prefer later (so use negative epoch for ascending sort)
    df = df.sort_values(
        by=["RID", "VISCODE2", "_pri_pipeline", "_pri_qc", "RUNDATE"],
        ascending=[True, True, True, True, False],
        na_position="last",
    )
    deduped = df.drop_duplicates(subset=["RID", "VISCODE2"], keep="first")
    return deduped.drop(columns=["_pri_pipeline", "_pri_qc"])


def rename_st_to_friendly(df: pd.DataFrame) -> pd.DataFrame:
    """Rename ST* columns to readable names matching the project schema."""
    rename = {}
    for c in df.columns:
        friendly = st_col_to_friendly(c)
        if friendly is not None:
            rename[c] = friendly
    return df.rename(columns=rename)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--zip",     type=Path, default=DEFAULT_ZIP)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    p.add_argument("--dry-run", action="store_true",
        help="Print row counts per source without writing.")
    args = p.parse_args()

    if not args.zip.exists():
        print(f"[ERROR] ZIP not found: {args.zip}", file=sys.stderr)
        return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[08] Reading FreeSurfer bundle: {args.zip}")
    with zipfile.ZipFile(args.zip, "r") as zf:
        names_in_zip = {e.filename: e for e in zf.infolist()}

        # ---- Cross-sectional -------------------------------------------------
        xs_parts = []
        for label, entry in XS_SOURCES:
            if entry not in names_in_zip:
                print(f"  [skip] {label}: not in ZIP -- {entry}", file=sys.stderr)
                continue
            raw = read_csv_from_zip(zf, entry)
            tagged = filter_and_label(raw, label)
            print(f"  [{label}] read {len(raw)} rows, kept {len(tagged)} after QC drop "
                  f"({raw.shape[1]} -> {tagged.shape[1]} cols)")
            xs_parts.append(tagged)
        if not xs_parts:
            print("[ERROR] No cross-sectional sources found in ZIP.", file=sys.stderr)
            return 2
        xs_concat = pd.concat(xs_parts, ignore_index=True)
        xs_dedup  = dedupe_priority(xs_concat, XS_PRIORITY)
        xs_named  = rename_st_to_friendly(xs_dedup)
        print(f"  [XS] concat={len(xs_concat)} -> dedup={len(xs_dedup)} "
              f"({xs_dedup['RID'].nunique()} unique RIDs)")

        # ---- Longitudinal ----------------------------------------------------
        long_parts = []
        for label, entry in LONG_SOURCES:
            if entry not in names_in_zip:
                print(f"  [skip] {label}: not in ZIP -- {entry}", file=sys.stderr)
                continue
            raw = read_csv_from_zip(zf, entry)
            tagged = filter_and_label(raw, label)
            print(f"  [{label}] read {len(raw)} rows, kept {len(tagged)} after QC drop "
                  f"({raw.shape[1]} -> {tagged.shape[1]} cols)")
            long_parts.append(tagged)
        if not long_parts:
            print("[ERROR] No longitudinal sources found in ZIP.", file=sys.stderr)
            return 3
        long_concat = pd.concat(long_parts, ignore_index=True)
        long_dedup  = dedupe_priority(long_concat, LONG_PRIORITY)
        long_named  = rename_st_to_friendly(long_dedup)
        print(f"  [Long] concat={len(long_concat)} -> dedup={len(long_dedup)} "
              f"({long_dedup['RID'].nunique()} unique RIDs)")

    # Coverage table by pipeline source.
    def cov(df: pd.DataFrame, name: str) -> None:
        print(f"\n  [{name}] pipeline coverage:")
        print(df.groupby("fs_pipeline").size().to_string())
        print(f"  [{name}] per VISCODE2 (top 10):")
        print(df["VISCODE2"].value_counts().head(10).to_string())
    cov(xs_named,   "XS")
    cov(long_named, "Long")

    if args.dry_run:
        print("\n[dry-run] Not writing outputs.")
        return 0

    xs_out   = args.out_dir / "fs_xs_per_scan.csv"
    long_out = args.out_dir / "fs_long_per_scan.csv"
    xs_named.to_csv(xs_out, index=False)
    long_named.to_csv(long_out, index=False)
    print(f"\n[ok] Wrote {xs_out}  ({len(xs_named)} rows, {xs_named.shape[1]} cols)")
    print(f"[ok] Wrote {long_out}  ({len(long_named)} rows, {long_named.shape[1]} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
