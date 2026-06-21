"""
09_join_freesurfer_to_mri_master.py
====================================
Join the per-scan FreeSurfer ROI tables (from `08_extract_adni_freesurfer.py`)
onto the post-exclusion MRI master, derive within-subject longitudinal deltas
(Delta from baseline + percent change + annualised slope), and write the
column-extended copy to the `_with_FSL/` sister directory.

Per the project plan:

  Anchor:
    D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/
    VISCODE_2_aligned_extended_post_exclusion/
    master_mri_clinical_matched_viscode2_extended_post_exclusion.csv
    -- 1456 rows, 571 unique subjects, one row per MRI scan.

  Outputs:
    D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/
    VISCODE_2_aligned_extended_post_exclusion_with_FSL/
    master_mri_clinical_matched_viscode2_extended_post_exclusion_with_FSL.csv

    mri_pipeline/outputs/freesurfer_join/coverage_report.csv
    mri_pipeline/outputs/freesurfer_join/coverage_report.png
    mri_pipeline/outputs/freesurfer_join/sample_subject_traces.png

Join key:
    (RID, adni_viscode) on the MRI side  ==  (RID, VISCODE2) on the FS side.
    RID is extracted from Patient_ID ('041_S_4629' -> 4629).

Tie-break when multiple FS rows match (already handled in 08 by pipeline
priority + QC + RUNDATE), so 08's output is already 1-per-(RID, VISCODE2).

Longitudinal change derivation:
    For each subject, sort scans by adni_viscode month-offset.
    Use the FS LONG row with VISCODE2 == 'bl' (or earliest available) as
    within-subject reference. Compute per ROI:
        delta = visit_value - bl_value
        pct   = 100 * delta / bl_value      (NaN if bl_value == 0 or NaN)
        slope = delta / years_from_bl       (NaN if years_from_bl <= 0)
    If the subject has no longitudinal bl row, all delta/pct/slope are NaN.

Usage:
    python clinical_pipeline/09_join_freesurfer_to_mri_master.py
    python clinical_pipeline/09_join_freesurfer_to_mri_master.py --dry-run
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
DEFAULT_IN_MRI = Path(
    r"D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/"
    r"VISCODE_2_aligned_extended_post_exclusion/"
    r"master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
)
DEFAULT_FS_DIR = Path(r"D:/ADNI_BIDS_project/derivatives/freesurfer_unified")
DEFAULT_OUT_DIR = Path(
    r"D:/ADNI_BIDS_project/derivatives/mri_clinical_matched/"
    r"VISCODE_2_aligned_extended_post_exclusion_with_FSL"
)
DEFAULT_OUT_FILE_NAME = (
    "master_mri_clinical_matched_viscode2_extended_post_exclusion_with_FSL.csv"
)
DEFAULT_REPORT_DIR = Path(__file__).resolve().parents[1] / "mri_pipeline" / "outputs" / "freesurfer_join"

# ROI short_names (must match 08's ROI_TABLE).
ROI_THICKNESS = [
    "entorhinal", "inf_temporal", "middle_temporal", "fusiform",
    "parahippocampal", "precuneus", "post_cingulate", "sup_frontal",
    "rost_mid_frontal", "caud_mid_frontal", "lat_orbitofrontal",
    "med_orbitofrontal", "frontal_pole", "pars_opercularis",
    "pars_triangularis", "pars_orbitalis",
]
ROI_VOLUME = ["hippocampus"]

# Jack/Schwarz AD-signature ROIs (mean of bilateral thickness).
AD_SIGNATURE_ROIS = ["entorhinal", "inf_temporal", "middle_temporal", "fusiform"]

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def ptid_to_rid(ptid: str) -> int | None:
    """'041_S_4629' or 'sub-041S4629' -> 4629. Returns None if unparseable."""
    if not isinstance(ptid, str):
        return None
    s = ptid.replace("sub-", "").replace("_", "").upper()
    m = re.search(r"S(\d+)$", s)
    if not m:
        return None
    return int(m.group(1))


def viscode_to_month(viscode: str) -> float:
    """'bl' -> 0, 'm06' -> 6, 'm12' -> 12, 'sc' -> -1 (pre-baseline).
    Returns NaN for unparseable.
    """
    if not isinstance(viscode, str):
        return float("nan")
    v = viscode.strip().lower()
    if v in ("bl", "baseline", "v01"):
        return 0.0
    if v in ("sc", "screening"):
        return -1.0
    m = re.match(r"^m(\d+)$", v)
    if m:
        return float(m.group(1))
    m = re.match(r"^v(\d+)$", v)
    if m:
        return float(m.group(1)) * 6  # rough fallback; v06 ~ m12
    return float("nan")


def rename_fs_block(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Apply `fs_<prefix>_` to FS-derived columns; keep RID + VISCODE2 raw for join."""
    rename = {
        "fs_pipeline": f"fs_{prefix}_pipeline",
        "IMAGEUID":    f"fs_{prefix}_imageuid",
        "EXAMDATE":    f"fs_{prefix}_examdate",
        "RUNDATE":     f"fs_{prefix}_rundate",
        "OVERALLQC":   f"fs_{prefix}_overallqc",
        "TEMPQC":      f"fs_{prefix}_tempqc",
        "FRONTQC":     f"fs_{prefix}_frontqc",
        "PARQC":       f"fs_{prefix}_parqc",
        "HIPPOQC":     f"fs_{prefix}_hippoqc",
        "LHIPQC":      f"fs_{prefix}_lhipqc",
        "RHIPQC":      f"fs_{prefix}_rhipqc",
        "st10cv":      f"fs_{prefix}_st10cv",
    }
    out = df.rename(columns={k: v for k, v in rename.items() if k in df.columns}).copy()
    # ROI cols already friendly-named (e.g. 'entorhinal_L_thickness'); add prefix.
    roi_cols = [c for c in out.columns
                if (any(c.startswith(f"{r}_") for r in ROI_THICKNESS + ROI_VOLUME))]
    out = out.rename(columns={c: f"fs_{prefix}_{c}" for c in roi_cols})
    return out


def add_ad_signature(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Compute Jack/Schwarz AD-signature mean thickness from the renamed cols."""
    cols = []
    for r in AD_SIGNATURE_ROIS:
        cols.append(f"fs_{prefix}_{r}_L_thickness")
        cols.append(f"fs_{prefix}_{r}_R_thickness")
    have = [c for c in cols if c in df.columns]
    if not have:
        df[f"fs_{prefix}_ad_signature_thickness"] = pd.NA
        return df
    df[f"fs_{prefix}_ad_signature_thickness"] = df[have].mean(axis=1, skipna=False)
    return df


def derive_long_deltas(merged: pd.DataFrame, fs_long_full: pd.DataFrame,
                       prefix: str = "long") -> pd.DataFrame:
    """Within-subject deltas relative to the subject's longitudinal-pipeline
    baseline row. Baseline is the EARLIEST EXAMDATE row per RID in the FS
    Long table -- robust to VISCODE alias (the FS Long table labels most
    ADNI1/GO baselines as 'scmri'/'sc' rather than 'bl', so a strict
    VISCODE2 == 'bl' filter misses ~70% of subjects).

    fs_long_full is the FULL fs_long_per_scan table (not just rows that
    matched the MRI master) so we can find bl values even if the MRI master
    doesn't include the bl scan.
    """
    # Build a lookup of baseline values per RID using earliest EXAMDATE.
    fs_long_copy = fs_long_full.copy()
    fs_long_copy["EXAMDATE"] = pd.to_datetime(fs_long_copy["EXAMDATE"], errors="coerce")
    # Drop rows with no EXAMDATE (can't establish ordering).
    fs_long_copy = fs_long_copy.dropna(subset=["EXAMDATE"])
    fs_long_copy = fs_long_copy.sort_values(["RID", "EXAMDATE"])
    bl_rows = fs_long_copy.drop_duplicates("RID", keep="first")

    long_value_cols = [c for c in merged.columns
                       if c.startswith(f"fs_{prefix}_")
                       and (c.endswith("_thickness") or c.endswith("_volume"))]

    # Map each prefixed col back to its raw ST-friendly name in fs_long_full.
    def strip_prefix(c: str) -> str:
        return c[len(f"fs_{prefix}_"):]

    # Make a per-RID bl-value dict for fast lookup.
    bl_lookup: dict[int, dict[str, float]] = {}
    for _, row in bl_rows.iterrows():
        rid = row.get("RID")
        if pd.isna(rid):
            continue
        bl_lookup[int(rid)] = {
            strip_prefix(c): row.get(strip_prefix(c)) for c in long_value_cols
        }
    bl_examdate_lookup = bl_rows.set_index("RID")["EXAMDATE"].to_dict()

    # Per-row Delta/pct/slope.
    rid_arr = merged.get("_rid", pd.Series([pd.NA] * len(merged))).to_numpy()
    examdate_arr = pd.to_datetime(merged[f"fs_{prefix}_examdate"],
                                  errors="coerce").to_numpy()

    for c in long_value_cols:
        base_name = strip_prefix(c)
        bl_vals = np.full(len(merged), np.nan)
        bl_dates = np.full(len(merged), np.datetime64("NaT"), dtype="datetime64[ns]")
        for i, rid in enumerate(rid_arr):
            if rid is None or pd.isna(rid):
                continue
            d = bl_lookup.get(int(rid))
            if not d:
                continue
            v = d.get(base_name)
            bl_vals[i] = float(v) if pd.notna(v) else np.nan
            dt = pd.to_datetime(bl_examdate_lookup.get(int(rid)), errors="coerce")
            bl_dates[i] = np.datetime64("NaT") if pd.isna(dt) else np.datetime64(dt)
        cur = pd.to_numeric(merged[c], errors="coerce").to_numpy()
        delta = cur - bl_vals
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.where(
                (bl_vals != 0) & np.isfinite(bl_vals),
                100.0 * delta / bl_vals,
                np.nan,
            )
        years = (examdate_arr - bl_dates) / np.timedelta64(365, "D")
        years = pd.to_numeric(pd.Series(years), errors="coerce").to_numpy()
        slope = np.where(years > 0, delta / years, np.nan)
        merged[c + "_delta"] = delta
        merged[c + "_pct"]   = pct
        merged[c + "_slope"] = slope
    return merged


def date_fallback_join(merged: pd.DataFrame, fs_pref: pd.DataFrame,
                       prefix: str, max_days: int = 30) -> pd.DataFrame:
    """Backfill rows where fs_<prefix>_pipeline is NaN by matching the
    closest FS row of the same RID whose EXAMDATE is within +-max_days of
    mri_date. Recovers ADNI1-3T rows that have VISCODE2 = NaN in the
    UCSFFSX51_ADNI1_3T table (data quality issue in that snapshot).
    """
    pipeline_col = f"fs_{prefix}_pipeline"
    if pipeline_col not in merged.columns:
        return merged
    unmatched_mask = merged[pipeline_col].isna()
    if not unmatched_mask.any():
        return merged

    date_col = f"fs_{prefix}_examdate"
    if date_col not in fs_pref.columns:
        print(f"  [warn] {date_col} not in fs_pref; skipping date-fallback.",
              file=sys.stderr)
        return merged

    # Columns to fill (all except the join keys we don't keep).
    fill_cols = [c for c in fs_pref.columns if c not in ("RID", "VISCODE2")]
    fs_dt = fs_pref.copy()
    fs_dt[date_col] = pd.to_datetime(fs_dt[date_col], errors="coerce")
    by_rid = {rid: g for rid, g in fs_dt.groupby("RID")}

    n_filled = 0
    mri_date_series = pd.to_datetime(merged["mri_date"], errors="coerce")
    for idx in merged.index[unmatched_mask]:
        rid = merged.at[idx, "_rid"]
        if pd.isna(rid):
            continue
        mri_dt = mri_date_series.iat[merged.index.get_loc(idx)]
        if pd.isna(mri_dt):
            continue
        cands = by_rid.get(int(rid))
        if cands is None or cands.empty:
            continue
        diffs = (cands[date_col] - mri_dt).abs()
        if diffs.isna().all():
            continue
        best_pos = diffs.values.argmin()
        if pd.isna(diffs.iloc[best_pos]) or diffs.iloc[best_pos].days > max_days:
            continue
        # Fill missing cells in `merged` row from the matched cand row.
        # Convert Timestamps to string before assignment to avoid pandas dtype
        # errors when the destination column is string-dtype (pandas >= 2.2
        # is strict about Timestamp -> str column inserts).
        cand_row = cands.iloc[best_pos]
        for c in fill_cols:
            if c in cand_row.index:
                v = cand_row[c]
                if isinstance(v, pd.Timestamp):
                    v = v.strftime("%Y-%m-%d")
                merged.at[idx, c] = v
        n_filled += 1
    print(f"  Date-fallback ({prefix}, <={max_days}d): filled {n_filled} rows")
    return merged


def compute_match_status(merged: pd.DataFrame) -> pd.Series:
    has_xs   = merged["fs_xs_pipeline"].notna()
    has_long = merged["fs_long_pipeline"].notna()
    out = pd.Series(["none"] * len(merged), index=merged.index, dtype="object")
    out[has_xs & has_long] = "both"
    out[has_xs & ~has_long] = "xs_only"
    out[~has_xs & has_long] = "long_only"
    return out


# --------------------------------------------------------------------------- #
# Coverage report
# --------------------------------------------------------------------------- #
def make_coverage_report(extended: pd.DataFrame, report_dir: Path) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for pipe_col, label in [("fs_xs_pipeline", "XS"), ("fs_long_pipeline", "Long")]:
        for vc, g_vc in extended.groupby("adni_viscode"):
            counts = g_vc[pipe_col].fillna("MISS").value_counts().to_dict()
            for pipe, n in counts.items():
                rows.append({"side": label, "VISCODE2": vc, "pipeline": pipe, "n": n})
    cov = pd.DataFrame(rows)
    cov_path = report_dir / "coverage_report.csv"
    cov.to_csv(cov_path, index=False)
    print(f"[ok] Wrote {cov_path}  ({len(cov)} cells)")

    # Heatmap PNG.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, side in zip(axes, ["XS", "Long"]):
            sub = cov[cov["side"] == side].pivot_table(
                index="VISCODE2", columns="pipeline", values="n", fill_value=0)
            if sub.empty:
                ax.set_title(f"{side}: no data"); continue
            sub = sub.reindex(sorted(sub.index, key=lambda v: (
                999 if v in (None, "nan") else
                0 if v.lower() == "bl" else
                int(re.match(r"^m(\d+)", v.lower()).group(1)) if re.match(r"^m(\d+)", v.lower()) else 999)))
            im = ax.imshow(sub.values, aspect="auto", cmap="viridis")
            ax.set_xticks(range(len(sub.columns)))
            ax.set_xticklabels(sub.columns, rotation=45, ha="right")
            ax.set_yticks(range(len(sub.index)))
            ax.set_yticklabels(sub.index)
            ax.set_title(f"{side}: match count by VISCODE2 x pipeline")
            for (i, j), v in np.ndenumerate(sub.values):
                ax.text(j, i, str(int(v)), ha="center", va="center",
                        color="white" if v < sub.values.max() / 2 else "black",
                        fontsize=8)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.suptitle("FreeSurfer match coverage  (n = MRI-master scans matched)")
        fig.tight_layout()
        png_path = report_dir / "coverage_report.png"
        fig.savefig(png_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"[ok] Wrote {png_path}")
    except Exception as e:
        print(f"[warn] Coverage heatmap failed: {e}", file=sys.stderr)


def make_sample_traces(extended: pd.DataFrame, report_dir: Path) -> None:
    """Plot fs_long_entorhinal_L_thickness trajectories for 5 sMCI->AD subjects."""
    if "conversion_group" not in extended.columns:
        print("[warn] No conversion_group col; skipping sample traces.", file=sys.stderr)
        return
    converters = extended[extended["conversion_group"].isin(["pMCI"])].copy()
    if converters.empty:
        print("[warn] No pMCI subjects; skipping sample traces.", file=sys.stderr)
        return
    # Pick 5 with the most longitudinal data.
    counts = converters.groupby("bids_sub").size().sort_values(ascending=False)
    picks = list(counts.head(5).index)
    if not picks:
        return

    col = "fs_long_entorhinal_L_thickness"
    if col not in extended.columns:
        print(f"[warn] No {col} col; skipping sample traces.", file=sys.stderr)
        return

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5))
        for sub in picks:
            df_s = extended[extended["bids_sub"] == sub].copy()
            df_s = df_s.sort_values("mri_date")
            df_s["mri_date"] = pd.to_datetime(df_s["mri_date"], errors="coerce")
            ax.plot(df_s["mri_date"], df_s[col], "o-", label=sub)
        ax.set_xlabel("MRI date")
        ax.set_ylabel(f"{col} (mm)")
        ax.set_title("Sample pMCI subjects: entorhinal L thickness trajectory "
                     "(within-subject longitudinal pipeline)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        png_path = report_dir / "sample_subject_traces.png"
        fig.savefig(png_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"[ok] Wrote {png_path}")
    except Exception as e:
        print(f"[warn] Sample-trace plot failed: {e}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #
def verify_bl_deltas_zero(extended: pd.DataFrame, fs_long_full: pd.DataFrame) -> None:
    """Sanity check: for any merged row whose fs_long_examdate equals the
    subject's FS-Long baseline EXAMDATE (earliest per RID in the FULL FS
    Long table), all fs_long_*_delta columns must be exactly 0.

    Non-baseline rows can (and should) have non-zero deltas.
    """
    fs_full = fs_long_full.copy()
    fs_full["EXAMDATE"] = pd.to_datetime(fs_full["EXAMDATE"], errors="coerce")
    fs_full = fs_full.dropna(subset=["EXAMDATE"])
    bl_examdate_per_rid = fs_full.groupby("RID")["EXAMDATE"].min().to_dict()

    # Compute the row's RID and compare its fs_long_examdate to the baseline.
    ext = extended.copy()
    ext["_rid_tmp"] = ext["Patient_ID"].map(ptid_to_rid)
    ext["fs_long_examdate_dt"] = pd.to_datetime(ext["fs_long_examdate"], errors="coerce")
    ext["_bl_dt"] = ext["_rid_tmp"].map(bl_examdate_per_rid)
    is_baseline_row = (ext["fs_long_examdate_dt"].notna()
                       & ext["_bl_dt"].notna()
                       & (ext["fs_long_examdate_dt"] == ext["_bl_dt"]))
    n_bl = int(is_baseline_row.sum())
    if n_bl == 0:
        print("[verify] No merged rows coincide with FS-Long baseline; skipping check.")
        return
    bl_rows = ext[is_baseline_row]
    delta_cols = [c for c in extended.columns
                  if c.startswith("fs_long_") and c.endswith("_delta")]
    bad = 0
    for c in delta_cols:
        v = pd.to_numeric(bl_rows[c], errors="coerce")
        n_nonzero = ((v != 0) & v.notna()).sum()
        if n_nonzero:
            print(f"  [verify] {c}: {n_nonzero} non-zero baseline-row deltas",
                  file=sys.stderr)
            bad += n_nonzero
    if bad:
        print(f"[verify] {bad} non-zero baseline-row deltas across "
              f"{len(delta_cols)} cols", file=sys.stderr)
    else:
        print(f"[verify] All {len(delta_cols)} baseline-row deltas are 0 "
              f"(across {n_bl} baseline rows in the merged output).")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in",      dest="in_mri", type=Path, default=DEFAULT_IN_MRI)
    p.add_argument("--fs-dir",  type=Path, default=DEFAULT_FS_DIR)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument("--out-name", default=DEFAULT_OUT_FILE_NAME)
    p.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    if not args.in_mri.exists():
        print(f"[ERROR] MRI master not found: {args.in_mri}", file=sys.stderr)
        return 1
    fs_xs_csv   = args.fs_dir / "fs_xs_per_scan.csv"
    fs_long_csv = args.fs_dir / "fs_long_per_scan.csv"
    for f in (fs_xs_csv, fs_long_csv):
        if not f.exists():
            print(f"[ERROR] Missing {f}. Run 08_extract_adni_freesurfer.py first.",
                  file=sys.stderr)
            return 2

    print(f"[09] MRI master:  {args.in_mri}")
    print(f"[09] FS XS:       {fs_xs_csv}")
    print(f"[09] FS Long:     {fs_long_csv}")

    # ---- Load -----------------------------------------------------------
    mri = pd.read_csv(args.in_mri, low_memory=False)
    fs_xs   = pd.read_csv(fs_xs_csv, low_memory=False)
    fs_long = pd.read_csv(fs_long_csv, low_memory=False)
    print(f"  MRI master: {mri.shape}; FS XS: {fs_xs.shape}; FS Long: {fs_long.shape}")

    # ---- RID derivation -------------------------------------------------
    mri["_rid"] = mri["Patient_ID"].map(ptid_to_rid).astype("Int64")
    n_bad = mri["_rid"].isna().sum()
    if n_bad:
        print(f"  [warn] {n_bad} rows with unparseable Patient_ID; will not match.",
              file=sys.stderr)

    # ---- Prepare FS sides ------------------------------------------------
    fs_xs["RID"]      = pd.to_numeric(fs_xs["RID"],   errors="coerce").astype("Int64")
    fs_long["RID"]    = pd.to_numeric(fs_long["RID"], errors="coerce").astype("Int64")
    fs_xs["VISCODE2"]   = fs_xs["VISCODE2"].astype(str).str.strip()
    fs_long["VISCODE2"] = fs_long["VISCODE2"].astype(str).str.strip()

    # Prefix XS / Long block columns (excluding RID + VISCODE2).
    fs_xs_pref   = rename_fs_block(fs_xs,   "xs")
    fs_long_pref = rename_fs_block(fs_long, "long")
    fs_xs_pref   = add_ad_signature(fs_xs_pref,   "xs")
    fs_long_pref = add_ad_signature(fs_long_pref, "long")

    # ---- Left-join (XS) --------------------------------------------------
    merged = mri.merge(
        fs_xs_pref,
        how="left",
        left_on=["_rid", "adni_viscode"],
        right_on=["RID",  "VISCODE2"],
        suffixes=("", "_fs_xs"),
    )
    print(f"  After XS join: {merged.shape}")
    # Drop the FS-side RID/VISCODE2 (already on merged via _rid / adni_viscode).
    merged = merged.drop(columns=[c for c in ("RID", "VISCODE2_fs_xs") if c in merged.columns])
    if "VISCODE2" in merged.columns and "VISCODE_2" in merged.columns:
        merged = merged.drop(columns=["VISCODE2"])
    # Date-based fallback for VISCODE2-NaN rows (ADNI1-3T data quality issue).
    # Exact-day only: empirically every row recovered by the fallback has
    # EXAMDATE == mri_date (diff=0), so a wider window adds risk without
    # adding recall.
    merged = date_fallback_join(merged, fs_xs_pref, "xs", max_days=0)

    # ---- Left-join (Long) ------------------------------------------------
    merged = merged.merge(
        fs_long_pref,
        how="left",
        left_on=["_rid", "adni_viscode"],
        right_on=["RID",  "VISCODE2"],
        suffixes=("", "_fs_long"),
    )
    print(f"  After Long join: {merged.shape}")
    merged = merged.drop(columns=[c for c in ("RID", "VISCODE2_fs_long") if c in merged.columns])
    if "VISCODE2" in merged.columns and "VISCODE_2" in merged.columns:
        merged = merged.drop(columns=["VISCODE2"])
    # Date-based fallback for the longitudinal join too (exact-day only).
    merged = date_fallback_join(merged, fs_long_pref, "long", max_days=0)

    # ---- Within-subject longitudinal deltas -----------------------------
    print("  Deriving within-subject deltas...")
    merged = derive_long_deltas(merged, fs_long, prefix="long")

    # ---- AD-signature delta variants ------------------------------------
    # Composite for delta/pct/slope: mean of bilateral AD-signature ROI deltas.
    for suffix in ("delta", "pct", "slope"):
        cols = []
        for r in AD_SIGNATURE_ROIS:
            cols.append(f"fs_long_{r}_L_thickness_{suffix}")
            cols.append(f"fs_long_{r}_R_thickness_{suffix}")
        have = [c for c in cols if c in merged.columns]
        if have:
            merged[f"fs_long_ad_signature_thickness_{suffix}"] = merged[have].mean(
                axis=1, skipna=False)

    # ---- Match status ----------------------------------------------------
    merged["fs_match_status"] = compute_match_status(merged)
    print("  Match status counts:")
    print(merged["fs_match_status"].value_counts().to_string())

    # ---- Sanity: row count preserved ------------------------------------
    if len(merged) != len(mri):
        print(f"[ERROR] Row count drifted: input {len(mri)} -> output {len(merged)}",
              file=sys.stderr)
        return 3
    print(f"  Row count preserved: {len(merged)}")

    # ---- Verify bl deltas zero ------------------------------------------
    verify_bl_deltas_zero(merged, fs_long)

    # ---- Drop internal --------------------------------------------------
    merged = merged.drop(columns=["_rid"])

    if args.dry_run:
        print("\n[dry-run] Not writing outputs.")
        return 0

    # ---- Write -----------------------------------------------------------
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_dir / args.out_name
    merged.to_csv(out_csv, index=False)
    print(f"\n[ok] Wrote {out_csv}")
    print(f"     shape: {merged.shape}")
    print(f"     cols added: {merged.shape[1] - mri.shape[1]} "
          f"(was {mri.shape[1]}, now {merged.shape[1]})")

    # ---- Coverage report -------------------------------------------------
    make_coverage_report(merged, args.report_dir)
    make_sample_traces(merged, args.report_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
