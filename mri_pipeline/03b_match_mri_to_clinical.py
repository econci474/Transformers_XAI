"""
Script 03b - Match preprocessed MRI scans to clinical visits
=============================================================
Builds a per-scan labelled CSV that pairs every preprocessed MRI volume with the
right row of `master_clinical_verbose.csv`. Also runs a transparent audit of
the two competing matching strategies so the user can see why scans match or
don't.

Why this script exists
----------------------
`mni_staging` (the "first scan per subject" sMRIprep tree on local Windows)
holds one preprocessed T1w per subject in `sub-XYZ/anat/...` with NO session
label in the path. `03_prepare_ViT.py` blindly tags those as `ses-bl`, but
that's wrong for many subjects (the first available MRI for an ADNI3 entrant
is often m60 or later, not bl). The HPC sessionwise tree
(`smriprep_sessionwise/`) has the correct ses-* dirs but we still need a
labels CSV that aligns each preprocessed scan with the clinical row produced
by 03_clinical preprocessing.

Inputs
------
- `bids/imaging/mri_visit_dates_snp_subjects.csv` (SPINE)
                                   Wide MRI dates per subject — one row per
                                   subject, one column per ADNI VISCODE
                                   (bl/m12/m24/...). 616 subjects, ~1,598
                                   non-null scan dates. This is the
                                   authoritative cohort-level scan inventory:
                                   one row per (subject, visit), no duplicate
                                   ImageUIDs.
- `derivatives/clinical/verbose/longitudinal/master_clinical_verbose.csv`
                                   Long-format clinical visits + label cols.
- `metadata/scan_selection.csv`   Optional: provides ImageUID_selected per
                                   (bids_sub, bids_ses) for traceability.
- `metadata/session_map.csv`      Optional: per-ImageUID inventory
                                   (1+ rows per scan session). Only used to
                                   attach ImageUIDs to the master.
- `derivatives/vit_inputs/vit_manifest.csv`
                                   What's actually preprocessed and on disk.

Why the spine is `mri_visit_dates_snp_subjects.csv` and NOT `session_map.csv`:
session_map has one row per ImageUID (e.g. sub-003S4119 m48 2016-08-23 might
have I1293320 + I1293321 — two processed derivatives of the same scan series).
Using session_map as the spine would 9x-inflate per-session counts.
mri_visit_dates_snp_subjects.csv is already deduplicated to one cell per
(subject, visit).

Outputs
-------
**Cohort-level master** (under `derivatives/mri_clinical_matched/`):
- `master_mri_clinical_matched_labels.csv`  One row per MRI scan in the
      SNP+MRI cohort (~4,619 scans across 616 subjects). Source of truth
      for any downstream consumer regardless of which scans are
      preprocessed locally vs on HPC.
      Columns: bids_sub, bids_ses, mri_date, ImageUID, clinical_date,
      days_diff, match_status, VISCODE_long, Label_*.

- `mri_clinical_match_audit.csv`        Per-VISCODE summary of |days_diff|
      median/p95/max, n_matched_viscode_exact, n_matched_within_14d,
      n_matched_within_extended, n_unmatched, recovery rates at multiple
      tolerances.

- `mri_clinical_strategy_comparison.csv`  Strategy A (VISCODE) vs B
      (nearest-by-date +-14d): counts and percentages incl. A OR B union.

- `mri_clinical_date_diffs_per_viscode.png`
      Histogram per ADNI session of |MRI date - clinical date| for the
      direct-VISCODE join. Reveals how tight ADNI's per-VISCODE scheduling
      really is (the protocol says +-14 days, but reality may differ).

**Preprocessed-subset view** (under `derivatives/vit_inputs/`):
- `mri_clinical_matched_labels.csv`     Master filtered to scans actually
      present in vit_manifest.csv. Convenience for 04 fine-tuning.

Match strategies
----------------
A) **Direct VISCODE join** (highest priority)
   Join `(bids_sub, bids_ses_adni)` from MRI side to
   `(Patient_ID, VISCODE_long)` from clinical side. ADNI assigned both labels
   from the same scheduling protocol, so this is the authoritative pairing.
   Pros: zero ambiguity, no leakage. Cons: if a clinical row at that VISCODE
   doesn't exist for the subject (missed assessment), we get no match.

B) **Nearest-by-date within +-14d** (optional fallback for A-unmatched)
   pd.merge_asof(direction="nearest") within +-14d. Used when the ADNI VISCODE
   has no matching clinical row but a clinical visit happened within the
   ADNI-2/3 protocol "+-2 weeks of target date" window. The +-14d cap is
   deliberate: at this window, the clinical visit and MRI are part of the
   same scheduled assessment, so Label_visit_diag is reliable.

DO NOT impute labels by forward/backward fill from a different VISCODE's
clinical row. Label_visit_diag is per-visit (a subject's diagnosis at m06
can differ from m12), so e.g. using m06's labels for an m12 scan would be
incorrect even though the rows belong to the same subject.

Audit-only (do NOT change match_status):
1. Wider date-tolerance recovery: at +-30/60/90/180/365 days, what fraction
   of scans Strategy B would match if we widened the window.
2. Baseline-offset consistency: for each scan, compute
   bl_offset_days = mri_date - subject_bl_clinical_date and check whether
   it lands within +-14d of the canonical ADNI interval for the MRI's ADNI
   VISCODE2 (m12=365d, m24=730d, m36=1095d, ...). This confirms ADNI's
   VISCODE label is consistent with the date, even for scans where the
   matching clinical row doesn't exist.

Final label assigned to each scan
---------------------------------
match_status priority: viscode_exact > nearest_within_14d > unmatched
The `--label-policy` flag controls fallback:
- viscode_only (default):       only A — strictest, no per-visit diagnosis drift
- viscode_then_nearest_14d:     A then B (+-14d nearest)

Usage
-----
    python mri_pipeline/03b_match_mri_to_clinical.py
    python mri_pipeline/03b_match_mri_to_clinical.py --label-policy viscode_then_nearest_14d
    python mri_pipeline/03b_match_mri_to_clinical.py --vit-manifest /alt/path/vit_manifest.csv
"""

import os
import re
import sys
import argparse
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def render_summary_table(summary_df, out_csv, out_png, out_pdf):
    """Render a per-session match-status table in baseline_model_table style.
    Cols: ses, n_scans, matched VISCODE, additionally +-14d nearest,
          unmatched (date-consistent), unmatched (date-inconsistent)."""
    import matplotlib
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    summary_df.to_csv(out_csv, index=False)

    HEADERS = ["Session", "N scans", "matched VISCODE",
               "+ matched ±14d", "unmatched\n(date-consistent)",
               "unmatched\n(date-inconsistent)"]
    cols    = ["ses", "n_scans", "n_matched_viscode",
               "n_added_nearest_14d", "n_unmatched_consistent",
               "n_unmatched_inconsistent"]

    rows = summary_df[cols].astype(str).values.tolist()
    n_rows = len(rows)
    n_cols = len(HEADERS)

    COL_W = [0.85, 0.85, 1.30, 1.20, 1.45, 1.55]   # inches
    ROW_H = 0.34
    HDR_H = 0.55  # header row taller for two-line wraps

    fig_w = 0.30 + sum(COL_W) + 0.30
    fig_h = 0.50 + HDR_H + ROW_H * n_rows + 0.20  # title + header + data + bottom margin

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)

    LEFT = 0.30
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(n_cols)]
    RIGHT = col_left[-1]

    TOP = fig_h - 0.20
    title_y = TOP
    hdr_top = TOP - 0.30
    hdr_bot = hdr_top - HDR_H
    row_tops = [hdr_bot - i * ROW_H for i in range(n_rows + 1)]

    def hline(y, lw=1.0):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw,
                solid_capstyle="butt", zorder=3)

    hline(hdr_top, lw=1.5)   # top border
    hline(hdr_bot, lw=1.0)   # below headers
    hline(row_tops[-1], lw=1.5)  # bottom border

    # Title (two lines: bold main + italic subtitle)
    ax.text((LEFT + RIGHT) / 2, title_y + 0.10,
            "MRI–Clinical Match Status by ADNI Session",
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text((LEFT + RIGHT) / 2, title_y - 0.07,
            "SNP+MRI cohort (n=616 subjects, 1,598 scans)",
            ha="center", va="center", fontsize=8, fontstyle="italic")

    # Headers (italic, multi-line)
    for j, h in enumerate(HEADERS):
        ax.text(col_cx[j], (hdr_top + hdr_bot) / 2, h,
                ha="center", va="center", fontsize=8, fontstyle="italic")

    # Data rows
    for i, row in enumerate(rows):
        y = (row_tops[i] + row_tops[i + 1]) / 2
        for j, val in enumerate(row):
            weight = "bold" if j == 0 else "normal"
            ax.text(col_cx[j], y, val, ha="center", va="center",
                    fontsize=8, fontweight=weight)

    # Outer box
    rect = plt.Rectangle((LEFT, row_tops[-1]),
                          RIGHT - LEFT, hdr_top - row_tops[-1],
                          facecolor="none", edgecolor="black", linewidth=1.5,
                          zorder=5)
    ax.add_patch(rect)

    # Vertical separators (dotted) between column groups
    for x in [col_left[2], col_left[4]]:  # after N scans, after +-14d
        ax.plot([x, x], [row_tops[-1], hdr_top],
                color="black", linewidth=0.6, linestyle=(0, (3, 3)), zorder=2)

    plt.tight_layout(pad=0.1)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def viscode_months(ses: str) -> int:
    """Chronological sort key: 'sc' (-1) -> 'bl' (0) -> 'm{N}' (N).
    Reflects ADNI visit timeline: screening precedes baseline by days-to-weeks.
    Unknown labels sort last."""
    s = str(ses).strip().lower()
    if s == "sc":
        return -1
    if s == "bl":
        return 0
    m = re.match(r"^m(\d+)$", s)
    if m:
        return int(m.group(1))
    return 9999


# ---- Paths -----------------------------------------------------------------
PROJECT_ROOT   = Path(r"D:\ADNI_BIDS_project")
MRI_DATES_CSV  = PROJECT_ROOT / "bids" / "imaging" / "mri_visit_dates_snp_subjects.csv"  # SPINE
MASTER_CLIN    = PROJECT_ROOT / "derivatives" / "clinical" / "verbose" / "longitudinal" / "master_clinical_verbose.csv"
SCAN_SELECTION = PROJECT_ROOT / "metadata" / "scan_selection.csv"     # optional: attach ImageUID
SESSION_MAP    = PROJECT_ROOT / "metadata" / "session_map.csv"        # optional: secondary lookup
SNP_TSV        = PROJECT_ROOT / "bids" / "genotype" / "subjects_with_snp_and_mri.tsv"
VIT_INPUTS     = PROJECT_ROOT / "derivatives" / "vit_inputs"
VIT_MANIFEST   = VIT_INPUTS / "vit_manifest.csv"

# Cohort-level (data-agnostic) outputs — referenced regardless of preprocessing
MATCHED_DIR   = PROJECT_ROOT / "derivatives" / "mri_clinical_matched"
OUT_MASTER    = MATCHED_DIR / "master_mri_clinical_matched_labels.csv"
OUT_AUDIT     = MATCHED_DIR / "mri_clinical_match_audit.csv"
OUT_HIST      = MATCHED_DIR / "mri_clinical_date_diffs_per_viscode.png"
OUT_COMPARE   = MATCHED_DIR / "mri_clinical_strategy_comparison.csv"
OUT_SUMMARY_CSV = MATCHED_DIR / "mri_clinical_match_summary.csv"
OUT_SUMMARY_PNG = MATCHED_DIR / "mri_clinical_match_summary.png"
OUT_SUMMARY_PDF = MATCHED_DIR / "mri_clinical_match_summary.pdf"

# Clinical-aligned variant: bids_ses replaced by VISCODE_long for matched scans
ALIGNED_DIR             = MATCHED_DIR / "clinical_aligned_intermediates"
OUT_ALIGNED_MASTER      = ALIGNED_DIR / "master_mri_clinical_matched_labels.csv"
OUT_ALIGNED_SUMMARY_CSV = ALIGNED_DIR / "mri_clinical_match_summary.csv"
OUT_ALIGNED_SUMMARY_PNG = ALIGNED_DIR / "mri_clinical_match_summary.png"
OUT_ALIGNED_SUMMARY_PDF = ALIGNED_DIR / "mri_clinical_match_summary.pdf"

MATCHED_STATUSES = ("viscode_exact", "nearest_within_14d")

# Preprocessed-subset view (derived from master via vit_manifest filter)
OUT_LABELS    = VIT_INPUTS / "mri_clinical_matched_labels.csv"

LABEL_COLS = [
    "Label_bl_multi", "Label_visit_diag", "Label_visit_granular",
    "Label_1y", "Label_2y", "Label_3y", "Label_4y", "Label_5y",
    "Label_6y", "Label_7y", "Label_8y", "Label_9y", "Label_10y",
]

EXTENDED_TOLERANCES_DAYS = [14, 30, 60, 90, 180, 365]


def normalise_subject(sid: str) -> str:
    """002_S_0413 -> 002S0413; sub-002S0413 -> 002S0413."""
    s = str(sid).strip()
    if s.startswith("sub-"):
        s = s[4:]
    return s.replace("_", "")


def main():
    ap = argparse.ArgumentParser(description="Match preprocessed MRI scans to clinical labels.")
    ap.add_argument("--label-policy",
                    choices=["viscode_only", "viscode_then_nearest_14d"],
                    default="viscode_then_nearest_14d",
                    help="How to assign labels per MRI scan. "
                         "viscode_only: only Strategy A (direct ADNI VISCODE2 join). "
                         "viscode_then_nearest_14d (default): A then nearest clinical visit "
                         "within +-14d (handles e.g. ADNI3 patients whose MRI was taken at the "
                         "screening visit but tagged ses-bl by ADNI). Wider date-based fallbacks "
                         "are not supported because per-visit diagnosis can drift between VISCODEs.")
    ap.add_argument("--vit-manifest", type=str, default=str(VIT_MANIFEST),
                    help="Override path to vit_manifest.csv (default: derivatives/vit_inputs/vit_manifest.csv).")
    ap.add_argument("--no-cohort-filter", action="store_true",
                    help="Audit ALL ADNI subjects in session_map (default: only SNP+MRI cohort, n=616).")
    args = ap.parse_args()

    print("=" * 70)
    print("Script 03b - Match preprocessed MRI scans to clinical visits")
    print("=" * 70)
    print(f"  session_map     : {SESSION_MAP}")
    print(f"  master_clinical : {MASTER_CLIN}")
    print(f"  vit_manifest    : {args.vit_manifest}")
    print(f"  label_policy    : {args.label_policy}")
    print(f"  cohort_filter   : {'ALL subjects' if args.no_cohort_filter else 'SNP+MRI only (n=616)'}")
    print(f"  master out dir  : {MATCHED_DIR}")
    print()
    MATCHED_DIR.mkdir(parents=True, exist_ok=True)

    # ---- Load: SPINE = mri_visit_dates_snp_subjects.csv (wide -> long) ------
    md = pd.read_csv(MRI_DATES_CSV)
    visit_cols = [c for c in md.columns if c != "participant_id"]
    md_long = md.melt(id_vars="participant_id", value_vars=visit_cols,
                      var_name="bids_ses", value_name="mri_date_str")
    # Drop empty cells ("n/a" or NaN)
    md_long = md_long[md_long["mri_date_str"].notna()
                      & (md_long["mri_date_str"].astype(str).str.lower() != "n/a")]
    md_long["mri_date"] = pd.to_datetime(md_long["mri_date_str"], errors="coerce")
    md_long = md_long.dropna(subset=["mri_date"])
    md_long["bids_sub"] = md_long["participant_id"].apply(normalise_subject)
    sm = md_long[["bids_sub", "bids_ses", "mri_date"]].drop_duplicates(
        ["bids_sub", "bids_ses"]).reset_index(drop=True)
    sm["scan_row_id"] = sm.index   # unique per (subject, visit)
    print(f"  Spine: mri_visit_dates_snp_subjects.csv")
    print(f"    subjects                : {md['participant_id'].nunique()}")
    print(f"    (subject, visit) scans  : {len(sm):,}")

    clin = pd.read_csv(MASTER_CLIN, low_memory=False)
    clin["bids_sub"]      = clin["Patient_ID"].apply(normalise_subject)
    clin["clinical_date"] = pd.to_datetime(clin["Date"], errors="coerce")
    clin = clin.dropna(subset=["clinical_date"])
    # Dedupe at (bids_sub, VISCODE_long) — 47 dup pairs in ADNI from visit retakes
    n_pre = len(clin)
    clin = clin.sort_values("clinical_date").drop_duplicates(
        subset=["bids_sub", "VISCODE_long"], keep="first")
    print(f"  master_clinical rows     : {n_pre:,} -> {len(clin):,} after dedupe (subject,VISCODE)")

    # Optional: attach the QC-selected ImageUID for each (sub, ses)
    image_uid_by_key = {}
    if SCAN_SELECTION.is_file():
        ss = pd.read_csv(SCAN_SELECTION, dtype=str)
        ss = ss[["bids_sub", "bids_ses", "ImageUID_selected"]].dropna(
            subset=["bids_sub", "bids_ses"])
        image_uid_by_key = {(r["bids_sub"], r["bids_ses"]): r["ImageUID_selected"]
                            for _, r in ss.iterrows()}
        sm["ImageUID"] = sm.apply(
            lambda r: image_uid_by_key.get((r["bids_sub"], r["bids_ses"]), ""), axis=1)
        n_ids = (sm["ImageUID"] != "").sum()
        print(f"  ImageUID attached (scan_selection.csv): {n_ids:,}/{len(sm):,}")

    if not args.no_cohort_filter:
        cohort = pd.read_csv(SNP_TSV, sep="\t")
        cohort_ids = set(cohort["participant_id"].apply(normalise_subject))
        before = len(sm)
        sm = sm[sm["bids_sub"].isin(cohort_ids)].reset_index(drop=True)
        sm["scan_row_id"] = sm.index
        clin = clin[clin["bids_sub"].isin(cohort_ids)]
        print(f"  Cohort filter (SNP+MRI)  : {before:,} -> {len(sm):,} MRI scans, "
              f"{len(clin):,} clinical rows")
    print()

    # ---- Strategy A: direct VISCODE join ------------------------------------
    print("[Strategy A] Direct VISCODE join (ADNI's authoritative pairing)")
    a = sm.merge(
        clin[["bids_sub", "VISCODE_long", "clinical_date"] + LABEL_COLS],
        left_on=["bids_sub", "bids_ses"],
        right_on=["bids_sub", "VISCODE_long"],
        how="left",
    )
    a["days_diff_A"] = (a["mri_date"] - a["clinical_date"]).abs().dt.days
    a["matched_A"]   = a["clinical_date"].notna()
    n_a = int(a["matched_A"].sum())
    print(f"  matched: {n_a:,}/{len(a):,} ({100 * n_a / len(a):.1f}%)")
    if n_a:
        d = a.loc[a["matched_A"], "days_diff_A"]
        print(f"  |days_diff|  median={d.median():.0f}  p95={d.quantile(0.95):.0f}  "
              f"max={d.max():.0f}  n>14d={(d > 14).sum()}")
    print()

    # ---- Strategy B: nearest-by-date within +-14d ---------------------------
    print("[Strategy B] Nearest clinical visit by date (+-14d)")
    sm_b = sm.sort_values("mri_date").reset_index(drop=True)
    cl_b = clin[["bids_sub", "VISCODE_long", "clinical_date"] + LABEL_COLS].copy()
    cl_b = cl_b.sort_values("clinical_date").reset_index(drop=True)
    b14 = pd.merge_asof(
        sm_b, cl_b,
        left_on="mri_date", right_on="clinical_date",
        by="bids_sub",
        direction="nearest",
        tolerance=pd.Timedelta(days=14),
    )
    b14["matched_B14"] = b14["clinical_date"].notna()
    n_b14 = int(b14["matched_B14"].sum())
    print(f"  matched +-14d: {n_b14:,}/{len(b14):,} ({100 * n_b14 / len(b14):.1f}%)")
    print()

    # ---- Strategy A vs B comparison (per-scan via scan_row_id) -------------
    cmp = pd.DataFrame({
        "scan_row_id": sm["scan_row_id"].values,
        "matched_A":   a["matched_A"].values,
    }).merge(
        b14[["scan_row_id", "matched_B14"]], on="scan_row_id", how="left"
    )
    cmp["matched_B14"] = cmp["matched_B14"].fillna(False)

    n_total   = len(cmp)
    n_both    = int(((cmp["matched_A"]) & (cmp["matched_B14"])).sum())
    n_only_a  = int(((cmp["matched_A"]) & (~cmp["matched_B14"])).sum())
    n_only_b  = int(((~cmp["matched_A"]) & (cmp["matched_B14"])).sum())
    n_neither = int(((~cmp["matched_A"]) & (~cmp["matched_B14"])).sum())
    n_union   = n_both + n_only_a + n_only_b
    print("[A vs B] (per-scan)")
    print(f"  matched by both A AND B (+-14d): {n_both:,}")
    print(f"  matched only by A (VISCODE)    : {n_only_a:,}  <- ADNI label exists but date >14d off")
    print(f"  matched only by B (+-14d)      : {n_only_b:,}  <- date close but ADNI VISCODE differs")
    print(f"  matched by A OR B (union)      : {n_union:,}/{n_total:,} ({100 * n_union / n_total:.1f}%)")
    print(f"  matched by neither             : {n_neither:,}")
    print()

    # ---- Extended-tolerance recovery (audit only) ---------------------------
    print("[Recovery] Strategy B at wider tolerances")
    recovery = {}
    for tol_days in EXTENDED_TOLERANCES_DAYS:
        b_tol = pd.merge_asof(
            sm_b, cl_b,
            left_on="mri_date", right_on="clinical_date",
            by="bids_sub",
            direction="nearest",
            tolerance=pd.Timedelta(days=tol_days),
        )
        n = int(b_tol["clinical_date"].notna().sum())
        recovery[tol_days] = n
        print(f"  +-{tol_days:>3d}d: {n:,} matched ({100 * n / len(sm_b):.1f}%)")
    print()

    # ---- Per-VISCODE audit table --------------------------------------------
    # Attach per-scan B+-14d match flag to a, indexed by scan_row_id
    a_aug = a.merge(cmp[["scan_row_id", "matched_B14"]], on="scan_row_id", how="left")
    audit_rows = []
    for vis, group in a_aug.groupby("bids_ses"):
        n_total      = len(group)
        n_match_a    = int(group["matched_A"].sum())
        days = group.loc[group["matched_A"], "days_diff_A"]
        n_match_b14  = int(group["matched_B14"].sum())

        audit_rows.append({
            "viscode_adni":           vis,
            "n_scans":                n_total,
            "n_matched_viscode":      n_match_a,
            "pct_matched_viscode":    round(100 * n_match_a / n_total, 1) if n_total else 0,
            "n_matched_nearest_14d":  n_match_b14,
            "days_diff_median":       round(days.median(), 1) if len(days) else np.nan,
            "days_diff_p95":          round(days.quantile(0.95), 1) if len(days) else np.nan,
            "days_diff_max":          int(days.max()) if len(days) else 0,
            "n_days_diff_gt_14":      int((days > 14).sum()) if len(days) else 0,
            "n_days_diff_gt_30":      int((days > 30).sum()) if len(days) else 0,
        })
    audit_df = pd.DataFrame(audit_rows)
    audit_df["_sort"] = audit_df["viscode_adni"].apply(viscode_months)
    audit_df = audit_df.sort_values("_sort").drop(columns="_sort")
    audit_df.to_csv(OUT_AUDIT, index=False)
    print(f"Saved audit table -> {OUT_AUDIT}")

    # ---- Strategy comparison summary ----------------------------------------
    total = len(cmp)
    cmp_metrics = [
        ("total_scans",                  total,    "All MRI scans in cohort"),
        ("A_matched_viscode",            n_a,      "Strategy A: direct ADNI VISCODE2 join"),
        ("B_matched_within_14d",         n_b14,    "Strategy B: nearest clinical visit within +-14d"),
        ("A_AND_B_intersection",         n_both,   "Both strategies agree (clean match)"),
        ("A_only_viscode_only",          n_only_a, "ADNI VISCODE matches, but date >14d off"),
        ("B_only_nearest14d_only",       n_only_b, "Date close (+-14d), but ADNI VISCODE differs"),
        ("A_OR_B_union",                 n_union,  "Matched by at least one strategy (use this for label coverage)"),
        ("neither",                      n_neither,"Unmatched by both strategies"),
    ]
    rows = [{"metric": m, "count": c, "pct_of_total": round(100 * c / total, 1),
             "description": d} for (m, c, d) in cmp_metrics]
    for t in EXTENDED_TOLERANCES_DAYS:
        rows.append({
            "metric":       f"B_recovery_within_{t}d",
            "count":        recovery[t],
            "pct_of_total": round(100 * recovery[t] / total, 1),
            "description":  f"Strategy B at wider +-{t}d (audit only; not used for labels)",
        })
    pd.DataFrame(rows).to_csv(OUT_COMPARE, index=False)
    print(f"Saved strategy comparison -> {OUT_COMPARE}")

    # ---- Histogram per VISCODE ----------------------------------------------
    common_vis = ["bl", "m06", "m12", "m18", "m24", "m36", "m48", "m60", "m72", "m84", "m96"]
    panels = [v for v in common_vis if v in a["bids_ses"].unique()]
    n_pan = len(panels)
    if n_pan:
        ncols = 4
        nrows = int(np.ceil(n_pan / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.2 * nrows),
                                 sharex=False)
        axes = np.atleast_2d(axes).flatten()
        for i, vis in enumerate(panels):
            ax = axes[i]
            d = a.loc[(a["bids_ses"] == vis) & a["matched_A"], "days_diff_A"]
            if not len(d):
                ax.set_title(f"ses-{vis} (no match)", fontsize=9)
                ax.axis("off")
                continue
            cap = min(int(d.quantile(0.99)) + 1, 365)
            ax.hist(d.clip(upper=cap), bins=30, color="#4878CF", alpha=0.85)
            ax.axvline(14, color="red", linestyle="--", lw=0.8,
                       label="ADNI +-14d protocol")
            ax.set_title(f"ses-{vis}   n={len(d)}   median={d.median():.0f}d",
                         fontsize=9)
            ax.set_xlabel("|MRI date - clinical date| (days)", fontsize=8)
            ax.set_ylabel("scans", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=7, loc="upper right")
        for j in range(n_pan, len(axes)):
            axes[j].axis("off")
        fig.suptitle("Direct-VISCODE join: |MRI date - clinical date| per ADNI session\n"
                     "(median = median absolute days between paired MRI and clinical visit)",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(OUT_HIST, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved histogram -> {OUT_HIST}")
    print()

    # ---- Build MASTER labels CSV (every cohort MRI scan) --------------------
    print("[Master] Building master_mri_clinical_matched_labels.csv")
    print(f"  Scope: every MRI scan in the {'SNP+MRI cohort' if not args.no_cohort_filter else 'full ADNI cohort'} ({len(sm):,} scans)")

    # Strategy A is already in `a` — it's the per-scan VISCODE2 join.
    master = a.copy()
    master["match_status"] = np.where(
        master["matched_A"], "viscode_exact", "unmatched")

    # Strategy B fallback: nearest clinical visit within +-14d (protocol window)
    if args.label_policy == "viscode_then_nearest_14d":
        unmatched_mask = master["match_status"] == "unmatched"
        n_unmatched = int(unmatched_mask.sum())
        if n_unmatched:
            sub_unmatched = (master.loc[unmatched_mask, ["scan_row_id", "bids_sub", "mri_date"]]
                                  .sort_values("mri_date"))
            cl_sorted = clin[["bids_sub", "VISCODE_long", "clinical_date"] + LABEL_COLS].sort_values("clinical_date")
            recover = pd.merge_asof(
                sub_unmatched, cl_sorted,
                left_on="mri_date", right_on="clinical_date",
                by="bids_sub", direction="nearest",
                tolerance=pd.Timedelta(days=14),
            )
            recovered_b = recover[recover["clinical_date"].notna()]
            print(f"  Strategy B (+-14d nearest) recovered: {len(recovered_b)}/{n_unmatched}")
            if not recovered_b.empty:
                rec = recovered_b.set_index("scan_row_id")
                sel = master["scan_row_id"].isin(rec.index)
                for col in ["VISCODE_long", "clinical_date"] + LABEL_COLS:
                    master.loc[sel, col] = master.loc[sel, "scan_row_id"].map(rec[col]).values
                master.loc[sel, "match_status"] = "nearest_within_14d"

    # Recompute days_diff (now reflects whichever strategy matched)
    master["days_diff"] = (master["mri_date"] - master["clinical_date"]).abs().dt.days

    # ---- Audit-only: baseline-offset consistency check ---------------------
    # For each scan, compute mri_date - subject's bl clinical date, then
    # check whether it lands within +-14d of the canonical ADNI interval
    # for the MRI's ADNI VISCODE2 (m12=365d, m24=730d, ...). Does NOT
    # change match_status; just adds two diagnostic columns.
    bl_dates = (clin[clin["VISCODE_long"] == "bl"][["bids_sub", "clinical_date"]]
                    .rename(columns={"clinical_date": "bl_clinical_date"}))
    master = master.merge(bl_dates, on="bids_sub", how="left")
    master["bl_offset_days"] = (master["mri_date"] - master["bl_clinical_date"]).dt.days
    master["bl_offset_expected_days"] = master["bids_ses"].apply(
        lambda s: viscode_months(s) * 30.4375 if viscode_months(s) < 9000 else np.nan)
    master["bl_offset_diff"] = (master["bl_offset_days"]
                                 - master["bl_offset_expected_days"]).abs()
    master["bl_offset_consistent"] = master["bl_offset_diff"] <= 14

    n_consist = int(master["bl_offset_consistent"].fillna(False).sum())
    n_have_bl = int(master["bl_clinical_date"].notna().sum())
    print(f"  Baseline-offset audit: {n_consist:,}/{n_have_bl:,} scans have ADNI VISCODE "
          f"consistent with date (+-14d of canonical interval)")
    unmatched_mask = master["match_status"] == "unmatched"
    if unmatched_mask.any():
        n_unmatched_consist = int(
            master.loc[unmatched_mask, "bl_offset_consistent"].fillna(False).sum())
        n_unmatched_have_bl = int(
            master.loc[unmatched_mask, "bl_clinical_date"].notna().sum())
        print(f"  Among {int(unmatched_mask.sum())} unmatched scans: "
              f"{n_unmatched_consist}/{n_unmatched_have_bl} VISCODE consistent with "
              f"bl-offset (i.e. ADNI's session label looks correct, just no clinical row "
              f"at that VISCODE)")

    master_cols = (
        ["bids_sub", "bids_ses", "mri_date", "ImageUID",
         "clinical_date", "days_diff", "match_status", "VISCODE_long",
         "bl_clinical_date", "bl_offset_days", "bl_offset_expected_days",
         "bl_offset_diff", "bl_offset_consistent"]
        + LABEL_COLS
    )
    master_cols = [c for c in master_cols if c in master.columns]
    master["_sort_ses"] = master["bids_ses"].apply(viscode_months)
    master = master.sort_values(["bids_sub", "_sort_ses", "mri_date"]).drop(columns="_sort_ses")
    master[master_cols].to_csv(OUT_MASTER, index=False)
    print(f"  Saved master CSV -> {OUT_MASTER}")
    master_status = master["match_status"].value_counts().to_dict()
    for k, v in master_status.items():
        print(f"    {k:<25s} {v:>6d}")
    print()

    # ---- Per-session summary table (CSV + PNG/PDF) --------------------------
    print("[Summary] Building per-session match-status summary")
    rows = []
    for ses, g in master.groupby("bids_ses"):
        n = len(g)
        n_a    = int((g["match_status"] == "viscode_exact").sum())
        n_b    = int((g["match_status"] == "nearest_within_14d").sum())
        n_um   = g["match_status"] == "unmatched"
        n_um_c = int((n_um & g["bl_offset_consistent"].fillna(False)).sum())
        n_um_i = int(n_um.sum() - n_um_c)
        rows.append({
            "ses": ses, "n_scans": n,
            "n_matched_viscode": n_a,
            "n_added_nearest_14d": n_b,
            "n_unmatched_consistent": n_um_c,
            "n_unmatched_inconsistent": n_um_i,
        })
    # Add a TOTAL row
    summary_df = pd.DataFrame(rows)
    summary_df["_sort"] = summary_df["ses"].apply(viscode_months)
    summary_df = summary_df.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    totals = {
        "ses": "TOTAL",
        "n_scans": int(summary_df["n_scans"].sum()),
        "n_matched_viscode": int(summary_df["n_matched_viscode"].sum()),
        "n_added_nearest_14d": int(summary_df["n_added_nearest_14d"].sum()),
        "n_unmatched_consistent": int(summary_df["n_unmatched_consistent"].sum()),
        "n_unmatched_inconsistent": int(summary_df["n_unmatched_inconsistent"].sum()),
    }
    summary_df = pd.concat([summary_df, pd.DataFrame([totals])], ignore_index=True)
    render_summary_table(summary_df, OUT_SUMMARY_CSV, OUT_SUMMARY_PNG, OUT_SUMMARY_PDF)
    print(f"  Saved summary -> {OUT_SUMMARY_CSV}")
    print(f"                -> {OUT_SUMMARY_PNG}")
    print(f"                -> {OUT_SUMMARY_PDF}")
    print()

    # ---- Clinical-aligned variant ------------------------------------------
    # For matched scans (Strategy A or B), replace bids_ses with VISCODE_long
    # (the matched clinical row's session). Preserves the original MRI tag in
    # bids_ses_mri. This tabulation includes intermediate clinical visits
    # (m06, m18, m30, ...) that have no MRI-side spine column. The user can
    # decide whether to prefer this view or the ADNI-MRI-VISCODE view.
    print("[Aligned] Building clinical-aligned variant under "
          f"{ALIGNED_DIR.name}/")
    ALIGNED_DIR.mkdir(parents=True, exist_ok=True)
    aligned = master.copy()
    aligned["bids_ses_mri"] = aligned["bids_ses"]
    matched_mask = aligned["match_status"].isin(MATCHED_STATUSES)
    aligned.loc[matched_mask, "bids_ses"] = aligned.loc[matched_mask, "VISCODE_long"]
    aligned["_sort_ses"] = aligned["bids_ses"].apply(viscode_months)
    aligned = aligned.sort_values(["bids_sub", "_sort_ses", "mri_date"]).drop(columns="_sort_ses")
    aligned_cols = master_cols + ["bids_ses_mri"]
    aligned_cols = [c for c in aligned_cols if c in aligned.columns]
    aligned[aligned_cols].to_csv(OUT_ALIGNED_MASTER, index=False)
    print(f"  Saved aligned master -> {OUT_ALIGNED_MASTER}")

    # Per-session summary for the aligned view
    rows_a = []
    for ses, g in aligned.groupby("bids_ses"):
        n = len(g)
        n_a    = int((g["match_status"] == "viscode_exact").sum())
        n_b    = int((g["match_status"] == "nearest_within_14d").sum())
        n_um   = g["match_status"] == "unmatched"
        n_um_c = int((n_um & g["bl_offset_consistent"].fillna(False)).sum())
        n_um_i = int(n_um.sum() - n_um_c)
        rows_a.append({
            "ses": ses, "n_scans": n,
            "n_matched_viscode": n_a,
            "n_added_nearest_14d": n_b,
            "n_unmatched_consistent": n_um_c,
            "n_unmatched_inconsistent": n_um_i,
        })
    aligned_summary_df = pd.DataFrame(rows_a)
    aligned_summary_df["_sort"] = aligned_summary_df["ses"].apply(viscode_months)
    aligned_summary_df = aligned_summary_df.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    totals_a = {
        "ses": "TOTAL",
        "n_scans":                  int(aligned_summary_df["n_scans"].sum()),
        "n_matched_viscode":        int(aligned_summary_df["n_matched_viscode"].sum()),
        "n_added_nearest_14d":      int(aligned_summary_df["n_added_nearest_14d"].sum()),
        "n_unmatched_consistent":   int(aligned_summary_df["n_unmatched_consistent"].sum()),
        "n_unmatched_inconsistent": int(aligned_summary_df["n_unmatched_inconsistent"].sum()),
    }
    aligned_summary_df = pd.concat([aligned_summary_df, pd.DataFrame([totals_a])],
                                    ignore_index=True)
    render_summary_table(aligned_summary_df, OUT_ALIGNED_SUMMARY_CSV,
                         OUT_ALIGNED_SUMMARY_PNG, OUT_ALIGNED_SUMMARY_PDF)
    print(f"  Saved aligned summary -> {OUT_ALIGNED_SUMMARY_CSV}")
    print(f"                       -> {OUT_ALIGNED_SUMMARY_PNG}")
    print(f"                       -> {OUT_ALIGNED_SUMMARY_PDF}")
    print()

    # ---- Derive preprocessed-subset view from MASTER ------------------------
    print("[Subset] Filtering master to preprocessed scans in vit_manifest.csv")
    if not Path(args.vit_manifest).is_file():
        print(f"  vit_manifest.csv not found at {args.vit_manifest}")
        print(f"  Skipping preprocessed-subset CSV. Master + audit are still produced.")
        return

    manifest = pd.read_csv(args.vit_manifest, dtype=str)
    manifest = manifest[manifest["status"].isin(["ok", "skipped_exists"])].copy()
    manifest["bids_sub"] = manifest["bids_sub"].astype(str)
    print(f"  preprocessed scans in manifest: {len(manifest)}")

    # 'first scan per subject' (mni_staging local layout): all rows tagged
    # ses-bl. Resolve to each subject's earliest scan in master, then keep one.
    if (manifest["bids_ses"] == "bl").all() and len(manifest):
        print("  Manifest tagged all-ses-bl ('first scan per subject' layout).")
        print("  Picking each subject's earliest master row as the preprocessed scan.")
        master_subset = (master.sort_values(["bids_sub", "mri_date"])
                              .drop_duplicates("bids_sub", keep="first"))
        local = manifest[["bids_sub", "output_path"]].merge(
            master_subset[master_cols], on="bids_sub", how="left")
    else:
        # Sessionwise layout (HPC): match master on (bids_sub, bids_ses).
        # If multiple master rows match, pick the earliest mri_date.
        master_subset = (master.sort_values(["bids_sub", "bids_ses", "mri_date"])
                              .drop_duplicates(["bids_sub", "bids_ses"], keep="first"))
        local = manifest[["bids_sub", "bids_ses", "output_path"]].merge(
            master_subset[master_cols], on=["bids_sub", "bids_ses"], how="left")

    local["_sort_ses"] = local["bids_ses"].apply(viscode_months)
    local = local.sort_values(["bids_sub", "_sort_ses", "mri_date"]).drop(columns="_sort_ses")
    local_cols = master_cols + ["output_path"]
    local[local_cols].to_csv(OUT_LABELS, index=False)
    print(f"  Saved subset CSV -> {OUT_LABELS}")
    n_with_labels = local["match_status"].notna().sum()
    print(f"  Of {len(local)} preprocessed scans: {n_with_labels} have a master row, "
          f"{len(local) - n_with_labels} not in cohort filter")
    print()

    # ---- Summary ------------------------------------------------------------
    print("=" * 70)
    print("Outputs:")
    print(f"  Master CSV    : {OUT_MASTER}")
    print(f"  Audit CSV     : {OUT_AUDIT}")
    print(f"  Comparison    : {OUT_COMPARE}")
    print(f"  Histogram     : {OUT_HIST}")
    print(f"  Subset CSV    : {OUT_LABELS}")
    print("=" * 70)


if __name__ == "__main__":
    main()
