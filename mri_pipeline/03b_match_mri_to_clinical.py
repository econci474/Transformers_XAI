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
by `clinical_pipeline/01b_build_clinical_csv.py`.

Inputs
------
- `metadata/session_map.csv`      MRI inventory with ADNI-assigned VISCODE2
                                   per scan (from bidsification/01_build_*.py).
- `bids/imaging/mri_visit_dates_snp_subjects.csv`
                                   Wide MRI dates per subject, one column per
                                   ADNI session.
- `derivatives/clinical/verbose/longitudinal/master_clinical_verbose.csv`
                                   Long-format clinical visits + label cols.
- `bids/genotype/subjects_with_snp_and_mri.tsv`
                                   616-subject SNP+MRI cohort filter.
- `derivatives/vit_inputs/vit_manifest.csv`
                                   What's actually preprocessed and on disk.

Outputs (under `derivatives/vit_inputs/`)
-----------------------------------------
- `mri_clinical_matched_labels.csv`     One row per preprocessed scan with:
      bids_sub, bids_ses_adni, mri_date, matched_viscode, matched_clinical_date,
      days_diff, match_status, all clinical Label_* columns, output_path.

- `mri_clinical_match_audit.csv`        Per-VISCODE summary of |days_diff|
      median/p95/max, n_matched_viscode_exact, n_matched_within_14d,
      n_matched_within_extended, n_unmatched, recovery rates at multiple
      tolerances.

- `mri_clinical_date_diffs_per_viscode.png`
      Histogram per ADNI session of |MRI date - clinical date| for the
      direct-VISCODE join. Reveals how tight ADNI's per-VISCODE scheduling
      really is (the protocol says +-14 days, but reality may differ).

Match strategies (and why both are reported)
--------------------------------------------
A) **Direct VISCODE join** (default labels source)
   Join `(bids_sub, bids_ses_adni)` from MRI side to
   `(Patient_ID, VISCODE_long)` from clinical side. ADNI assigned both labels
   from the same scheduling protocol, so this is the authoritative pairing.
   Pros: zero ambiguity, no leakage. Cons: if a clinical row at that VISCODE
   doesn't exist for the subject (missed assessment), we get no match.

B) **Nearest-by-date** (recovery only)
   pd.merge_asof(direction="nearest") within +-14d. Used as a fallback for
   scans whose ADNI-VISCODE has no matching clinical row. Tolerance is from
   the ADNI-2/3 protocol "all assessments within 2 weeks of target date".

C) **Extended-tolerance recovery** (audit only)
   For scans unmatched at +-14d, report whether they would match at +-30/60/
   90/180 days. The script writes the counts but does NOT use these as labels
   (would risk leakage / mis-attribution).

Final label assigned to each scan
---------------------------------
match_status priority: viscode_exact > nearest_within_14d > unmatched
The `--label-policy` flag controls whether to fall back from A->B at all.
Default = "viscode_only" (strict ADNI labels), use --label-policy=viscode_then_nearest_14d
to enable the +-14d fallback.

Usage
-----
    python mri_pipeline/03b_match_mri_to_clinical.py
    python mri_pipeline/03b_match_mri_to_clinical.py --label-policy viscode_then_nearest_14d
    python mri_pipeline/03b_match_mri_to_clinical.py --vit-manifest /alt/path/vit_manifest.csv
"""

import os
import sys
import argparse
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---- Paths -----------------------------------------------------------------
PROJECT_ROOT  = Path(r"D:\ADNI_BIDS_project")
SESSION_MAP   = PROJECT_ROOT / "metadata" / "session_map.csv"
MRI_DATES_CSV = PROJECT_ROOT / "bids" / "imaging" / "mri_visit_dates_snp_subjects.csv"
MASTER_CLIN   = PROJECT_ROOT / "derivatives" / "clinical" / "verbose" / "longitudinal" / "master_clinical_verbose.csv"
SNP_TSV       = PROJECT_ROOT / "bids" / "genotype" / "subjects_with_snp_and_mri.tsv"
VIT_INPUTS    = PROJECT_ROOT / "derivatives" / "vit_inputs"
VIT_MANIFEST  = VIT_INPUTS / "vit_manifest.csv"

OUT_LABELS    = VIT_INPUTS / "mri_clinical_matched_labels.csv"
OUT_AUDIT     = VIT_INPUTS / "mri_clinical_match_audit.csv"
OUT_HIST      = VIT_INPUTS / "mri_clinical_date_diffs_per_viscode.png"
OUT_COMPARE   = VIT_INPUTS / "mri_clinical_strategy_comparison.csv"

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
                    default="viscode_only",
                    help="viscode_only: only assign labels via ADNI VISCODE2 join (strictest, default). "
                         "viscode_then_nearest_14d: fall back to nearest clinical visit within +-14d "
                         "for scans whose ADNI VISCODE has no matching clinical row.")
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
    print()

    # ---- Load ---------------------------------------------------------------
    sm = pd.read_csv(SESSION_MAP, low_memory=False)
    sm = sm[["SubjectID", "bids_sub", "bids_ses", "adni_viscode", "StudyDate"]].copy()
    sm = sm.dropna(subset=["StudyDate"])
    sm["mri_date"] = pd.to_datetime(sm["StudyDate"], errors="coerce")
    sm = sm.dropna(subset=["mri_date"])
    print(f"  session_map rows         : {len(sm):,}")

    clin = pd.read_csv(MASTER_CLIN, low_memory=False)
    clin["bids_sub"]      = clin["Patient_ID"].apply(normalise_subject)
    clin["clinical_date"] = pd.to_datetime(clin["Date"], errors="coerce")
    clin = clin.dropna(subset=["clinical_date"])
    print(f"  master_clinical rows     : {len(clin):,}")

    if not args.no_cohort_filter:
        cohort = pd.read_csv(SNP_TSV, sep="\t")
        cohort_ids = set(cohort["participant_id"].apply(normalise_subject))
        before = len(sm)
        sm = sm[sm["bids_sub"].isin(cohort_ids)]
        print(f"  Cohort filter (SNP+MRI)  : {before:,} -> {len(sm):,} MRI scans")
        clin = clin[clin["bids_sub"].isin(cohort_ids)]
        print(f"  Cohort filter (clinical) : {len(clin):,} clinical rows in cohort")
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

    # ---- Strategy A vs B comparison -----------------------------------------
    a_keyed   = a.set_index(["bids_sub", "bids_ses"]).index
    b14_keyed = b14.set_index(["bids_sub", "bids_ses"])
    matched_in_a   = set(a.loc[a["matched_A"]].set_index(["bids_sub", "bids_ses"]).index)
    matched_in_b14 = set(b14.loc[b14["matched_B14"]].set_index(["bids_sub", "bids_ses"]).index)

    only_a   = matched_in_a - matched_in_b14
    only_b14 = matched_in_b14 - matched_in_a
    both     = matched_in_a & matched_in_b14
    neither  = set(a_keyed) - matched_in_a - matched_in_b14
    print("[A vs B]")
    print(f"  matched by both (A & B+-14d) : {len(both):,}")
    print(f"  matched only by A (VISCODE)  : {len(only_a):,}  <- ADNI label exists but date >14d off")
    print(f"  matched only by B (+-14d)    : {len(only_b14):,}  <- date close but ADNI VISCODE differs")
    print(f"  matched by neither           : {len(neither):,}")
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
    audit_rows = []
    for vis, group in a.groupby("bids_ses"):
        n_total      = len(group)
        n_match_a    = int(group["matched_A"].sum())
        days = group.loc[group["matched_A"], "days_diff_A"]

        # how many of this VISCODE bucket also matched in B+-14d
        keys = group.set_index(["bids_sub", "bids_ses"]).index
        n_match_b14 = sum(1 for k in keys if k in matched_in_b14)

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
    audit_df = pd.DataFrame(audit_rows).sort_values("viscode_adni")
    audit_df.to_csv(OUT_AUDIT, index=False)
    print(f"Saved audit table -> {OUT_AUDIT}")

    # ---- Strategy comparison summary ----------------------------------------
    pd.DataFrame({
        "metric": [
            "total_scans", "matched_A_viscode", "matched_B_within_14d",
            "matched_both", "matched_only_A_viscode", "matched_only_B_14d",
            "matched_neither",
        ] + [f"recovery_within_{t}d" for t in EXTENDED_TOLERANCES_DAYS],
        "count": [
            len(sm), n_a, n_b14, len(both), len(only_a), len(only_b14), len(neither),
        ] + [recovery[t] for t in EXTENDED_TOLERANCES_DAYS],
    }).to_csv(OUT_COMPARE, index=False)
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
            ax.axvline(14, color="red", linestyle="--", lw=0.8, label="14d")
            ax.set_title(f"ses-{vis}  n={len(d)}  med={d.median():.0f}d", fontsize=9)
            ax.set_xlabel("|MRI - clinical| (days)", fontsize=8)
            ax.set_ylabel("scans", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.legend(fontsize=7, loc="upper right")
        for j in range(n_pan, len(axes)):
            axes[j].axis("off")
        fig.suptitle("Direct-VISCODE join: |MRI date - clinical date| per ADNI session",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(OUT_HIST, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved histogram -> {OUT_HIST}")
    print()

    # ---- Build per-scan labels CSV (for the preprocessed scans) -------------
    print("[Labels] Resolving labels for preprocessed scans in vit_manifest.csv")
    if not Path(args.vit_manifest).is_file():
        print(f"  WARNING: vit_manifest.csv not found at {args.vit_manifest}")
        print(f"  Skipping per-scan labels CSV. Audit outputs above are still valid.")
        return

    manifest = pd.read_csv(args.vit_manifest, dtype=str)
    manifest = manifest[manifest["status"].isin(["ok", "skipped_exists"])].copy()
    manifest["bids_sub"] = manifest["bids_sub"].astype(str)
    print(f"  preprocessed scans in manifest: {len(manifest)}")

    # The "first scan per subject" mni_staging tags everything as ses-bl, but
    # the ADNI-correct VISCODE is the subject's earliest available session in
    # session_map. Resolve that here.
    if (manifest["bids_ses"] == "bl").all():
        print("  Manifest looks like 'first scan per subject' (all ses-bl).")
        print("  Resolving each scan's true ADNI VISCODE2 from session_map...")
        # earliest scan per subject in session_map (this matches what the
        # 'first' sMRIprep run would have picked up)
        earliest = sm.sort_values(["bids_sub", "mri_date"]).drop_duplicates(
            "bids_sub", keep="first")
        earliest = earliest[["bids_sub", "bids_ses", "adni_viscode", "mri_date"]]
        earliest = earliest.rename(columns={"bids_ses": "bids_ses_resolved"})
        manifest = manifest.merge(earliest, on="bids_sub", how="left")
        manifest["bids_ses"] = manifest["bids_ses_resolved"].fillna(manifest["bids_ses"])
        n_resolved = manifest["bids_ses_resolved"].notna().sum()
        print(f"  resolved {n_resolved}/{len(manifest)} scans")
    else:
        manifest["mri_date"] = pd.NaT
        # sessionwise input: trust the bids_ses already in the manifest
        sm_lookup = sm.set_index(["bids_sub", "bids_ses"])["mri_date"].to_dict()
        manifest["mri_date"] = manifest.apply(
            lambda r: sm_lookup.get((r["bids_sub"], r["bids_ses"])), axis=1)

    # Pull labels via direct VISCODE join (Strategy A). Optionally fall back
    # to nearest-by-date within 14d if --label-policy=viscode_then_nearest_14d.
    label_a = manifest.merge(
        clin[["bids_sub", "VISCODE_long", "clinical_date"] + LABEL_COLS],
        left_on=["bids_sub", "bids_ses"],
        right_on=["bids_sub", "VISCODE_long"],
        how="left",
    )
    label_a["match_status"] = np.where(
        label_a["clinical_date"].notna(), "viscode_exact", "unmatched")
    n_matched_a = int((label_a["match_status"] == "viscode_exact").sum())
    print(f"  Strategy A matched: {n_matched_a}/{len(label_a)}")

    if args.label_policy == "viscode_then_nearest_14d":
        unmatched_mask = label_a["match_status"] == "unmatched"
        if unmatched_mask.any():
            sub_unmatched = label_a.loc[unmatched_mask, ["bids_sub", "bids_ses", "mri_date"]].copy()
            sub_unmatched["mri_date"] = pd.to_datetime(sub_unmatched["mri_date"], errors="coerce")
            sub_unmatched = sub_unmatched.dropna(subset=["mri_date"]).sort_values("mri_date")
            cl_b = clin[["bids_sub", "VISCODE_long", "clinical_date"] + LABEL_COLS].sort_values("clinical_date")
            recover = pd.merge_asof(
                sub_unmatched, cl_b,
                left_on="mri_date", right_on="clinical_date",
                by="bids_sub", direction="nearest",
                tolerance=pd.Timedelta(days=14),
            )
            recovered = recover[recover["clinical_date"].notna()]
            print(f"  Strategy B fallback recovered: {len(recovered)} additional scans within +-14d")
            # Splice recovered labels back in
            for col in ["VISCODE_long", "clinical_date"] + LABEL_COLS:
                label_a.loc[recovered.index, col] = recovered[col].values
            label_a.loc[recovered.index, "match_status"] = "nearest_within_14d"

    label_a["days_diff"] = (
        pd.to_datetime(label_a["mri_date"], errors="coerce")
        - pd.to_datetime(label_a["clinical_date"], errors="coerce")
    ).abs().dt.days

    out_cols = (
        ["bids_sub", "bids_ses", "mri_date", "clinical_date", "days_diff",
         "match_status", "VISCODE_long"] + LABEL_COLS + ["output_path"]
    )
    out_cols = [c for c in out_cols if c in label_a.columns]
    label_a[out_cols].to_csv(OUT_LABELS, index=False)
    print(f"  Saved labels CSV -> {OUT_LABELS}")
    print()

    # ---- Summary ------------------------------------------------------------
    print("=" * 70)
    final_status = label_a["match_status"].value_counts().to_dict()
    print("Per-scan label resolution:")
    for k, v in final_status.items():
        print(f"  {k:<30s} {v:>6d}")
    print("=" * 70)
    print(f"Audit CSV:     {OUT_AUDIT}")
    print(f"Comparison:    {OUT_COMPARE}")
    print(f"Histogram:     {OUT_HIST}")
    print(f"Labels CSV:    {OUT_LABELS}")
    print("=" * 70)


if __name__ == "__main__":
    main()
