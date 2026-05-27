"""
40_compute_aao_labels.py
========================
Build the AAO-augmented converters-only splits + per-subgroup demographic
summary table.

Inputs (per seed, per split):
  D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/
      tabular/baseline/seed_{0,1,2}/{train,val,test}.csv

Filter to converters: AD_bl=1 OR pMCI=1 OR CN_to_AD=1.

Per row, compute:
  age_at_baseline  = fractional years between Date and YOB-07-01 (mid-year)
  years_to_AD      = 0 if AD_bl=1 else raw years_to_AD column
  AAO              = age_at_baseline + years_to_AD            (canonical AAO)
  aao_residual     = AAO - mean(age_at_baseline | train fold) (per-seed)
  gender           = Sex (alias; already 'Male'/'Female' strings)

Outputs (per seed):
  D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion_aao/
      seed_{0,1,2}/{train,val,test}.csv
        — ALL source columns + the 5 new columns above

Summary stats (one-off, single table across all 9 seed × split combos):
  D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion_aao/
      summary/aao_demographic_summary.{tsv,png}

Usage (local Windows, env: snp):
  python snp_pipeline/40_compute_aao_labels.py
  python snp_pipeline/40_compute_aao_labels.py --dry-run
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SRC_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                "no_cdr_stratified_post_exclusion/tabular/baseline")
OUT_ROOT = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                "no_cdr_stratified_post_exclusion_aao")
SEEDS = (0, 1, 2)
SPLITS = ("train", "val", "test")
SUBGROUPS = ("AD_bl", "pMCI", "CN_to_AD")          # converters only


def _compute_aao_cols(df: pd.DataFrame, train_mean_age: float | None) -> pd.DataFrame:
    """Append AAO + derived columns. If train_mean_age is None, return df with
    aao_residual=NaN — caller is expected to compute train mean first then re-run
    for val/test using that train mean."""
    df = df.copy()

    # age_at_baseline: use 1 July of YOB as proxy birthday (only YOB is known).
    yob = pd.to_numeric(df["YOB"], errors="coerce").astype("Int64")
    date = pd.to_datetime(df["Date"], errors="coerce")
    birth = pd.to_datetime(yob.astype("Int64").astype(str) + "-07-01",
                             errors="coerce")
    df["age_at_baseline"] = (date - birth).dt.days / 365.25

    # years_to_AD: 0 for AD_bl, raw value otherwise. NaN if not AD_bl and the raw
    # column is missing — meaning the subgroup filter should have already
    # removed these rows; assertion-style check below.
    yta_raw = pd.to_numeric(df["years_to_AD"], errors="coerce")
    is_ad_bl = (df["AD_bl"].astype(str) == "1")
    df["years_to_AD"] = np.where(is_ad_bl, 0.0, yta_raw)

    # Canonical AAO
    df["AAO"] = df["age_at_baseline"] + df["years_to_AD"]

    # AAO residual (only filled if train_mean_age provided)
    if train_mean_age is not None:
        df["aao_residual"] = df["AAO"] - train_mean_age
    else:
        df["aao_residual"] = np.nan

    # gender alias
    df["gender"] = df["Sex"]
    return df


def _subgroup_label(row) -> str:
    for k in SUBGROUPS:
        if str(row[k]) == "1":
            return k
    return "OTHER"


def _per_subgroup_stats(df: pd.DataFrame) -> dict:
    """Return summary stats for a single (subgroup) cohort."""
    return {
        "n":            int(len(df)),
        "pct_female":   round(100 * (df["gender"] == "Female").mean(), 1) if len(df) else None,
        "age_mean":     round(df["age_at_baseline"].mean(), 1) if len(df) else None,
        "age_std":      round(df["age_at_baseline"].std(),  1) if len(df) else None,
        "aao_mean":     round(df["AAO"].mean(), 1)          if len(df) else None,
        "aao_std":      round(df["AAO"].std(),  1)          if len(df) else None,
        "yta_mean":     round(df["years_to_AD"].mean(), 1)  if len(df) else None,
        "yta_std":      round(df["years_to_AD"].std(),  1)  if len(df) else None,
        "apoe4_mean":   round(pd.to_numeric(df["APOE4_Dosage"], errors="coerce").mean(), 2)
                            if len(df) else None,
        "pct_apoe4_carrier": round(100 * (pd.to_numeric(df["APOE4_Dosage"], errors="coerce") > 0).mean(), 1)
                            if len(df) else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="Print counts + sample rows; do not write CSVs.")
    args = ap.parse_args()

    if not args.dry_run:
        OUT_ROOT.mkdir(parents=True, exist_ok=True)
        (OUT_ROOT / "summary").mkdir(parents=True, exist_ok=True)
        for s in SEEDS:
            (OUT_ROOT / f"seed_{s}").mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for s in SEEDS:
        # Pass 1: read train, compute train_mean_age from train converters only
        train_src = pd.read_csv(SRC_ROOT / f"seed_{s}/train.csv", dtype=str)
        train_conv = train_src[(train_src["AD_bl"] == "1") |
                               (train_src["pMCI"]  == "1") |
                               (train_src["CN_to_AD"] == "1")].copy()
        train_conv = _compute_aao_cols(train_conv, train_mean_age=None)
        train_mean_age = float(train_conv["age_at_baseline"].mean())
        print(f"\nseed {s}: train converters n={len(train_conv)}  "
              f"mean(age_at_baseline)={train_mean_age:.2f}")

        # Pass 2: train/val/test, applying train_mean_age for aao_residual
        for sp in SPLITS:
            src = pd.read_csv(SRC_ROOT / f"seed_{s}/{sp}.csv", dtype=str)
            mask = (src["AD_bl"] == "1") | (src["pMCI"] == "1") | (src["CN_to_AD"] == "1")
            conv = src[mask].copy()
            conv = _compute_aao_cols(conv, train_mean_age=train_mean_age)
            # Sanity: AAO must be non-null for every converter
            n_bad = conv["AAO"].isna().sum()
            if n_bad:
                print(f"  ** seed{s} {sp}: {n_bad} converters have NaN AAO "
                      f"(unexpected — investigate)")

            # Emit per-subgroup summary
            for sg in SUBGROUPS:
                sub = conv[conv[sg] == "1"]
                stats = _per_subgroup_stats(sub)
                summary_rows.append({"seed": s, "split": sp, "subgroup": sg, **stats})
            # ALL converters row
            stats_all = _per_subgroup_stats(conv)
            summary_rows.append({"seed": s, "split": sp, "subgroup": "ALL_converters", **stats_all})

            out_p = OUT_ROOT / f"seed_{s}/{sp}.csv"
            if not args.dry_run:
                conv.to_csv(out_p, index=False)
                print(f"  wrote {out_p}  ({len(conv)} converters)")
            else:
                print(f"  [DRY] would write {out_p}  ({len(conv)} converters)")

    # ── Summary table ────────────────────────────────────────────────────
    df_sum = pd.DataFrame(summary_rows)
    df_sum = df_sum[["seed", "split", "subgroup", "n", "pct_female",
                       "age_mean", "age_std", "aao_mean", "aao_std",
                       "yta_mean", "yta_std", "apoe4_mean", "pct_apoe4_carrier"]]
    tsv_p = OUT_ROOT / "summary/aao_demographic_summary.tsv"
    if not args.dry_run:
        df_sum.to_csv(tsv_p, sep="\t", index=False)
        print(f"\nwrote {tsv_p}")

    # PNG render
    fig_h = max(2.0, 0.28 * len(df_sum) + 1.4)
    fig, ax = plt.subplots(figsize=(13.5, fig_h))
    ax.axis("off")
    ax.set_title("AAO converters-only cohort — demographics per "
                 "(seed × split × subgroup)",
                 fontsize=12, fontweight="bold", pad=14, loc="left")
    display = df_sum.copy()
    for c in ("age_mean", "age_std", "aao_mean", "aao_std",
              "yta_mean", "yta_std", "apoe4_mean"):
        display[c] = display[c].apply(lambda v: ("" if pd.isna(v) else f"{v:.1f}"))
    display["pct_female"] = display["pct_female"].apply(
        lambda v: ("" if pd.isna(v) else f"{v:.1f}%"))
    display["pct_apoe4_carrier"] = display["pct_apoe4_carrier"].apply(
        lambda v: ("" if pd.isna(v) else f"{v:.1f}%"))

    tbl = ax.table(cellText=display.values.tolist(),
                   colLabels=display.columns.tolist(),
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.1)
    for j in range(len(display.columns)):
        c = tbl[(0, j)]; c.set_facecolor("#222"); c.set_text_props(color="white", weight="bold")
    # Highlight ALL_converters rows
    for i, row in enumerate(display.itertuples(index=False), start=1):
        if row.subgroup == "ALL_converters":
            for j in range(len(display.columns)):
                tbl[(i, j)].set_facecolor("#fff3b0")
    png_p = OUT_ROOT / "summary/aao_demographic_summary.png"
    if not args.dry_run:
        fig.savefig(png_p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {png_p}")
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
