r"""
22_build_conversion_labels.py
=============================
Build CLEAN, strict per-patient AD/MCI conversion labels + continuous
time-to-event from the per-visit longitudinal master, replacing the
methodologically weak `Label_ever_convert` (which counted prevalent
baseline-MCI and transient/reverting non-CN visits as "converters").

Source (one row per visit):
  D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\
  longitudinal\master_clinical_tabular.csv
  cols used: Patient_ID, VISCODE_long, Date, Label_bl_multi, Label_visit_diag

Per patient (visits sorted by Date; baseline = VISCODE_long=='bl', else
earliest visit; last_dx/last_date = final non-null Label_visit_diag):
  bl_dx                  baseline diagnosis (AD here = prevalent)
  last_dx, followup_years (last_date − baseline)/365.25
  first_conversion_MCI   yr to first post-baseline MCI visit  (NaN if none)
  first_conversion_AD    yr to first post-baseline AD  visit  (NaN if none)
  reached_AD_ever        ≥1 post-baseline AD visit EVER (raw/noisy; transient
                         counts) — diagnostic only, NOT a strict label. The
                         strict `ever_conversion_AD` ⊆ this; the difference is
                         exactly the AD→MCI/CN reverters (touched AD but not
                         AD at last visit).
  -- survival endpoints (incident only; prevalent baseline-AD EXCLUDED) --
  ever_conversion_AD      bl∈{CN,MCI} & reaches AD & last_dx==AD   (strict)
  ever_conversion_MCI     bl==CN & reaches MCI/AD & last_dx∈{MCI,AD}(strict)
  ever_conversion_AD_or_MCI  bl∈{CN,MCI} & (ever_conversion_AD|ever_conversion_MCI)
  progressive_MCI         bl==MCI & reaches AD & last_dx==AD
  stable_MCI              bl==MCI & last_dx==MCI & never AD  (MCI at BOTH
                          baseline and last visit; MCI→CN is reversion, NOT
                          stable — see reversion_MCI_to_CN)
  stable_CN               bl==CN & last_dx==CN   (CN at BOTH baseline and
                          last visit — the cognitively-stable controls)
  reversion_MCI_to_CN     bl==MCI & last_dx==CN   (excludable noise)
  -- classification targets (genotype@baseline → has/gets AD; prevalent
     baseline-AD INCLUDED as positives alongside the converters) --
  ad_case                 bl==AD  OR  ever_conversion_AD==1
  ad_or_mci_case          bl==AD  OR  ever_conversion_AD_or_MCI==1

"Confirmed at last visit" removes diagnostic reversion (e.g. CN→MCI→CN,
~40/378 ever_convert=1 rows in seed_0 train ended CN). Output feeds the
Cox survival models (23, incident endpoints) and the strict-label
genetic-baseline classification re-run (20, ad_case targets).

Usage:  conda run -n snp python snp_pipeline/22_build_conversion_labels.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

LONG = Path("D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/"
            "tabular/longitudinal/master_clinical_tabular.csv")
SPLITS = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
              "no_cdr_stratified_ever_convert/tabular/baseline")
OUT = Path("D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels")
DX = {"CN", "MCI", "AD"}


def build() -> pd.DataFrame:
    df = pd.read_csv(LONG, low_memory=False,
                     usecols=["Patient_ID", "VISCODE_long", "Date",
                              "Label_bl_multi", "Label_visit_diag"])
    df["Patient_ID"] = df["Patient_ID"].astype(str)
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df = df.dropna(subset=["Date"]).sort_values(["Patient_ID", "Date"])

    rows = []
    for pid, g in df.groupby("Patient_ID", sort=False):
        g = g.sort_values("Date")
        bl = g[g["VISCODE_long"].astype(str) == "bl"]
        bl_date = (bl["Date"].iloc[0] if len(bl) else g["Date"].iloc[0])
        bl_dx = str((bl["Label_bl_multi"].iloc[0] if len(bl)
                     else g["Label_bl_multi"].iloc[0])).strip().upper()
        dxg = g[g["Label_visit_diag"].astype(str).str.upper().isin(DX)].copy()
        dxg["dx"] = dxg["Label_visit_diag"].astype(str).str.upper()
        if dxg.empty:
            continue
        last_dx = dxg["dx"].iloc[-1]
        last_date = dxg["Date"].iloc[-1]
        post = dxg[dxg["Date"] > bl_date]
        fm = post.loc[post["dx"] == "MCI", "Date"]
        fa = post.loc[post["dx"] == "AD", "Date"]
        first_mci = fm.min() if len(fm) else pd.NaT
        first_ad = fa.min() if len(fa) else pd.NaT
        yr = lambda d: ((d - bl_date).days / 365.25
                        if pd.notna(d) else np.nan)
        reaches_ad = pd.notna(first_ad)
        reaches_mci = pd.notna(first_mci) or reaches_ad
        ev_ad = int(bl_dx in ("CN", "MCI") and reaches_ad and last_dx == "AD")
        ev_mci = int(bl_dx == "CN" and reaches_mci and last_dx in ("MCI", "AD"))
        ev_admci = int(bl_dx in ("CN", "MCI") and (ev_ad or ev_mci))
        prog_mci = int(bl_dx == "MCI" and reaches_ad and last_dx == "AD")
        # stable = MCI at BOTH baseline and last visit, never AD. MCI→CN is a
        # reversion (flagged separately so it can be excluded), NOT stable.
        stab_mci = int(bl_dx == "MCI" and last_dx == "MCI" and not reaches_ad)
        # CN at BOTH baseline and last visit — the cognitively-stable
        # controls (same convention as clinical_pipeline/05).
        stab_cn = int(bl_dx == "CN" and last_dx == "CN")
        rev_mci_cn = int(bl_dx == "MCI" and last_dx == "CN")
        rows.append({
            "Patient_ID": pid, "bl_dx": bl_dx, "last_dx": last_dx,
            "n_visits": int(g["VISCODE_long"].nunique()),
            "followup_years": round(yr(last_date), 4),
            "first_conversion_MCI": (round(yr(first_mci), 4)
                                     if pd.notna(first_mci) else np.nan),
            "first_conversion_AD": (round(yr(first_ad), 4)
                                    if pd.notna(first_ad) else np.nan),
            "reached_AD_ever": int(reaches_ad),
            # incident survival endpoints (prevalent baseline-AD excluded)
            "ever_conversion_AD": ev_ad,
            "ever_conversion_MCI": ev_mci,
            "ever_conversion_AD_or_MCI": ev_admci,
            "progressive_MCI": prog_mci,
            "stable_MCI": stab_mci,
            "stable_CN": stab_cn,
            "reversion_MCI_to_CN": rev_mci_cn,
            # classification targets — prevalent baseline-AD INCLUDED as
            # positives alongside the incident converters (genotype@baseline
            # → has/gets AD; case/control, prevalent = cleanest positive)
            "ad_case": int(bl_dx == "AD" or ev_ad == 1),
            "ad_or_mci_case": int(bl_dx == "AD" or ev_admci == 1),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    lab = build()
    fout = args.out / "conversion_labels.tsv"
    lab.to_csv(fout, sep="\t", index=False)
    print(f"[out] {fout}  ({len(lab)} patients)")

    # ── QC ──────────────────────────────────────────────────────────────────
    nAD = int((lab["bl_dx"] == "AD").sum())
    print("\n[QC] baseline-dx counts (all patients in longitudinal master):")
    print(lab["bl_dx"].value_counts().to_dict())
    print(f"\n[QC] incident SURVIVAL endpoints (prevalent baseline-AD={nAD} "
          f"excluded; at-risk = {len(lab)-nAD}):")
    for c in ("ever_conversion_AD", "ever_conversion_MCI",
              "ever_conversion_AD_or_MCI", "progressive_MCI",
              "stable_MCI", "stable_CN", "reversion_MCI_to_CN"):
        print(f"  {c:26s} positives = {int(lab[c].sum())}")
    print("\n[QC] CLASSIFICATION targets (prevalent baseline-AD INCLUDED "
          "as positives):")
    inc_ad = int(((lab["bl_dx"] != "AD") & (lab["ever_conversion_AD"] == 1)).sum())
    inc_am = int(((lab["bl_dx"] != "AD")
                  & (lab["ever_conversion_AD_or_MCI"] == 1)).sum())
    print(f"  ad_case        positives = {int(lab['ad_case'].sum())}  "
          f"(= {nAD} prevalent + {inc_ad} incident ever_conversion_AD)")
    print(f"  ad_or_mci_case positives = {int(lab['ad_or_mci_case'].sum())}  "
          f"(= {nAD} prevalent + {inc_am} incident ever_conversion_AD_or_MCI)")

    mci = lab[lab["bl_dx"] == "MCI"]
    p, s, r = (int(mci["progressive_MCI"].sum()),
               int(mci["stable_MCI"].sum()),
               int(mci["reversion_MCI_to_CN"].sum()))
    tr = len(mci) - p - s - r
    print(f"\n[QC] baseline-MCI partition: {len(mci)} = progressive {p} + "
          f"stable {s} + reversion→CN {r} + transient/other {tr} "
          f"(transient = touched AD then ended MCI, or no usable follow-up)")
    rev = lab[(lab["reached_AD_ever"] == 1) & (lab["ever_conversion_AD"] == 0)
              & (lab["bl_dx"] != "AD")]
    print(f"\n[QC] reached_AD_ever vs ever_conversion_AD: "
          f"reached_AD_ever={int((lab['reached_AD_ever']==1).sum())}, "
          f"strict ever_conversion_AD={int(lab['ever_conversion_AD'].sum())}; "
          f"delta = {len(rev)} non-prevalent patients touched an AD visit but "
          f"were NOT AD at last visit (AD→MCI/CN reverters — the strict "
          f"label's noise, kept only as the reached_AD_ever diagnostic).")

    # ── spot-checks (reproducible) ──────────────────────────────────────────
    li = lab.set_index("Patient_ID")
    cols = ["bl_dx", "last_dx", "ever_conversion_AD", "ever_conversion_MCI",
            "progressive_MCI", "stable_MCI", "reversion_MCI_to_CN",
            "ad_case", "first_conversion_AD"]
    print("\n[QC] spot-checks:")
    for pid in ("068_S_2184", "129_S_0778"):
        if pid in li.index:
            print(f"  {pid}: " + li.loc[pid, cols].to_dict().__repr__())
    adp = li[li["bl_dx"] == "AD"]
    if len(adp):
        pid = adp.index[0]
        print(f"  prevalent-AD {pid}: " + li.loc[pid, cols].to_dict().__repr__()
              + "  (expect ad_case=1, ever_conversion_AD=0)")

    # Per-seed: how many split patients have a label (membership unchanged)
    for sd in (0, 1, 2):
        try:
            ids = set()
            for sp in ("train", "val", "test"):
                ids |= set(pd.read_csv(SPLITS / f"seed_{sd}" / f"{sp}.csv",
                                       usecols=["Patient_ID"])["Patient_ID"]
                           .astype(str))
        except FileNotFoundError:
            continue
        cov = lab["Patient_ID"].isin(ids)
        print(f"[QC] seed_{sd}: {cov.sum()}/{len(ids)} split patients "
              f"have a conversion label "
              f"({len(ids)-cov.sum()} missing from longitudinal master)")


if __name__ == "__main__":
    main()
