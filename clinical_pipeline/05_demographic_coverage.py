r"""
05_demographic_coverage.py
==========================
Per stratified-baseline seed (0/1/2) AND split (train/val/test) demographic
comparison across the clinical conversion groups, rendered as a figure
table in the diagnostic_coverage house style (PNG + PDF + CSV).

Cohort / fold membership + demographics from the no_cdr_stratified baseline
splits:
  D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\baseline
    seed_{0,1,2}\{train,val,test}.csv
  (Patient_ID, Date, YOB, Sex, APOE_Alleles, APOE4_Dosage)

Group labels joined from script-22 conversion_labels.tsv:
  ever_conversion_AD        CN/MCI → AD, sustained (AD at last visit)
  ad_case                   prevalent baseline-AD OR ever_conversion_AD
  ever_conversion_AD_or_MCI strict →AD OR sustained CN→MCI
  ad_or_mci_case            prevalent baseline-AD OR ever_conversion_AD_or_MCI
  stable_MCI                MCI → MCI (bl & last both MCI, never AD)
  stable_CN                 CN → CN  (bl & last both CN)
  + full at-risk cohort (reference)

Per (group × seed × split): N · age(mean±sd) · %Female · mean ε2/ε3/ε4
dosage · %ε4-carrier.

Outputs to Transformers_XAI\clinical_pipeline\outputs\demographic_coverage\:
  demographic_coverage_table.{png,pdf,csv}   (figure + tidy long CSV)

Usage:  conda run -n snp python clinical_pipeline/05_demographic_coverage.py
"""
from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

SPLIT = Path("D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/"
             "tabular/baseline")
LABELS = Path("D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/"
              "conversion_labels.tsv")
OUT = Path(__file__).parent / "outputs" / "demographic_coverage"
SEEDS = [0, 1, 2]
SPLITS = ["train", "val", "test"]
GROUPS = [
    ("full_cohort", "Full at-risk cohort"),
    ("ever_conversion_AD", "ever_conversion_AD (progressive MCI)"),
    ("ad_case", "ad_case = ever_conversion_AD + AD_bl"),
    ("ever_conversion_MCI", "ever_conversion_MCI (CN→MCI)"),
    ("ever_conversion_AD_or_MCI", "ever_conversion_AD_or_MCI"),
    ("ad_or_mci_case", "ad_or_mci_case = ever_conversion_AD_or_MCI + AD_bl"),
    ("stable_MCI", "stable MCI (MCI→MCI)"),
    ("stable_CN", "stable CN (CN→CN)"),
    # ── "Other" = cohort − (ad_case ∪ stable_MCI ∪ stable_CN ∪
    #    sustained-CN→MCI-no-AD), the genuine residual reverters ─────────
    ("oth_rev_MCI_CN", "Other: MCI→CN reversion"),
    ("oth_AD_reverted", "Other: reached AD but reverted"),
]


def _demog(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["Patient_ID"] = d["Patient_ID"].astype(str)
    yr = pd.to_datetime(d["Date"], errors="coerce").dt.year
    al = d.get("APOE_Alleles", pd.Series([""] * len(d))).astype(str)
    has = al.str.contains(r"\d", na=False)
    return pd.DataFrame({
        "Patient_ID": d["Patient_ID"].values,
        "age": (yr - pd.to_numeric(d["YOB"], errors="coerce")).values,
        "female": (d["Sex"].astype(str).str.strip().str.lower()
                   .map({"female": 1, "f": 1, "2": 1,
                         "male": 0, "m": 0, "1": 0})).values,
        "apoe2": al.str.count("2").where(has).values,
        "apoe3": al.str.count("3").where(has).values,
        "apoe4": al.str.count("4").where(has).values,
    }).drop_duplicates("Patient_ID")


def _summ(g: pd.DataFrame) -> dict:
    n = len(g)
    if n == 0:
        return {"N": 0, "age": "", "F": "", "e2": "", "e3": "", "e4": "",
                "e4c": ""}
    age = pd.to_numeric(g["age"], errors="coerce").dropna()
    fem = pd.to_numeric(g["female"], errors="coerce")
    a2 = pd.to_numeric(g["apoe2"], errors="coerce")
    a3 = pd.to_numeric(g["apoe3"], errors="coerce")
    a4 = pd.to_numeric(g["apoe4"], errors="coerce")
    return {
        "N": n,
        "age": f"{age.mean():.1f}±{age.std(ddof=1):.1f}" if len(age) else "",
        "F": f"{100*fem.mean():.0f}%" if fem.notna().any() else "",
        "e2": f"{a2.mean():.2f}" if a2.notna().any() else "",
        "e3": f"{a3.mean():.2f}" if a3.notna().any() else "",
        "e4": f"{a4.mean():.2f}" if a4.notna().any() else "",
        "e4c": f"{100*(a4 >= 1).mean():.0f}%" if a4.notna().any() else "",
    }


def _cell(s: dict) -> str:
    if s["N"] == 0:
        return "—"
    return (f"N={s['N']}  age {s['age']}  F {s['F']}\n"
            f"ε2/3/4 {s['e2']}/{s['e3']}/{s['e4']}  ε4+ {s['e4c']}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    lab = pd.read_csv(LABELS, sep="\t")
    lab["Patient_ID"] = lab["Patient_ID"].astype(str)
    bl = lab["bl_dx"].astype(str).str.upper()
    ls = lab["last_dx"].astype(str).str.upper()
    lab["stable_CN"] = ((bl == "CN") & (ls == "CN")).astype(int)

    def _num(c):
        return pd.to_numeric(lab.get(c), errors="coerce").fillna(0)
    in_main = ((_num("ad_case") == 1) | (_num("stable_MCI") == 1)
               | (lab["stable_CN"] == 1))
    rev = (_num("reversion_MCI_to_CN") == 1)
    adrev = (_num("reached_AD_ever") == 1)
    # Sustained CN→MCI that never reached AD is now its own descriptive
    # group (ever_conversion_MCI) — exclude it from "Other" so the residual
    # is the genuine reverters / reached-AD-reverted only.
    cnmci_noad = ((bl == "CN") & (_num("ever_conversion_MCI") == 1)
                  & ~adrev & ~rev)
    covered = in_main | cnmci_noad
    # "Other" reason (priority order), only for rows not otherwise covered:
    reason = pd.Series("", index=lab.index)
    oth = ~covered
    reason[oth & rev] = "oth_rev_MCI_CN"
    reason[oth & (reason == "") & adrev] = "oth_AD_reverted"
    reason[oth & (reason == "")] = "oth_unclassified"
    lab["other_reason"] = reason
    nu = int((reason == "oth_unclassified").sum())
    print(f"[reconcile] cohort = ad_case ∪ stable_MCI ∪ stable_CN ∪ "
          f"sustained-CN→MCI-no-AD ({int(cnmci_noad.sum())}) ∪ "
          f"Other[MCI→CN rev {int((reason=='oth_rev_MCI_CN').sum())}, "
          f"reached-AD-reverted {int((reason=='oth_AD_reverted').sum())}]"
          + (f"  · UNCLASSIFIED={nu} (not rendered)" if nu else ""))
    lab = lab.set_index("Patient_ID")

    # demog[seed][split] -> per-PTID demographics joined to labels
    demog = {}
    for sd in SEEDS:
        demog[sd] = {}
        for sp in SPLITS:
            f = SPLIT / f"seed_{sd}" / f"{sp}.csv"
            if not f.exists():
                continue
            dm = _demog(pd.read_csv(f, low_memory=False)).join(
                lab, on="Patient_ID")
            demog[sd][sp] = dm

    def _grp(dm, key):
        if dm is None or len(dm) == 0:
            return dm.iloc[0:0] if dm is not None else None
        if key == "full_cohort":
            return dm
        if key.startswith("oth_"):
            return dm[dm.get("other_reason") == key]
        return dm[pd.to_numeric(dm.get(key), errors="coerce") == 1]

    # ── tidy long CSV ──────────────────────────────────────────────────────
    rows = []
    for (gk, gp) in GROUPS:
        for sd in SEEDS:
            for sp in SPLITS:
                dm = demog.get(sd, {}).get(sp)
                if dm is None:
                    continue
                rows.append({"group": gp, "seed": sd, "split": sp,
                             **_summ(_grp(dm, gk))})
    long_df = pd.DataFrame(rows)
    long_df.to_csv(args.out / "demographic_coverage_table.csv", index=False)
    print(f"[out] {args.out/'demographic_coverage_table.csv'} "
          f"({len(long_df)} rows)")

    # ── figure (diagnostic_coverage house style) ───────────────────────────
    COL_W = [3.0, 3.0, 3.0, 3.0]            # train | val | test | (group N)
    LEFT, RIGHT_PAD = 0.25, 0.25
    GROUP_COL_W = 2.6
    widths = [GROUP_COL_W] + COL_W[:3]
    fig_w = LEFT + sum(widths) + RIGHT_PAD
    TITLE_H, SUB_H, H1, H2 = 0.66, 0.46, 0.30, 0.28
    R_GRP, R_SEED = 0.30, 0.52
    TOP_PAD, BOT_PAD, FOOT_H = 0.12, 0.12, 0.95
    n_grp, n_seed = len(GROUPS), len(GROUPS) * len(SEEDS)
    fig_h = (TOP_PAD + TITLE_H + SUB_H + H1 + H2 + R_GRP * n_grp
             + R_SEED * n_seed + FOOT_H + BOT_PAD)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    cl = [LEFT]
    for w in widths:
        cl.append(cl[-1] + w)
    RIGHT = cl[-1]
    cx = [(cl[i] + cl[i + 1]) / 2 for i in range(len(widths))]

    def hline(y, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw,
                linestyle=ls, solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD
    y_top = y
    y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, (y_top + y) / 2,
            "ADNI Demographic Coverage by Conversion Group, Seed and Split",
            ha="center", va="center", fontsize=12, fontweight="bold")
    yt = y
    y -= SUB_H
    ax.text((LEFT + RIGHT) / 2, (yt + y) / 2,
            "Cell = N subjects · age (mean±sd) · %Female · "
            "ε2/ε3/ε4 (mean APOE ε-allele dosage 0–2) · "
            "ε4+ (% ε4-carrier, ε4 dosage ≥ 1).   "
            "no_cdr_stratified 80/10/10 splits (seeds 0/1/2); groups from "
            "conversion_labels.tsv.   seed label shows Σ = train+val+test N.",
            ha="center", va="center", fontsize=8.0, fontstyle="italic")
    hline(y_top, lw=1.5)
    hline(y, lw=0.8)

    # header
    yt = y
    y -= H1
    ax.text(cx[0], (yt + y) / 2, "Group / Seed", ha="center", va="center",
            fontsize=10, fontweight="bold")
    for j, sp in enumerate(SPLITS):
        ax.text(cx[1 + j], (yt + y) / 2, sp, ha="center", va="center",
                fontsize=10, fontweight="bold", fontstyle="italic")
    hline(y, lw=1.2)
    y_data_top = y

    for (gk, gp) in GROUPS:
        yt = y
        y -= R_GRP
        ax.text(cl[0] + 0.06, (yt + y) / 2, gp, ha="left", va="center",
                fontsize=9.5, fontweight="bold")
        # full-cohort N for this group (seed 0 union)
        if gk != "full_cohort":
            hline(yt, lw=0.5, ls=(0, (3, 3)))
        for sd in SEEDS:
            yt = y
            y -= R_SEED
            ym = (yt + y) / 2
            tot = sum(len(_grp(demog.get(sd, {}).get(sp), gk))
                      for sp in SPLITS
                      if demog.get(sd, {}).get(sp) is not None)
            ax.text(cl[0] + 0.15, ym, f"seed {sd}  (Σ={tot})",
                    ha="left", va="center", fontsize=7.5,
                    fontstyle="italic", color="#444444")
            for j, sp in enumerate(SPLITS):
                dm = demog.get(sd, {}).get(sp)
                txt = _cell(_summ(_grp(dm, gk))) if dm is not None else "—"
                ax.text(cx[1 + j], ym, txt, ha="center", va="center",
                        fontsize=7.0, color="#222222", linespacing=1.35)
    BOTTOM = y
    # verticals
    ax.plot([cl[1], cl[1]], [BOTTOM, y_data_top], color="black",
            linewidth=0.7, linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT,
                               y_top - BOTTOM, facecolor="none",
                               edgecolor="black", linewidth=1.5, zorder=5))
    foot = (
        "Each cell: N subjects · mean±sd age · %Female · mean APOE ε2/ε3/ε4 "
        "allele dosage · %ε4-carrier (ε4 dosage ≥ 1), for that group within "
        "that seed's train/val/test partition; seed-label Σ = combined "
        "train+val+test N. Groups: ever_conversion_AD = CN/MCI→AD sustained "
        "(AD at last visit; here MCI→AD progressive MCI dominates); "
        "ad_case = ever_conversion_AD + prevalent baseline-AD (AD_bl); "
        "ever_conversion_MCI = baseline-CN that sustainably progress to MCI "
        "(CN→MCI; overlaps ad_case for those who go on to AD); "
        "ad_or_mci_case = ever_conversion_AD_or_MCI + prevalent baseline-AD "
        "(AD_bl) — note stable baseline-MCI are NOT added; "
        "stable MCI = MCI at baseline & last visit; "
        "stable CN = CN at baseline & last visit. The two \"Other\" rows "
        "are the genuine residual = full cohort MINUS (ad_case ∪ stable "
        "MCI ∪ stable CN ∪ sustained-CN→MCI-without-AD): MCI→CN reversion "
        "(reversion_MCI_to_CN); reached AD but reverted (touched AD, last "
        "visit ≠ AD ⇒ not ad_case, not stable). full_cohort = all at-risk; "
        "the descriptive groups overlap (ad_case ⊂ ad_or_mci_case; "
        "ever_conversion_MCI overlaps ad_case), but ad_case + stable MCI + "
        "stable CN + sustained-CN→MCI-without-AD + the two Other rows "
        "partition the cohort. Splits 80/10/10 by Patient_ID, seeds 0/1/2."
    )
    ax.text(LEFT, BOTTOM - 0.16,
            "\n".join(textwrap.wrap(foot, width=int((RIGHT - LEFT) * 16))),
            ha="left", va="top", fontsize=7.5)

    plt.tight_layout(pad=0.1)
    fig.savefig(args.out / "demographic_coverage_table.png",
                bbox_inches="tight", dpi=300)
    fig.savefig(args.out / "demographic_coverage_table.pdf",
                bbox_inches="tight", dpi=300)
    plt.close()
    print(f"[out] {args.out/'demographic_coverage_table.png'} + .pdf")


if __name__ == "__main__":
    main()
