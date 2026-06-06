r"""
05b_demographic_coverage_reconciled.py
=======================================
Reconciled demographic coverage table — same source data and house style
as `05_demographic_coverage.py`, but with the user's reconciled label set
that splits the old `ever_conversion_AD` (mis-labelled "progressive MCI"
in the legacy table, which actually conflated MCI→AD with CN→AD) and
unifies all reversions into ONE row (no double counting).

Reconciled groups (in display order):

  Full at-risk cohort           reference — n=616
  AD_bl                         baseline-AD (n=36)
  pMCI                          baseline-MCI → AD sustained (n=121)
  CN_to_AD                      baseline-CN → AD sustained (n=25)
  ever_conversion_AD            pMCI ∪ CN_to_AD (n=146; composite)
  ad_cases                      AD_bl ∪ ever_conversion_AD (n=182)
  CN_to_MCI                     baseline-CN → MCI ever (= ever_conversion_MCI;
                                n=60; INCLUDES the 25 who went CN→MCI→AD)
  ever_conversion_MCI_AD        AD_bl ∪ ever_conversion_AD ∪ CN_to_MCI
                                (= ad_or_mci_case, n=217)
  sMCI                          stable_MCI, MCI→MCI (n=203)
  sCN                           stable_CN, CN→CN (n=152)
  Reversions                    ALL reversions in one category, no double
                                counting (= Excluded flag from
                                conversion_labels.tsv; n=44)

Note on overlap: the table is **not a strict partition**; AD_bl, pMCI,
CN_to_AD, ever_conversion_AD, ad_cases overlap with each other by design
(they're composites of each other). The "partition rows" (= rows whose
counts sum to the full 616-subject cohort) are:
  AD_bl + pMCI + CN_to_AD + pCN_to_MCI_only + sMCI + sCN + Reversions
= 36 + 121 + 25 + 35 + 203 + 152 + 44 = 616  ✓
where pCN_to_MCI_only = CN_to_MCI − CN_to_AD = 60 − 25 = 35.
(The 25 CN→MCI→AD overlap shows up once in CN_to_AD and counts in the
overlapping descriptive groups CN_to_MCI / ever_conversion_MCI_AD.)

Outputs to Transformers_XAI\clinical_pipeline\outputs\demographic_coverage\:
  demographic_coverage_table_reconciled.{png,pdf,csv}

Usage:  conda run -n snp python clinical_pipeline/05b_demographic_coverage_reconciled.py
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

# Repo root for bidsification.exclusions import
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT))
from bidsification.exclusions import is_excluded_subject  # noqa: E402

SPLIT_PRE  = Path("D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified/"
                  "tabular/baseline")
SPLIT_POST = Path("D:/ADNI_BIDS_project/derivatives/clinical/"
                  "no_cdr_stratified_post_exclusion/tabular/baseline")
LABELS = Path("D:/ADNI_SNP_Omni2.5M_20140220/conversion_labels/"
              "conversion_labels.tsv")
OUT = Path(__file__).parent / "outputs" / "demographic_coverage"
SEEDS = [0, 1, 2]
SPLITS = ["train", "val", "test"]

# (group_key, pretty_label).  group_key is either:
#   - column name in conversion_labels.tsv (boolean / 0-1 column)
#   - "full_cohort" special-cased
GROUPS = [
    ("full_cohort",            "Full at-risk cohort"),
    ("AD_bl",                  "AD_bl  (baseline AD)"),
    ("pMCI",                   "pMCI  (MCI → AD sustained)"),
    ("CN_to_AD",               "CN_to_AD  (CN → AD sustained)"),
    ("ever_conversion_AD",     "ever_conversion_AD  (pMCI ∪ CN_to_AD)"),
    ("ad_case",                "ad_cases  (AD_bl ∪ ever_conversion_AD)"),
    ("CN_to_MCI",              "CN_to_MCI  (= ever_conversion_MCI)"),
    ("ad_or_mci_case",         "ever_conversion_MCI_AD  (CN→MCI/AD ∪ MCI→AD ∪ AD_bl)"),
    ("stable_MCI",             "sMCI  (stable MCI, MCI→MCI)"),
    ("stable_CN",              "sCN  (stable CN, CN→CN)"),
    ("Excluded",               "Reversions  (all reversions, no double counting)"),
]

# Cleaned copy: non-overlapping partition (drops the composite overlap rows
# ever_conversion_AD / ad_cases / ever_conversion_MCI_AD). pCN_to_MCI_only =
# baseline-CN -> MCI but NOT -> AD (= CN_to_MCI - CN_to_AD), so the CN->AD cases
# are not double-counted with the pCN->AD row. Sums to 616.
CLEANED_GROUPS = [
    ("full_cohort",      "Full Cohort"),
    ("AD_bl",            "AD at baseline visit"),
    ("pMCI",             "progressive MCI (pMCI)"),
    ("stable_MCI",       "stable MCI (sMCI)"),
    ("CN_to_AD",         "progressive CN (pCN, to AD)"),
    ("pCN_to_MCI_only",  "progressive CN (pCN, to MCI)"),
    ("stable_CN",        "stable CN (sCN)"),
    ("Excluded",         "Reversions (excluded)"),
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
                "e2c": "", "e4c": ""}
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
        "e2c": f"{100*(a2 >= 1).mean():.0f}%" if a2.notna().any() else "",
        "e4c": f"{100*(a4 >= 1).mean():.0f}%" if a4.notna().any() else "",
    }


def _cell(s: dict, cleaned: bool = False) -> str:
    if s["N"] == 0:
        return "—"
    if cleaned:
        return (f"N={s['N']}  age {s['age']}  F {s['F']}\n"
                f"ε2+ {s['e2c']}  ε4+ {s['e4c']}")
    return (f"N={s['N']}  age {s['age']}  F {s['F']}\n"
            f"ε2/3/4 {s['e2']}/{s['e3']}/{s['e4']}  ε4+ {s['e4c']}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--post-exclusion", action="store_true",
                    help="Use no_cdr_stratified_post_exclusion splits, drop "
                         "diagnostic-reversion + corrupted-MRI subjects, and "
                         "suffix output filenames with _post_exclusion.")
    ap.add_argument("--cleaned", action="store_true",
                    help="Render a cleaned copy: non-overlapping partition rows, "
                         "ε2+/ε4+ carrier %% (no dosage), trimmed title/subtitle, "
                         "no footnote -> demographic_cleaned.{png,pdf,csv}.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    cleaned = bool(args.cleaned)
    groups = CLEANED_GROUPS if cleaned else GROUPS
    post = bool(args.post_exclusion)
    suffix = "_post_exclusion" if post else ""
    SPLIT = SPLIT_POST if post else SPLIT_PRE
    print(f"[mode] post_exclusion={post}  SPLIT_DIR={SPLIT}  SUFFIX='{suffix}'")

    lab = pd.read_csv(LABELS, sep="\t")
    lab["Patient_ID"] = lab["Patient_ID"].astype(str)
    if post:
        keep_mask = ~lab["Patient_ID"].apply(
            lambda p: is_excluded_subject(p, include_diagnostic_reversion=True))
        n_drop = int((~keep_mask).sum())
        lab = lab[keep_mask].reset_index(drop=True)
        print(f"  POST_EXCLUSION: dropped {n_drop} subjects from labels → {len(lab)} remain")
    # Sanity: confirm the reconciled columns are present (script 07 must have
    # been run; the user updated conversion_labels.tsv earlier).
    required = ["AD_bl", "pMCI", "CN_to_AD", "ever_conversion_AD",
                "ad_case", "CN_to_MCI", "ad_or_mci_case",
                "stable_MCI", "stable_CN", "Excluded"]
    miss = [c for c in required if c not in lab.columns]
    if miss:
        raise SystemExit(
            f"[FATAL] conversion_labels.tsv missing columns: {miss}.\n"
            "Run clinical_pipeline/07_conversion_group_extended.py first.")
    lab = lab.set_index("Patient_ID")

    # Per-cohort sanity print
    n_total = len(lab)
    print(f"\n[cohort] N total in conversion_labels.tsv = {n_total}")
    for col in ("AD_bl", "pMCI", "CN_to_AD", "ever_conversion_AD",
                "ad_case", "CN_to_MCI", "ad_or_mci_case",
                "stable_MCI", "stable_CN", "Excluded"):
        print(f"  {col:<22s} = {int(lab[col].sum()):4d}")
    # Partition arithmetic check (rows that sum to N total)
    pCN_to_MCI_only = int(
        (lab["CN_to_MCI"].astype(int) & ~lab["CN_to_AD"].astype(bool)).sum())
    print(f"  pCN_to_MCI_only (= CN_to_MCI − CN_to_AD) = {pCN_to_MCI_only}")
    partition_sum = (int(lab["AD_bl"].sum())
                      + int(lab["pMCI"].sum())
                      + int(lab["CN_to_AD"].sum())
                      + pCN_to_MCI_only
                      + int(lab["stable_MCI"].sum())
                      + int(lab["stable_CN"].sum())
                      + int(lab["Excluded"].sum()))
    print(f"  partition Σ = AD_bl + pMCI + CN_to_AD + pCN_to_MCI_only + sMCI "
          f"+ sCN + Reversions = {partition_sum}  "
          f"({'✓' if partition_sum == n_total else 'MISMATCH'})")

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
        if key == "pCN_to_MCI_only":
            cn_mci = pd.to_numeric(dm.get("CN_to_MCI"), errors="coerce") == 1
            cn_ad  = pd.to_numeric(dm.get("CN_to_AD"), errors="coerce") == 1
            return dm[cn_mci & ~cn_ad]
        return dm[pd.to_numeric(dm.get(key), errors="coerce") == 1]

    # ── tidy long CSV ──────────────────────────────────────────────────────
    rows = []
    for (gk, gp) in groups:
        for sd in SEEDS:
            for sp in SPLITS:
                dm = demog.get(sd, {}).get(sp)
                if dm is None:
                    continue
                rows.append({"group": gp, "seed": sd, "split": sp,
                             **_summ(_grp(dm, gk))})
    long_df = pd.DataFrame(rows)
    _csv_stem = "demographic_cleaned" if cleaned else "demographic_coverage_table_reconciled"
    long_df.to_csv(args.out / f"{_csv_stem}{suffix}.csv", index=False)
    print(f"[out] {args.out/(_csv_stem+suffix+'.csv')} ({len(long_df)} rows)")

    # ── figure (diagnostic_coverage house style) ───────────────────────────
    COL_W = [3.0, 3.0, 3.0, 3.0]
    LEFT, RIGHT_PAD = 0.25, 0.25
    GROUP_COL_W = 3.4   # wider: reconciled labels are longer
    widths = [GROUP_COL_W] + COL_W[:3]
    fig_w = LEFT + sum(widths) + RIGHT_PAD
    TITLE_H, SUB_H, H1, H2 = 0.66, 0.46, 0.30, 0.28
    R_GRP, R_SEED = 0.30, 0.52
    TOP_PAD, BOT_PAD = 0.12, 0.12
    FOOT_H = 0.10 if cleaned else 1.20
    n_grp, n_seed = len(groups), len(groups) * len(SEEDS)
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
    _title = ("ADNI Demographic Coverage by Conversion Group, Seed and Split"
              if cleaned else
              "ADNI Demographic Coverage by Conversion Group, Seed and Split  "
              "— Reconciled")
    ax.text((LEFT + RIGHT) / 2, (y_top + y) / 2, _title,
            ha="center", va="center", fontsize=12, fontweight="bold")
    yt = y
    y -= SUB_H
    _subtitle = (
        "N subjects · age (mean±sd) · %Female · ε2+ (% ε2-carrier) · "
        "ε4+ (% ε4-carrier, ε4 dosage ≥ 1).   "
        "Cohort stratified by baseline diagnosis CN/MCI/AD and "
        "80/10/10 splits (seeds 0/1/2)."
        if cleaned else
        "Cell = N subjects · age (mean±sd) · %Female · "
        "ε2/ε3/ε4 (mean APOE ε-allele dosage 0–2) · "
        "ε4+ (% ε4-carrier, ε4 dosage ≥ 1).   "
        "no_cdr_stratified 80/10/10 splits (seeds 0/1/2).  "
        "Reversions row is ALL reversions in ONE category — no double counting.")
    ax.text((LEFT + RIGHT) / 2, (yt + y) / 2, _subtitle,
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

    for (gk, gp) in groups:
        yt = y
        y -= R_GRP
        ax.text(cl[0] + 0.06, (yt + y) / 2, gp, ha="left", va="center",
                fontsize=9.5, fontweight="bold")
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
                txt = _cell(_summ(_grp(dm, gk)), cleaned) if dm is not None else "—"
                ax.text(cx[1 + j], ym, txt, ha="center", va="center",
                        fontsize=7.0, color="#222222", linespacing=1.35)
    BOTTOM = y
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
        "train+val+test N. "
        "GROUPS — AD_bl = baseline AD; pMCI = baseline-MCI → AD sustained; "
        "CN_to_AD = baseline-CN → AD sustained (= sustained CN→MCI→AD); "
        "ever_conversion_AD = pMCI ∪ CN_to_AD (CN/MCI→AD sustained); "
        "ad_cases = AD_bl ∪ ever_conversion_AD; "
        "CN_to_MCI = baseline-CN that ever reached MCI sustained "
        "(= ever_conversion_MCI; INCLUDES the CN_to_AD subset that went "
        "via MCI — overlaps ever_conversion_AD); "
        "ever_conversion_MCI_AD = AD_bl ∪ ever_conversion_AD ∪ CN_to_MCI "
        "(= ad_or_mci_case; baseline-AD or any sustained progression); "
        "sMCI = stable MCI (MCI baseline AND last visit); "
        "sCN = stable CN (CN baseline AND last visit); "
        "Reversions = all reversions in ONE category, no double counting "
        "(MCI→CN reversion ∪ reached-AD-reverted ∪ other non-sustained "
        "transitions; = `Excluded` flag in conversion_labels.tsv). "
        "Note the table is NOT a strict partition: AD_bl, pMCI, CN_to_AD, "
        "ever_conversion_AD, ad_cases, CN_to_MCI, ever_conversion_MCI_AD "
        "are composites that overlap each other. "
        "Partition rows that sum to N total: "
        "AD_bl + pMCI + CN_to_AD + (CN_to_MCI − CN_to_AD) + sMCI + sCN + "
        "Reversions  =  36 + 121 + 25 + 35 + 203 + 152 + 44  =  616. "
        "Splits 80/10/10 by Patient_ID, seeds 0/1/2."
    )
    if not cleaned:
        ax.text(LEFT, BOTTOM - 0.16,
                "\n".join(textwrap.wrap(foot, width=int((RIGHT - LEFT) * 16))),
                ha="left", va="top", fontsize=7.5)

    plt.tight_layout(pad=0.1)
    stem = "demographic_cleaned" if cleaned else "demographic_coverage_table_reconciled"
    fig.savefig(args.out / f"{stem}{suffix}.png", bbox_inches="tight", dpi=300)
    fig.savefig(args.out / f"{stem}{suffix}.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"[out] {args.out/(stem+suffix+'.png')} + .pdf")


if __name__ == "__main__":
    main()
