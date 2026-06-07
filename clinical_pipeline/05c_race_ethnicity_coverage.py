r"""
05c_race_ethnicity_coverage.py
================================
Race + ethnicity breakdown for the canonical 571-patient clinical cohort
(post-exclusion tabular baseline). One standalone descriptive table —
NOT split by seed/conversion-group (this is a whole-cohort composition table).

Two sources, both for the SAME 571 subjects:
  - RACE      : the `Ethnicity` column of the tabular baseline split CSVs,
                which despite its name is mapped from ADNI `PTRACCAT` (RACE),
                see 01b_build_clinical_csv.py:113.
  - ETHNICITY : ADNI `PTETHCAT` (Hispanic/Latino) joined from PTDEMOG.csv by
                Patient_ID (earliest VISDATE row per subject).

Each subject is classified JOINTLY into ONE mutually-exclusive category so the
breakdown is a true partition that sums to 571: Hispanic/Latino (any race) is
pulled out as its own category, and everyone else is counted by race. The 1
subject with Unknown race was Hispanic, so the partition has no Unknown-race
row. Underrepresented = everyone except non-Hispanic White (= non-White race
OR Hispanic/Latino).

Out: clinical_pipeline/outputs/race/race_ethnicity_571.{csv,png,pdf}

Run:  python clinical_pipeline/05c_race_ethnicity_coverage.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "outputs" / "race"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COHORT_BASE = Path(
    "D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/"
    "tabular/baseline/seed_0")
PTDEMOG = Path("D:/ADNI_BIDS_project/sourcedata/clinical/PTDEMOG.csv")

# Display order for the JOINT race/ethnicity partition (sums to 571).
# "Hispanic or Latino (any race)" is pulled out first; everyone else by race.
PARTITION_ORDER = [
    "White (non-Hispanic)",
    "Hispanic or Latino (any race)",
    "Black or African American",
    "Asian",
    "More than one race",
    "American Indian or Alaskan Native",
    "Native Hawaiian or Other Pacific Islander",
    "Unknown race",
]
HISPANIC = "Hispanic or Latino (any race)"


def load_cohort() -> pd.DataFrame:
    """571 unique subjects from the post-exclusion tabular baseline (seed_0)."""
    parts = [pd.read_csv(COHORT_BASE / f"{sp}.csv", low_memory=False)
             for sp in ("train", "val", "test")]
    df = pd.concat(parts, ignore_index=True).drop_duplicates("Patient_ID")
    df["Patient_ID"] = df["Patient_ID"].astype(str).str.strip()
    return df


def join_ethnicity(df: pd.DataFrame) -> pd.Series:
    """PTETHCAT (Hispanic/Latino) per subject, earliest VISDATE row."""
    demo = pd.read_csv(PTDEMOG, low_memory=False)
    demo["PTID"] = demo["PTID"].astype(str).str.strip()
    demo["_v"] = pd.to_datetime(demo.get("VISDATE"), errors="coerce")
    demo = demo.sort_values("_v").drop_duplicates("PTID", keep="first").set_index("PTID")
    eth = demo.reindex(df["Patient_ID"])["PTETHCAT"].astype(str).str.strip()
    eth.index = df["Patient_ID"].values
    return eth


def main() -> None:
    df = load_cohort()
    n = len(df)
    race = df["Ethnicity"].astype(str).str.strip().values
    eth = join_ethnicity(df).values
    hisp = (eth == "Hispanic or Latino")

    # Joint per-subject classification → ONE mutually-exclusive category.
    # Hispanic/Latino (any race) pulled out first; everyone else by race.
    cat = np.where(hisp, HISPANIC, race)
    cat = np.where(cat == "White", "White (non-Hispanic)", cat)
    cat = np.where(cat == "Unknown", "Unknown race", cat)
    counts = pd.Series(cat).value_counts()

    def pct(k):
        return 100 * k / n

    underrep = n - int(counts.get("White (non-Hispanic)", 0))  # non-White or Hispanic

    # ── tidy CSV ────────────────────────────────────────────────────────────
    rows = []
    for c in PARTITION_ORDER:
        k = int(counts.get(c, 0))
        if c == "Unknown race" and k == 0:
            continue                                   # no Unknown-race subjects remain
        rows.append({"Category": c, "N": k, "Pct_of_571": round(pct(k), 1)})
    rows.append({"Category": "Non-White",
                 "N": underrep, "Pct_of_571": round(pct(underrep), 1)})
    rows.append({"Category": "Total cohort", "N": n, "Pct_of_571": 100.0})
    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_DIR / "race_ethnicity_571.csv", index=False)
    print(f"Saved CSV → {OUT_DIR / 'race_ethnicity_571.csv'}")
    print(out_df.to_string(index=False))
    part_sum = sum(int(counts.get(c, 0)) for c in PARTITION_ORDER)
    print(f"[check] partition Σ = {part_sum}  ({'✓' if part_sum == n else 'MISMATCH'})")

    # ── styled figure (house style: serif, bordered, ruled — matches 01e/04d) ─
    # Build a flat render list: (label, n_str, pct_str, kind)
    #   kind: 'data' | 'subtotal' | 'total'
    render = []
    for c in PARTITION_ORDER:
        k = int(counts.get(c, 0))
        if c == "Unknown race" and k == 0:
            continue
        render.append((c, str(k), f"{pct(k):.1f}%", "data"))
    render.append(("Non-White", str(underrep), f"{pct(underrep):.1f}%", "subtotal"))
    render.append(("Total cohort", str(n), "100%", "total"))

    COL_W = [5.60, 1.05, 1.45]                 # Category | N | % of 571
    LEFT, RIGHT_PAD = 0.30, 0.30
    TITLE_H, HEAD_H, ROW_H = 0.62, 0.36, 0.40
    TOP_PAD, BOT_PAD = 0.12, 0.12
    n_rows = len(render)
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * n_rows + BOT_PAD

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    col_left = [LEFT]
    for w in COL_W:
        col_left.append(col_left[-1] + w)
    RIGHT = col_left[-1]
    col_cx = [(col_left[i] + col_left[i + 1]) / 2 for i in range(len(COL_W))]

    def hline(y, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [y, y], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD
    y_title_top = y
    y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, (y_title_top + y) / 2,
            "Clinical Cohort Ethnicity (n=571)",
            ha="center", va="center", fontsize=11, fontweight="bold", linespacing=1.4)
    hline(y_title_top, lw=1.5)
    hline(y, lw=1.2)

    # header row
    headers = ["Category", "N", "% of 571"]
    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    ax.text(col_left[0] + 0.08, ymid, headers[0], ha="left", va="center",
            fontsize=9.5, fontstyle="italic")
    for j in (1, 2):
        ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for (lab, nstr, pstr, kind) in render:
        if kind in ("subtotal", "total"):
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y
        y -= ROW_H
        ymid = (yr_top + y) / 2
        bold = kind in ("subtotal", "total")
        ax.text(col_left[0] + 0.08, ymid, lab, ha="left", va="center",
                fontsize=9.5, fontweight="bold" if bold else "normal")
        for j, val in zip((1, 2), (nstr, pstr)):
            ax.text(col_cx[j], ymid, val, ha="center", va="center",
                    fontsize=9.5, fontweight="bold" if bold else "normal")
    BOTTOM = y

    # dashed vertical divider after the Category column
    ax.plot([col_left[1], col_left[1]], [BOTTOM, y_data_top],
            color="black", linewidth=0.7, linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5,
                               zorder=5))

    fig.savefig(OUT_DIR / "race_ethnicity_571.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "race_ethnicity_571.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved PNG → {OUT_DIR / 'race_ethnicity_571.png'}")


if __name__ == "__main__":
    main()
