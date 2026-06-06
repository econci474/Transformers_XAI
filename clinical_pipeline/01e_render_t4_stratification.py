r"""
01e_render_t4_stratification.py
================================
Render the T4 converter-cohort stratification table (CSV + PNG/PDF), reusing the
label logic in 01e_build_T4_labels_and_splits.py so the counts always reflect the
real master (not a hardcoded docstring).

One row per T4 class (Label_T4 = 3-class time-to-AD horizon), columns:
  Label_T4 (+ horizon) | n (% of 146) | pMCI | pCN_to_AD     [+ Total row]

Stratification note: the 80/10/10 splits are stratified on Label_T4 ALONE (the 3
horizon classes), NOT on conversion_group — pCN_to_AD is too small (1 subject in
the <3y bucket) to stratify jointly, so its pMCI/pCN_to_AD mix per split is
whatever the horizon-stratification yields.

Out: clinical_pipeline/outputs/encoder_only_post_exclusions/t4_stratified.{csv,png,pdf}

Run:  python clinical_pipeline/01e_render_t4_stratification.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "outputs" / "encoder_only_post_exclusions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Bins match label_T4 exactly: 0 = years<3, 1 = 3<=years<7, 2 = years>=7.
WINDOW = {0: "< 3 yrs", 1: "≥ 3, < 7 yrs", 2: "≥ 7 yrs"}


def _load_01e():
    spec = importlib.util.spec_from_file_location(
        "build_t4", HERE / "01e_build_T4_labels_and_splits.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)            # __main__ guard prevents main() from running
    return mod


def main() -> None:
    b = _load_01e()
    master = pd.read_csv(b.MRI_MASTER, low_memory=False)
    subj = b.build_subject_table(b.compute_label_T4(master))
    n_total = len(subj)                                       # 146

    ct = pd.crosstab(subj["Label_T4"], subj["conversion_group"])
    for g in b.T4_COHORT_GROUPS:                              # ensure both columns exist
        if g not in ct.columns:
            ct[g] = 0

    # ── tidy CSV ────────────────────────────────────────────────────────────
    rows = []
    for cls in sorted(WINDOW):
        n = int((subj["Label_T4"] == cls).sum())
        rows.append({"Label_T4": cls, "Window": WINDOW[cls], "N": n,
                     "Pct_of_146": round(100 * n / n_total, 1),
                     "pMCI": int(ct.loc[cls, "pMCI"]),
                     "pCN_to_AD": int(ct.loc[cls, "pCN_to_AD"])})
    rows.append({"Label_T4": "Total", "Window": "", "N": n_total, "Pct_of_146": 100.0,
                 "pMCI": int(ct["pMCI"].sum()), "pCN_to_AD": int(ct["pCN_to_AD"].sum())})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "t4_stratified.csv", index=False)
    print(f"Saved CSV → {OUT_DIR / 't4_stratified.csv'}")
    print(df.to_string(index=False))

    # ── styled figure (house style: serif, bordered, ruled — matches 04d) ─────
    headers = ["Label_T4", "n (% of 146)", "pMCI", "pCN_to_AD"]
    data = []
    for r in rows:
        is_total = r["Label_T4"] == "Total"
        lab = f"{r['Label_T4']}  ({r['Window']})" if r["Window"] else f"{r['Label_T4']}"
        npct = f"{r['N']} (100%)" if is_total else f"{r['N']} ({r['Pct_of_146']:.1f}%)"
        data.append((lab, npct, str(r["pMCI"]), str(r["pCN_to_AD"]), is_total))

    COL_W = [2.30, 2.10, 1.50, 1.70]
    LEFT, RIGHT_PAD = 0.30, 0.30
    TITLE_H, HEAD_H, ROW_H = 0.96, 0.36, 0.42
    TOP_PAD, BOT_PAD = 0.12, 0.12
    n_rows = len(data)
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
            "Clinical Models: T4 Converter Cohort (n=146, pMCI + pCN_to_AD)\n"
            "splits stratified on the conversion window (Label_T4)\n"
            "80/10/10 splits, seeds 0/1/2",
            ha="center", va="center", fontsize=11, fontweight="bold", linespacing=1.4)
    hline(y_title_top, lw=1.5)
    hline(y, lw=1.2)

    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    ax.text(col_left[0] + 0.08, ymid, headers[0], ha="left", va="center",
            fontsize=9.5, fontstyle="italic")
    for j in (1, 2, 3):
        ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for (lab, npct, pmci, pcn, is_total) in data:
        if is_total:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y
        y -= ROW_H
        ymid = (yr_top + y) / 2
        ax.text(col_left[0] + 0.08, ymid, lab, ha="left", va="center",
                fontsize=9.5, fontweight="bold")
        for j, val in zip((1, 2, 3), (npct, pmci, pcn)):
            ax.text(col_cx[j], ymid, val, ha="center", va="center",
                    fontsize=9.5, fontweight="bold" if is_total else "normal")
    BOTTOM = y

    # dashed vertical divider after the Label_T4 column (mirrors the Strategy divider)
    ax.plot([col_left[1], col_left[1]], [BOTTOM, y_data_top],
            color="black", linewidth=0.7, linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))

    fig.savefig(OUT_DIR / "t4_stratified.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / "t4_stratified.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved PNG → {OUT_DIR / 't4_stratified.png'}")


if __name__ == "__main__":
    main()
