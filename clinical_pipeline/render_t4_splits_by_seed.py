r"""
render_t4_splits_by_seed.py
===========================
Small house-style figure of the T4 converter-cohort stratified splits, broken down BY SEED and FOLD
(train/val/test), per conversion window (Label_T4: <3y / 3-7y / >=7y). Reads the canonical split CSVs
(no recompute).

Out: clinical_pipeline/outputs/tasks/t4_stratification_splits.{csv,png,pdf}
Run: python clinical_pipeline/render_t4_splits_by_seed.py
"""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"

SPLIT_BASE = Path("D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/"
                  "tabular/baseline_T4")
OUT_DIR = Path(__file__).resolve().parent / "outputs" / "tasks"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SEEDS = [0, 1, 2]
FOLDS = ["train", "val", "test"]
WINDOW = {0: "< 3 yrs", 1: "≥ 3, < 7 yrs", 2: "≥ 7 yrs"}


def collect():
    rows = []
    for s in SEEDS:
        for fold in FOLDS:
            d = pd.read_csv(SPLIT_BASE / f"seed_{s}" / f"{fold}.csv")
            vc = d["Label_T4"].value_counts().to_dict()
            rows.append({"seed": s, "fold": fold,
                         "w0": int(vc.get(0, 0)), "w1": int(vc.get(1, 0)),
                         "w2": int(vc.get(2, 0)), "n": len(d)})
    return pd.DataFrame(rows)


def render(df, out_path):
    headers = ["Seed", "Fold", WINDOW[0], WINDOW[1], WINDOW[2], "n"]
    LEFT_COLS = {0, 1}
    body, rule_above, show_seed = [], [], []
    prev_s = None
    for _, r in df.iterrows():
        body.append([f"seed {r.seed}" if r.seed != prev_s else "", r.fold,
                     str(r.w0), str(r.w1), str(r.w2), str(r.n)])
        rule_above.append(prev_s is not None and r.seed != prev_s)
        show_seed.append(r.seed != prev_s)
        prev_s = r.seed

    COL_W = [1.05, 0.95, 1.30, 1.70, 1.30, 0.65]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H = 1.18, 0.38, 0.42
    TOP_PAD, BOT_PAD = 0.12, 0.12
    n_rows = len(body)
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
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H * 0.76,
            "T4 Converter Cohort: Stratified Splits by Seed & Fold",
            ha="center", va="center", fontsize=11.5, fontweight="bold")
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H * 0.30,
            "80/10/10 stratified on the conversion window (Label_T4)\n"
            "146 converters (pMCI + pCN_to_AD)\n"
            "per-seed totals: < 3 yrs = 63  ·  ≥ 3, < 7 yrs = 47  ·  ≥ 7 yrs = 36",
            ha="center", va="center", fontsize=8.3, fontstyle="italic", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=9.5, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rule_above[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                w = "bold" if (j == 0 and show_seed[i]) else "normal"
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center",
                        fontsize=9.0, fontweight=w)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=9.0)
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def main():
    df = collect()
    df.rename(columns={"w0": WINDOW[0], "w1": WINDOW[1], "w2": WINDOW[2]}).to_csv(
        OUT_DIR / "t4_stratification_splits.csv", index=False)
    render(df, OUT_DIR / "t4_stratification_splits.png")
    print(f"  CSV: {OUT_DIR / 't4_stratification_splits.csv'}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
