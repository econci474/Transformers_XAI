"""
_render_int_t2_simple_filtered.py
=================================
Render a house-style PNG/PDF twin of `outputs/int_t2_simple_filtered.tex` (the hand-curated
"best fusion method vs. its single-modality references" T2 table). Standalone — does NOT import
the pipeline module `04_summary_top_methods.py`, so it cannot perturb that pipeline's behaviour;
it only mirrors the visual conventions of its `_draw` house-style renderer (DejaVu Serif, bordered,
ruled, italic headers, dashed block dividers, bold best-per-column).

Rows are hardcoded to match the .tex exactly, including the two newly-added BrainMVP-ft (T1a) /
(T1b) MRI-only binary-detector references (genuine two-class macro-F1 derived from
brainmvp_debug/aug_stochastic per-class F1: T1a 0.55, T1b 0.68). Bolding mirrors the .tex
(bACC best = row 0; macro-F1 best = row 1).

Writes: outputs/int_t2_simple_filtered.png (+ .pdf).
Run (env with matplotlib, e.g. snp):  python integration/T2_classification/_render_int_t2_simple_filtered.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "DejaVu Serif"

OUT = Path(__file__).resolve().parent / "outputs"

# block, Method, Clinical (Task), MRI (Task), Aggregation, bACC, bACC_std, macro-AUC, macro-F1
ROWS = [
    ("fusion", "CL at baseline + MRI at 12 months", "BioClin-MBERT-L-ft (T2)",
     "BrainDINO-frozen (T2)", "weighted-avg (val bACC)", 0.86, 0.08, 0.91, 0.75),
    ("fusion", "CL longitudinal + MRI at 12 months",
     "BioClin-MBERT-B-ft (T1a) + ModernBERT-B-ft (T1b)", "BrainMVP-ft (T1a/T1b)",
     "class-wise, mean", 0.84, 0.04, 0.95, 0.84),
    ("fusion", "CL longitudinal + MRI at 12 months", "BioClin-MBERT-L-ft (T2)",
     "BrainMVP-ft (T1a/T1b)", "class-wise, mean", 0.82, 0.03, 0.95, 0.82),
    ("clin", "CL longitudinal", "BioClin-MBERT-B-ft (T1a) + ModernBERT-B-ft (T1b)",
     "—", "class-wise, mean", 0.84, 0.04, 0.95, 0.84),
    ("clin", "CL at baseline", "BioClin-MBERT-L-ft (T2)", "—", "—", 0.84, 0.06, 0.93, 0.74),
    ("mri", "MRI at 12 months", "—", "BrainDINO-frozen (T2)", "—", 0.57, 0.05, 0.67, 0.39),
    ("mri", "MRI at 12 months", "—", "BrainMVP-ft (T2)", "—", 0.45, 0.08, 0.64, 0.44),
    ("mri", "MRI at 12 months", "—", "BrainMVP-ft (T1a)", "—", 0.54, 0.06, 0.58, 0.52),
    ("mri", "MRI at 12 months", "—", "BrainMVP-ft (T1b)", "—", 0.74, 0.06, 0.81, 0.76),
]

# Bold cells, mirroring the .tex: (row_index, metric_index) where metric 0=bACC, 1=macro-AUC, 2=macro-F1
BOLD = {(0, 0), (1, 2)}

BLOCK_LABEL = {"fusion": "Clinical ⊕ MRI late fusion",
               "clin": "Clinical-only references", "mri": "MRI-only references"}
LEAD = ["Method", "Clinical (Task)", "MRI (Task)", "Aggregation"]
METRIC_HEAD = ["Test bACC", "macro-AUC", "macro-F1"]


def main():
    headers = LEAD + METRIC_HEAD
    N_LEAD = len(LEAD)

    body, rules = [], []
    prev_b = None
    for r in ROWS:
        block, method, clin, mri, agg, bacc, bstd, auc, f1 = r
        cells = [method, clin, mri, agg,
                 f"{bacc:.2f} ± {bstd:.2f}", f"{auc:.2f}", f"{f1:.2f}"]
        body.append(cells)
        rules.append(prev_b is not None and block != prev_b)
        prev_b = block

    COL_W = [4.05, 5.55, 2.95, 2.95, 1.55, 1.20, 1.15]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H = 1.10, 0.40, 0.44
    TOP_PAD, BOT_PAD = 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
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
    y_title_top = y; y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + TITLE_H / 2,
            "T2 Late Fusion — best method vs. its single-modality references\n"
            "mean ± std across seeds 0/1/2 · fusion on top, references below\n"
            "BrainMVP-ft (T1a)/(T1b) = standalone binary detectors on the same 12-month cohort "
            "(T1a CN vs MCI+AD n≈48; T1b CN vs AD n≈23); macro-F1 = two-class macro-averaged F1",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    LEFT_COLS = {0, 1, 2, 3}
    y_head_top = y; y -= HEAD_H
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
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        # block sub-header label sits at the left of the first row of each block
        yr_top = y; y -= ROW_H
        ymid = (yr_top + y) / 2
        for j in range(len(headers)):
            metric_j = j - N_LEAD
            bold = (metric_j >= 0 and (i, metric_j) in BOLD)
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ymid, cells[j], ha="left", va="center", fontsize=9.0)
            else:
                ax.text(col_cx[j], ymid, cells[j], ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bold else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))

    out_png = OUT / "int_t2_simple_filtered.png"
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_png}")
    print(f"  PDF: {out_png.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
