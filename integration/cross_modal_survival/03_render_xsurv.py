r"""
03_render_xsurv.py  (env: snp — matplotlib/pandas)
==================================================
House-style cleaned table of the cross-modal survival comparison (val+test pooled, mean ± std over 3
seeds, bold best per column). Rows: Fused (BrainMVP-T1d / BrainDINO-T1d × alignment arm) then the
Clinical-only and MRI-only references. Reads outputs/xsurv/master_xsurv_comparison.csv (no re-fit).

Run:  conda run -n snp python integration/cross_modal_survival/03_render_xsurv.py
Out:  outputs/xsurv/survival_xmodal_table.{csv,png,pdf}
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = Path(__file__).resolve().parent
MASTER = HERE / "outputs" / "xsurv" / "master_xsurv_comparison.csv"
OUT = HERE / "outputs" / "xsurv"
MRI_NAME = {"mvp": "BrainMVP-T1d", "dino": "BrainDINO-T1d"}
ALIGN_LABEL = {"none": "none", "joint": "joint (SigLIP)", "prealign": "pre-align (SigLIP)"}

# (model key, group, Modality, Imaging, Alignment)
ROWS_SPEC = (
    [(f"{v}__{a}", f"fused_{v}", "Fused", MRI_NAME[v], ALIGN_LABEL[a])
     for v in ("mvp", "dino") for a in ("none", "joint", "prealign")]
    + [("mvp__clinical", "ref", "Clinical-only", "—", "—"),
       ("mvp__mri", "ref", "MRI-only", MRI_NAME["mvp"], "—"),
       ("dino__mri", "ref", "MRI-only", MRI_NAME["dino"], "—")]
)
METRICS = [("c_index", "C-index", "max"), ("c_index_ipcw", "IPCW C", "max"), ("ibs", "IBS", "min"),
           ("auc_3y", "AUC@3y", "max"), ("auc_5y", "AUC@5y", "max"),
           ("auc_7y", "AUC@7y", "max"), ("auc_10y", "AUC@10y", "max")]
TITLE = "Cross-modal survival — bidirectional cross-attention (T2-long clinical ⊕ T1d MRI) → Weibull-piecewise (CN+MCI → AD)"
SUBTITLE = (
    "held-out = val+test pooled  ·  mean ± std over 3 seeds  ·  bold = best per column  ·  risk = 1 − S(t_ref=4y)\n"
    "Fused = clinical CLS (mean-pooled bl/m06/m12) ↔ T1d-MRI (m12) bidirectional cross-attention → concat → Weibull-piecewise head\n"
    "Alignment = patient-SigLIP (joint = throughout; pre-align = SigLIP-only warmup, then survival NLL)  ·  references on the SAME cohort")


def aggregate():
    d = pd.read_csv(MASTER)
    d = d[d["split"] == "valtest"]
    recs = []
    for key, group, modality, imaging, align in ROWS_SPEC:
        g = d[d["model"] == key]
        if g.empty:
            print(f"  [warn] missing valtest rows for {key} — skipping")
            continue
        rec = {"key": key, "group": group, "Modality": modality, "Imaging": imaging,
               "Alignment": align, "n": g["n"].mean(), "events": g["events"].mean()}
        for m, _, _ in METRICS:
            rec[f"{m}_m"] = g[m].mean(); rec[f"{m}_s"] = g[m].std(ddof=1)
        recs.append(rec)
    return pd.DataFrame(recs).reset_index(drop=True)


def render(df, out_path):
    headers = ["Modality", "Imaging", "Alignment", "N (events)"] + [lbl for _, lbl, _ in METRICS]
    LEFT_COLS = {0, 1, 2}
    best = {}
    for m, _, mode in METRICS:
        col = df[f"{m}_m"].to_numpy(float)
        best[m] = int(np.nanargmin(col) if mode == "min" else np.nanargmax(col))

    body, rules, bolds = [], [], []
    prev_g = None
    for i, r in df.iterrows():
        if prev_g is None:
            rule = None
        elif r.group == "ref" and prev_g != "ref":
            rule = "strong"
        elif r.group != prev_g:
            rule = "light"
        else:
            rule = None
        cells = [r.Modality, r.Imaging, r.Alignment, f"{r['n']:.0f} ({r['events']:.0f})"]
        for m, _, _ in METRICS:
            cells.append(f"{r[f'{m}_m']:.3f} ± {r[f'{m}_s']:.3f}")
        body.append(cells)
        bolds.append({4 + j: (best[m] == i) for j, (m, _, _) in enumerate(METRICS)})
        rules.append(rule); prev_g = r.group

    COL_W = [1.55, 1.65, 1.85, 1.20] + [1.62] * len(METRICS)
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_LINES = SUBTITLE.count("\n") + 1
    SUB_H = 0.28 * SUB_LINES
    TITLE_H, HEAD_H, ROW_H = 0.55 + SUB_H, 0.42, 0.46
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

    y = fig_h - TOP_PAD; y_title_top = y; y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) * 0.55, TITLE,
            ha="center", va="center", fontsize=11.5, fontweight="bold")
    ax.text(cx, y + SUB_H / 2, SUBTITLE, ha="center", va="center",
            fontsize=7.6, fontstyle="italic", linespacing=1.5)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

    y_head_top = y; y -= HEAD_H; ymid = (y_head_top + y) / 2
    for j in range(len(headers)):
        if j in LEFT_COLS:
            ax.text(col_left[j] + 0.06, ymid, headers[j], ha="left", va="center",
                    fontsize=9.0, fontstyle="italic")
        else:
            ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                    fontsize=9.0, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i] == "strong":
            hline(y, lw=0.9, ls=(0, (4, 3)))
        elif rules[i] == "light":
            hline(y, lw=0.5, ls=(0, (1, 3)))
        yt = y; y -= ROW_H; ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center", fontsize=8.8)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=8.8,
                        fontweight="bold" if bolds[i].get(j, False) else "normal")
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
    df = aggregate()
    out_cols = ["Modality", "Imaging", "Alignment", "n", "events"] \
        + [f"{m}_m" for m, _, _ in METRICS] + [f"{m}_s" for m, _, _ in METRICS]
    df[out_cols].to_csv(OUT / "survival_xmodal_table.csv", index=False, encoding="utf-8")
    render(df, OUT / "survival_xmodal_table.png")
    print(f"  CSV: {OUT / 'survival_xmodal_table.csv'}  ({len(df)} rows)")
    print(df[["Modality", "Imaging", "Alignment", "c_index_m", "c_index_ipcw_m", "ibs_m"]].to_string(index=False))


if __name__ == "__main__":
    main()
