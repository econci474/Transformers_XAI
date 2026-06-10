"""
03_summary_xmodal.py
====================
Aggregate the cross-modal fusion grid (outputs/**/metrics.json) into a house-style summary
table (PNG + PDF + .tex + .csv). One body row per (clin_arm, MRI variant, architecture, loss),
mean +/- std over seeds 0/1/2; clinical-only references (same 478 cohort) shown below the fusion
block. Bold = best per metric column. The prior late-fusion best (0.855, n~48) and MRI-only m12
are noted in the caption for context (different cohort).

Run:  python integration/cross_modal_attention/03_summary_xmodal.py
Out:  outputs/summary_xmodal.{csv,tex,png,pdf}
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
LOSS_DISP = {"ce": "CE", "ce_patient": "CE + patient-SigLIP",
             "ce_patient_label": "CE + patient + label-SigLIP"}
ARCH_DISP = {"mlp_concat": "MLP-concat", "cross_attn": "Cross-attn", "cross_attn_x": "X-attn (CLS)"}
MRI_DISP = {"A": "BrainDINO-T2", "B": "BrainMVP-T1b"}
METRICS = [("balanced_acc", "Test bACC"), ("auc_roc_ovr", "macro-AUC"),
           ("macro_f1", "macro-F1")]


def collect():
    rows = []
    for mf in glob.glob(os.path.join(OUT, "**", "metrics.json"), recursive=True):
        b = json.load(open(mf)); c = b["config"]; t = b["test_metrics"]; v = b["val_metrics"]
        rows.append({
            "clin_arm": c["clin_arm"], "clin_only": c.get("clin_only", False),
            "variant": c["variant"], "arch": c["arch"], "loss": c["loss"], "seed": c["seed"],
            "clin_pool": c.get("clin_pool", "mean"), "n_test": c["n_test"],
            "test_balanced_acc": t["balanced_acc"], "test_auc_roc_ovr": t["auc_roc_ovr"],
            "test_macro_f1": t["macro_f1"], "val_balanced_acc": v["balanced_acc"],
        })
    return pd.DataFrame(rows)


def aggregate(df):
    keys = ["clin_arm", "clin_only", "variant", "arch", "loss", "clin_pool"]
    g = df.groupby(keys)
    agg = g.agg(
        bacc_m=("test_balanced_acc", "mean"), bacc_s=("test_balanced_acc", "std"),
        auc_m=("test_auc_roc_ovr", "mean"), auc_s=("test_auc_roc_ovr", "std"),
        f1_m=("test_macro_f1", "mean"), f1_s=("test_macro_f1", "std"),
        valbacc_m=("val_balanced_acc", "mean"), n_test=("n_test", "mean"),
        n_seeds=("seed", "nunique"),
    ).reset_index()
    return agg


def _order(agg):
    """Fusion rows first (BL then LONG, by descending test bACC), then clinical-only refs."""
    fus = agg[~agg.clin_only].copy().sort_values(["clin_arm", "bacc_m"], ascending=[True, False])
    ref = agg[agg.clin_only].copy().sort_values(["clin_arm", "variant", "arch"])
    fus["block"] = "fusion"; ref["block"] = "clin"
    return pd.concat([fus, ref], ignore_index=True)


def _arch_disp(r):
    a = ARCH_DISP[r["arch"]]
    if r["arch"] == "cross_attn_x" and r["clin_arm"] == "LONG":
        a += f" · {r['clin_pool']}"          # distinguish mean vs frozen-Mamba pooling
    return a


def _label(r):
    if r["clin_only"]:
        return (f"Clinical-only ({r['clin_arm']})", "—", _arch_disp(r), "—")
    return (f"Fusion ({r['clin_arm']})", MRI_DISP[r["variant"]], _arch_disp(r),
            LOSS_DISP[r["loss"]])


def render_png(tab, out_png):
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    LEAD = ["Method", "MRI", "Architecture", "Loss"]
    headers = LEAD + [lbl for _, lbl in METRICS] + ["Val bACC", "n"]
    body, numeric, rules = [], [], []
    prev = None
    for _, r in tab.iterrows():
        cells = list(_label(r))
        nums = []
        for k in ("bacc", "auc", "f1"):
            m, s = r[f"{k}_m"], r[f"{k}_s"]
            cells.append(f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}")
            nums.append(float(m))
        cells.append(f"{r['valbacc_m']:.3f}")
        cells.append(f"{r['n_test']:.0f}")
        body.append(cells); numeric.append(nums)
        rules.append(prev is not None and r["block"] != prev); prev = r["block"]
    numeric = np.array(numeric, float)
    best = [int(np.nanargmax(numeric[:, j])) for j in range(len(METRICS))]

    COL_W = [3.05, 2.35, 2.05, 3.15, 1.70, 1.70, 1.70, 1.30, 0.55]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H, TOP_PAD, BOT_PAD = 1.30, 0.40, 0.44, 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    cl = [LEFT]
    for w in COL_W:
        cl.append(cl[-1] + w)
    RIGHT = cl[-1]
    ccx = [(cl[i] + cl[i + 1]) / 2 for i in range(len(COL_W))]
    LEFT_COLS = {0, 1, 2, 3}

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD; y_top = y; y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + TITLE_H / 2,
            "Cross-modal feature fusion — 12-month T2 diagnosis (CN/MCI/AD)\n"
            "mean ± std across seeds 0/1/2 · fusion ranked by Test bACC; clinical-only refs below\n"
            "Clinical embedding = BioClin-ModernBERT-L (T2, full_ft); MRI = m12 scan · common cohort n≈478",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    hline(y_top, lw=1.5); hline(y, lw=1.2)
    yh = y; y -= HEAD_H; ym = (yh + y) / 2
    for j, h in enumerate(headers):
        if j in LEFT_COLS:
            ax.text(cl[j] + 0.06, ym, h, ha="left", va="center", fontsize=9.5, fontstyle="italic")
        else:
            ax.text(ccx[j], ym, h, ha="center", va="center", fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)
    y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H; ym = (yt + y) / 2
        for j, txt in enumerate(cells):
            mj = j - len(LEAD)
            bold = (0 <= mj < len(METRICS) and best[mj] == i)
            if j in LEFT_COLS:
                ax.text(cl[j] + 0.06, ym, txt, ha="left", va="center", fontsize=9.0)
            else:
                ax.text(ccx[j], ym, txt, ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bold else "normal")
    BOTTOM = y
    for x in cl[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    fig.savefig(out_png.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_png}")


def write_tex(tab, out_tex):
    def esc(s):
        return str(s).replace("—", "--").replace("±", r"$\pm$").replace("_", r"\_")
    L = [r"% Cross-modal feature fusion — 12-month T2 (regenerated by 03_summary_xmodal.py)",
         r"\begin{table}[ht]", r"\centering", r"\small", r"\begin{tabular}{llll ccc c c}",
         r"\toprule",
         r"\textbf{Method} & \textbf{MRI} & \textbf{Architecture} & \textbf{Loss} & "
         r"\textbf{Test bACC} & \textbf{macro-AUC} & \textbf{macro-F1} & \textbf{Val bACC} & "
         r"\textbf{n} \\", r"\midrule"]
    prev = None
    for _, r in tab.iterrows():
        if prev is not None and r["block"] != prev:
            L.append(r"\midrule")
        prev = r["block"]
        meth, mri, arch, loss = _label(r)
        cells = [esc(meth), esc(mri), esc(arch), esc(loss),
                 f"${r['bacc_m']:.3f} \\pm {r['bacc_s']:.3f}$" if pd.notna(r['bacc_s']) else f"${r['bacc_m']:.3f}$",
                 f"${r['auc_m']:.3f} \\pm {r['auc_s']:.3f}$" if pd.notna(r['auc_s']) else f"${r['auc_m']:.3f}$",
                 f"${r['f1_m']:.3f} \\pm {r['f1_s']:.3f}$" if pd.notna(r['f1_s']) else f"${r['f1_m']:.3f}$",
                 f"{r['valbacc_m']:.3f}", f"{r['n_test']:.0f}"]
        L.append(" & ".join(cells) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\caption{Cross-modal \emph{feature} fusion for the 12-month three-way diagnosis "
          r"(CN/MCI/AD). Clinical embedding = BioClinical-ModernBERT-L (T2, full\_ft); MRI = m12 "
          r"scan. BL = baseline clinical token; LONG = bl+m06+m12 clinical tokens (single "
          r"T2\_long model). Common cohort n$\approx$478 across all rows. Auxiliary SigLIP losses "
          r"at $\lambda=\gamma=0.01$. Mean $\pm$ std over seeds 0/1/2. For reference, the prior "
          r"\emph{late} (probability) fusion reached 0.855 test bACC (n$\approx$48 cohort).}",
          r"\label{tab:xmodal_t2_m12}", r"\end{table}"]
    open(out_tex, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"  TEX: {out_tex}")


def main():
    df = collect()
    if df.empty:
        print("No metrics.json found — run 02_run_grid.py first."); return
    agg = aggregate(df)
    tab = _order(agg)
    tab.to_csv(os.path.join(OUT, "summary_xmodal.csv"), index=False)
    print(f"  CSV: {os.path.join(OUT, 'summary_xmodal.csv')}  ({len(tab)} rows, "
          f"{df['seed'].nunique()} seeds)")
    render_png(tab, os.path.join(OUT, "summary_xmodal.png"))
    write_tex(tab, os.path.join(OUT, "summary_xmodal.tex"))
    # quick console leaderboard
    top = tab[tab.block == "fusion"].nlargest(5, "bacc_m")
    print("\n  Top-5 fusion by Test bACC:")
    for _, r in top.iterrows():
        print(f"    {r['clin_arm']:4s} {MRI_DISP[r['variant']]:13s} {ARCH_DISP[r['arch']]:11s} "
              f"{LOSS_DISP[r['loss']]:28s} bACC {r['bacc_m']:.3f}±{r['bacc_s']:.3f} "
              f"AUC {r['auc_m']:.3f} F1 {r['f1_m']:.3f}")
    ref = tab[tab.block == "clin"]
    print("  Clinical-only refs (same cohort):")
    for _, r in ref.iterrows():
        print(f"    {r['clin_arm']:4s} {ARCH_DISP[r['arch']]:11s} (var {r['variant']}) "
              f"bACC {r['bacc_m']:.3f}±{r['bacc_s']:.3f} AUC {r['auc_m']:.3f} F1 {r['f1_m']:.3f}")


if __name__ == "__main__":
    main()
