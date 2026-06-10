"""
03b_summary_lambda_sweep.py
===========================
Summarise the LONG-arm contrastive-weight sweep (warmup=10, min_epochs=30) against the
lambda=0.01 main-grid baseline and the clinical-only LONG reference. Answers: does giving the
SigLIP losses a fair chance (warmup + larger lambda) ever beat clinical-only LONG (~0.832)?

Rows: LONG fusion (variant x arch x loss) at lambda in {0.01, 0.1, 0.5, 1.0}. Reference rows:
clinical-only LONG and LONG fusion with CE only (no aux). Mean +/- std over seeds 0/1/2.

Run:  python integration/cross_modal_attention/03b_summary_lambda_sweep.py
Out:  outputs/summary_lambda_sweep.{csv,tex,png,pdf}
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
ARCH_DISP = {"mlp_concat": "MLP-concat", "cross_attn": "Cross-attn", "cross_attn_x": "X-attn (CLS)"}
MRI_DISP = {"A": "BrainDINO-T2", "B": "BrainMVP-T1b"}
LOSS_DISP = {"ce_patient": "CE + patient", "ce_patient_label": "CE + patient + label"}


def _read(mf):
    b = json.load(open(mf)); c = b["config"]; t = b["test_metrics"]; v = b["val_metrics"]
    return {"clin_arm": c["clin_arm"], "clin_only": c.get("clin_only", False),
            "variant": c["variant"], "arch": c["arch"], "loss": c["loss"], "seed": c["seed"],
            "clin_pool": c.get("clin_pool", "mean"),
            "lam": c.get("lambda", 0.0), "warmup": c.get("warmup", 0), "n_test": c["n_test"],
            "bacc": t["balanced_acc"], "auc": t["auc_roc_ovr"], "f1": t["macro_f1"],
            "valbacc": v["balanced_acc"]}


def collect():
    rows = []
    # sweep (lambda 0.1/0.5/1.0, warmup 10)
    rows += [_read(m) for m in glob.glob(os.path.join(OUT, "sweep", "**", "metrics.json"),
                                         recursive=True)]
    # lambda=0.01 baseline + CE + clinical-only, all LONG, from the main grid (exclude sweep dir)
    for m in glob.glob(os.path.join(OUT, "LONG*", "**", "metrics.json"), recursive=True):
        if "sweep" in m:
            continue
        rows.append(_read(m))
    return pd.DataFrame(rows)


def agg(df, keys):
    g = df.groupby(keys)
    return g.agg(bacc_m=("bacc", "mean"), bacc_s=("bacc", "std"),
                 auc_m=("auc", "mean"), auc_s=("auc", "std"),
                 f1_m=("f1", "mean"), f1_s=("f1", "std"),
                 valbacc_m=("valbacc", "mean"), n_test=("n_test", "mean")).reset_index()


def build_table(df):
    fus = df[(df.clin_arm == "LONG") & (~df.clin_only) &
             (df.loss.isin(["ce_patient", "ce_patient_label"]))].copy()
    a = agg(fus, ["variant", "arch", "clin_pool", "loss", "lam"])
    a["block"] = "sweep"
    a = a.sort_values(["loss", "arch", "variant", "lam"]).reset_index(drop=True)

    # references
    refs = []
    co = df[(df.clin_arm == "LONG") & (df.clin_only)]
    if not co.empty:
        r = agg(co, ["arch"]); r["block"] = "ref"; r["loss"] = "clinical-only"; r["variant"] = "—"; r["lam"] = np.nan
        refs.append(r)
    ceonly = df[(df.clin_arm == "LONG") & (~df.clin_only) & (df.loss == "ce")]
    if not ceonly.empty:
        r = agg(ceonly, ["variant", "arch", "clin_pool"]); r["block"] = "ref"
        r["loss"] = "fusion-CE"; r["lam"] = 0.0
        refs.append(r)
    out = pd.concat([a] + refs, ignore_index=True)
    return out


def _arch_disp(r):
    a = ARCH_DISP[r["arch"]]
    if r["arch"] == "cross_attn_x":
        a += f" · {r.get('clin_pool', 'mean')}"   # LONG: mean vs frozen-Mamba pooling
    return a


def _label(r):
    if r["block"] == "ref" and r["loss"] == "clinical-only":
        return ("Clinical-only LONG", "—", ARCH_DISP[r["arch"]], "—", "—")
    if r["block"] == "ref" and r["loss"] == "fusion-CE":
        return ("Fusion LONG (CE only)", MRI_DISP[r["variant"]], _arch_disp(r), "CE", "—")
    return ("Fusion LONG", MRI_DISP[r["variant"]], _arch_disp(r),
            LOSS_DISP[r["loss"]], f"{r['lam']:g}")


def render_png(tab, out_png):
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    LEAD = ["Method", "MRI", "Architecture", "Loss", "λ=γ"]
    headers = LEAD + ["Test bACC", "macro-AUC", "macro-F1", "Val bACC", "n"]
    body, numeric, rules = [], [], []
    prev = None
    for _, r in tab.iterrows():
        cells = list(_label(r))
        nums = []
        for k in ("bacc", "auc", "f1"):
            m, s = r[f"{k}_m"], r[f"{k}_s"]
            cells.append(f"{m:.3f} ± {s:.3f}" if pd.notna(s) else f"{m:.3f}")
            nums.append(float(m))
        cells.append(f"{r['valbacc_m']:.3f}"); cells.append(f"{r['n_test']:.0f}")
        body.append(cells); numeric.append(nums)
        rules.append(prev is not None and r["block"] != prev); prev = r["block"]
    numeric = np.array(numeric, float)
    best = [int(np.nanargmax(numeric[:, j])) for j in range(3)]

    COL_W = [2.95, 2.30, 2.00, 2.85, 0.90, 1.70, 1.70, 1.70, 1.30, 0.55]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H, TOP_PAD, BOT_PAD = 1.30, 0.40, 0.44, 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    cl = [LEFT]
    for w in COL_W:
        cl.append(cl[-1] + w)
    RIGHT = cl[-1]; ccx = [(cl[i] + cl[i + 1]) / 2 for i in range(len(COL_W))]
    LEFT_COLS = {0, 1, 2, 3}

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD; y_top = y; y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H / 2,
            "Contrastive-weight sweep — LONG arm, 12-month T2 (CN/MCI/AD)\n"
            "warmup=10 (aux-only) + min_epochs=30 · mean ± std seeds 0/1/2 · ranked within block by Test bACC\n"
            "Does giving SigLIP a fair chance beat clinical-only LONG (0.832)?",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    hline(y_top, lw=1.5); hline(y, lw=1.2)
    yh = y; y -= HEAD_H; ym = (yh + y) / 2
    for j, h in enumerate(headers):
        if j in LEFT_COLS:
            ax.text(cl[j] + 0.06, ym, h, ha="left", va="center", fontsize=9.5, fontstyle="italic")
        else:
            ax.text(ccx[j], ym, h, ha="center", va="center", fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2); y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H; ym = (yt + y) / 2
        for j, txt in enumerate(cells):
            mj = j - len(LEAD)
            bold = (0 <= mj < 3 and best[mj] == i)
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
        return str(s).replace("—", "--").replace("λ", r"$\lambda$").replace("γ", r"$\gamma$")
    L = [r"% LONG-arm contrastive sweep (03b_summary_lambda_sweep.py)",
         r"\begin{table}[ht]", r"\centering", r"\small",
         r"\begin{tabular}{lllll ccc c c}", r"\toprule",
         r"\textbf{Method} & \textbf{MRI} & \textbf{Architecture} & \textbf{Loss} & "
         r"$\boldsymbol{\lambda{=}\gamma}$ & \textbf{Test bACC} & \textbf{macro-AUC} & "
         r"\textbf{macro-F1} & \textbf{Val bACC} & \textbf{n} \\", r"\midrule"]
    prev = None
    for _, r in tab.iterrows():
        if prev is not None and r["block"] != prev:
            L.append(r"\midrule")
        prev = r["block"]
        meth, mri, arch, loss, lam = _label(r)
        cells = [esc(meth), esc(mri), esc(arch), esc(loss), lam,
                 f"${r['bacc_m']:.3f} \\pm {r['bacc_s']:.3f}$" if pd.notna(r['bacc_s']) else f"${r['bacc_m']:.3f}$",
                 f"${r['auc_m']:.3f} \\pm {r['auc_s']:.3f}$" if pd.notna(r['auc_s']) else f"${r['auc_m']:.3f}$",
                 f"${r['f1_m']:.3f} \\pm {r['f1_s']:.3f}$" if pd.notna(r['f1_s']) else f"${r['f1_m']:.3f}$",
                 f"{r['valbacc_m']:.3f}", f"{r['n_test']:.0f}"]
        L.append(" & ".join(cells) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\caption{LONG-arm contrastive-weight sweep for the 12-month three-way diagnosis. "
          r"SigLIP given a fair chance via warmup=10 (aux-only epochs) + min\_epochs=30. "
          r"$\lambda{=}0.01$ is the main-grid baseline (no warmup). Clinical-only LONG and "
          r"fusion-with-CE-only shown as references. Mean $\pm$ std over seeds 0/1/2; n$\approx$478.}",
          r"\label{tab:xmodal_lambda_sweep}", r"\end{table}"]
    open(out_tex, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"  TEX: {out_tex}")


def main():
    df = collect()
    sweep_n = (df.warmup > 0).sum()
    if sweep_n == 0:
        print("No sweep runs found under outputs/sweep — run 02b first."); return
    tab = build_table(df)
    tab.to_csv(os.path.join(OUT, "summary_lambda_sweep.csv"), index=False)
    print(f"  CSV: summary_lambda_sweep.csv ({len(tab)} rows; {sweep_n} sweep runs)")
    render_png(tab, os.path.join(OUT, "summary_lambda_sweep.png"))
    write_tex(tab, os.path.join(OUT, "summary_lambda_sweep.tex"))
    best = tab[tab.block == "sweep"].nlargest(6, "bacc_m")
    print("\n  Top sweep rows by Test bACC (bar to beat: clinical-only LONG 0.832):")
    for _, r in best.iterrows():
        print(f"    {MRI_DISP[r['variant']]:13s} {ARCH_DISP[r['arch']]:11s} {LOSS_DISP[r['loss']]:22s} "
              f"λ={r['lam']:<4g} bACC {r['bacc_m']:.3f}±{r['bacc_s']:.3f} AUC {r['auc_m']:.3f}")


if __name__ == "__main__":
    main()
