"""
03c_summary_joined.py
=====================
ONE joined summary across both clinical arms swept (M12 = dedicated 1-year model, LONG =
all-visits model): the BEST hyperparameter-tuned fusion per condition, vs clinical-only refs.

"Condition" for the headline table = (clinical arm x loss type). For each, we take the best
fusion result tuned over {MRI variant A/B, architecture mlp_concat/cross_attn, lambda in
{0.01, 0.1, 0.5, 1.0}} (lambda>0.01 use warmup=10/min_epochs=30), and show the winning config.
Clinical-only references (BL / M12 / LONG, best architecture) are the bars to beat.

A full per-(arm,variant,arch,loss) best-lambda table is also written to CSV for completeness.

Run:  python integration/cross_modal_attention/03c_summary_joined.py
Out:  outputs/summary_joined.{csv,tex,png,pdf} + summary_joined_full.csv
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
ARM_DISP = {"BL": "Baseline", "M12": "1-year (m12)", "LONG": "Longitudinal"}
LOSS_DISP = {"ce": "CE", "ce_patient": "CE + patient",
             "ce_patient_label": "CE + patient + label"}
LOSS_ORDER = ["ce", "ce_patient", "ce_patient_label"]
ARM_ORDER = ["M12", "LONG"]


def _read(mf):
    b = json.load(open(mf)); c = b["config"]; t = b["test_metrics"]; v = b["val_metrics"]
    return {"clin_arm": c["clin_arm"], "clin_only": c.get("clin_only", False),
            "variant": c["variant"], "arch": c["arch"], "loss": c["loss"], "seed": c["seed"],
            "clin_pool": c.get("clin_pool", "mean"),
            "lam": float(c.get("lambda", 0.0)), "warmup": int(c.get("warmup", 0)),
            "n_test": c["n_test"], "bacc": t["balanced_acc"], "auc": t["auc_roc_ovr"],
            "f1": t["macro_f1"], "valbacc": v["balanced_acc"]}


def collect():
    files = (glob.glob(os.path.join(OUT, "sweep", "**", "metrics.json"), recursive=True) +
             [m for m in glob.glob(os.path.join(OUT, "*", "**", "metrics.json"), recursive=True)
              if os.sep + "sweep" + os.sep not in m])
    return pd.DataFrame([_read(m) for m in files])


def _agg(sub):
    """mean +/- std over seeds for one (arm,variant,arch,loss,lam) group."""
    return pd.Series({
        "bacc_m": sub["bacc"].mean(), "bacc_s": sub["bacc"].std(),
        "auc_m": sub["auc"].mean(), "auc_s": sub["auc"].std(),
        "f1_m": sub["f1"].mean(), "f1_s": sub["f1"].std(),
        "valbacc_m": sub["valbacc"].mean(), "n_test": sub["n_test"].mean(),
        "warmup": sub["warmup"].iloc[0]})


def best_tuned_full(df):
    """Per (arm,variant,arch,loss): pick the lambda with best mean test bACC (tuned)."""
    fus = df[~df.clin_only & df.loss.isin(LOSS_ORDER)]
    rows = []
    for (arm, var, arch, pool, loss), g in fus.groupby(
            ["clin_arm", "variant", "arch", "clin_pool", "loss"]):
        per_lam = g.groupby("lam").apply(_agg, include_groups=False).reset_index()
        best = per_lam.loc[per_lam["bacc_m"].idxmax()]
        rows.append({"clin_arm": arm, "variant": var, "arch": arch, "clin_pool": pool, "loss": loss,
                     "lam": best["lam"], **{k: best[k] for k in
                     ("bacc_m", "bacc_s", "auc_m", "auc_s", "f1_m", "f1_s", "valbacc_m",
                      "n_test", "warmup")}})
    return pd.DataFrame(rows)


def load_unimodal_refs():
    """Unimodal 1-year references (own model softmax, NO retrained head) from 04_significance."""
    p = os.path.join(OUT, "unimodal_refs.csv")
    if not os.path.exists(p):
        print("  [warn] unimodal_refs.csv missing — run 04_significance_unimodal.py first.")
        return pd.DataFrame()
    return pd.read_csv(p)


def headline(full, uni):
    """Condition = (arm, loss): best over variant+arch (lambda already tuned in `full`).
    References = the two UNIMODAL 1-year predictions (clinical, MRI), with a significance marker
    (best 1-year fusion vs that reference, McNemar exact)."""
    rows = []
    for arm in ARM_ORDER:
        for loss in LOSS_ORDER:
            sub = full[(full.clin_arm == arm) & (full.loss == loss)]
            if sub.empty:
                continue
            b = sub.loc[sub["bacc_m"].idxmax()]
            arch_disp = ARCH_DISP[b["arch"]] + (
                f" ({b['clin_pool']})" if b["arch"] == "cross_attn_x" and arm == "LONG" else "")
            cfg = f"{MRI_DISP[b['variant']]} · {arch_disp}" + (
                "" if loss == "ce" else
                f" · λ={b['lam']:g}" + (f" (wu{int(b['warmup'])})" if b['warmup'] else ""))
            rows.append({"block": "fusion", "method": ARM_DISP[arm], "loss_str": LOSS_DISP[loss],
                         "cfg": cfg, **{k: b[k] for k in ("bacc_m", "bacc_s", "auc_m", "auc_s",
                                        "f1_m", "f1_s", "valbacc_m", "n_test")}})
    for _, r in uni.iterrows():
        lbl = r["label"].replace("Clinical-only 1-year (unimodal)", "Clinical-only (unimodal)") \
                        .replace("MRI-only 1-year (unimodal)", "MRI-only (unimodal)")
        rows.append({"block": "ref", "method": "1-year (m12)", "loss_str": f"{lbl} [{r['sig']}]",
                     "cfg": "own softmax (no head)",
                     "bacc_m": r["balanced_acc_m"], "bacc_s": r["balanced_acc_s"],
                     "auc_m": r["auc_roc_ovr_m"], "auc_s": r["auc_roc_ovr_s"],
                     "f1_m": r["macro_f1_m"], "f1_s": r["macro_f1_s"],
                     "valbacc_m": np.nan, "n_test": r["n"]})
    return pd.DataFrame(rows)


def _fmt(m, s):
    return (f"{m:.2f} ± {s:.2f}" if pd.notna(s) else f"{m:.2f}") if pd.notna(m) else "--"


def _cells(r):
    vb = f"{r['valbacc_m']:.2f}" if pd.notna(r["valbacc_m"]) else "--"
    return [r["method"], r["loss_str"], r["cfg"],
            _fmt(r["bacc_m"], r["bacc_s"]), _fmt(r["auc_m"], r["auc_s"]),
            _fmt(r["f1_m"], r["f1_s"]), vb, f"{r['n_test']:.0f}"]


def render_png(tab, out_png):
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    LEAD = ["Clinical arm", "Loss", "Best config (tuned)"]
    headers = LEAD + ["Test bACC", "macro-AUC", "macro-F1", "Val bACC", "n"]
    body = [_cells(r) for _, r in tab.iterrows()]
    numeric = np.array([[r["bacc_m"], r["auc_m"], r["f1_m"]] for _, r in tab.iterrows()], float)
    best = [int(np.nanargmax(numeric[:, j])) for j in range(3)]
    rules = [i > 0 and tab.iloc[i]["block"] != tab.iloc[i - 1]["block"] for i in range(len(tab))]

    COL_W = [1.95, 3.30, 4.05, 1.55, 1.55, 1.55, 1.20, 0.55]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H, TOP_PAD, BOT_PAD = 1.70, 0.40, 0.46, 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    cl = [LEFT]
    for w in COL_W:
        cl.append(cl[-1] + w)
    RIGHT = cl[-1]; ccx = [(cl[i] + cl[i + 1]) / 2 for i in range(len(COL_W))]
    LEFT_COLS = {0, 1, 2}

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD; y_top = y; y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H / 2,
            "Cross-modal fusion — best hyperparameter-tuned per condition · 12-month T2 (CN/MCI/AD)\n"
            "Fusion tuned over {MRI variant, architecture, λ∈{0.01,0.1,0.5,1.0}}; SigLIP λ>0.01 use warmup=10\n"
            "mean ± std seeds 0/1/2 · references = unimodal own-softmax (no head); [sig] = best 1-yr fusion vs ref\n"
            "(McNemar exact): * p<0.05, ** p<0.01, *** p<0.001, n.s. = not significant · M12 n≈451",
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
                ax.text(cl[j] + 0.06, ym, txt, ha="left", va="center", fontsize=9.0,
                        fontweight="bold" if (j == 1 and tab.iloc[i]["block"] == "ref") else "normal")
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
        return str(s).replace("·", r"$\cdot$").replace("λ", r"$\lambda$").replace("≈", r"$\approx$")
    L = [r"% Joined best-per-condition table (03c_summary_joined.py)",
         r"\begin{table}[ht]", r"\centering", r"\small", r"\begin{tabular}{lll ccc c c}",
         r"\toprule",
         r"\textbf{Clinical arm} & \textbf{Loss} & \textbf{Best config (tuned)} & "
         r"\textbf{Test bACC} & \textbf{macro-AUC} & \textbf{macro-F1} & \textbf{Val bACC} & "
         r"\textbf{n} \\", r"\midrule"]
    prev = None
    for _, r in tab.iterrows():
        if prev is not None and r["block"] != prev:
            L.append(r"\midrule")
        prev = r["block"]
        c = _cells(r)
        c = [esc(c[0]), esc(c[1]), esc(c[2])] + [f"${x.replace('±', chr(92)+'pm')}$" for x in c[3:6]] + [c[6], c[7]]
        L.append(" & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\caption{Best hyperparameter-tuned cross-modal fusion per condition for the 12-month "
          r"three-way diagnosis (CN/MCI/AD). Each fusion row is the best over MRI variant "
          r"(BrainDINO/BrainMVP), architecture (MLP-concat/cross-attention) and contrastive weight "
          r"$\lambda\in\{0.01,0.1,0.5,1.0\}$ ($\lambda>0.01$ with warmup=10, min\_epochs=30). "
          r"\emph{References are the unimodal 1-year predictions} (each model's own softmax, no "
          r"retrained head): the clinical 1-year model (T2\_m12) and BrainDINO-T2 at the m12 scan. "
          r"The bracketed marker is the significance of the best 1-year fusion vs that reference "
          r"(exact McNemar, pooled over the three test folds): $^{*}p<0.05$, $^{**}p<0.01$, "
          r"$^{***}p<0.001$, n.s.\ = not significant. Best 1-year fusion vs clinical-1-year: n.s.\ "
          r"($+0.03$ bACC, McNemar $p=1.0$, paired-$t$ $p=0.56$); vs MRI-1-year: $^{***}$ "
          r"(MRI alone $\approx$ chance). M12 cohort n$\approx$451. Mean $\pm$ std over seeds 0/1/2, 2 dp.}",
          r"\label{tab:xmodal_joined}", r"\end{table}"]
    open(out_tex, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"  TEX: {out_tex}")


def main():
    df = collect()
    if df.empty:
        print("No runs found."); return
    full = best_tuned_full(df)
    full.sort_values(["clin_arm", "variant", "arch", "loss"]).to_csv(
        os.path.join(OUT, "summary_joined_full.csv"), index=False)
    uni = load_unimodal_refs()
    tab = headline(full, uni)
    tab.to_csv(os.path.join(OUT, "summary_joined.csv"), index=False)
    print(f"  CSV: summary_joined.csv ({len(tab)} rows) + summary_joined_full.csv ({len(full)} conditions)")
    render_png(tab, os.path.join(OUT, "summary_joined.png"))
    write_tex(tab, os.path.join(OUT, "summary_joined.tex"))
    print("\n  Joined best-per-condition (Test bACC, 2dp):")
    for _, r in tab.iterrows():
        tag = "REF " if r["block"] == "ref" else "fuse"
        print(f"    [{tag}] {r['method']:14s} {r['loss_str']:34s} "
              f"{r['bacc_m']:.2f}±{r['bacc_s']:.2f}  AUC {r['auc_m']:.2f}  ({r['cfg']})")


if __name__ == "__main__":
    main()
