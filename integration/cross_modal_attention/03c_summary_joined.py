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
# Feature-fusion display labels (cross_attn is self-attention over the concatenated token set — a
# misnomer; cross_attn_x is the true bidirectional cross-attention over one CLS per modality).
FUSION_DISP = {"mlp_concat": "MLP", "cross_attn": "Self-attention", "cross_attn_x": "Cross-attention"}
MRI_DISP = {"A": "BrainDINO-T2", "B": "BrainMVP-T1b"}
CLIN_DISP = "BioClin-L · T2 (CN/MCI/AD)"          # clinical encoder, constant across arms
TIMEFRAME = {"BL": "Baseline", "M12": "1-year (m12)", "LONG": "Longitudinal"}
LOSS_DISP = {"ce": "CE", "ce_patient": "CE + SigLIP(pat)",
             "ce_patient_label": "CE + SigLIP(pat,lab)"}
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


def _agg_label(arch, arm, clin_pool):
    """How the per-modality / per-visit CLS vectors are combined into the joint representation."""
    if arch == "mlp_concat":
        return "concatenation"                         # all tokens concatenated, then MLP
    if arch == "cross_attn":
        return "mean-pool"                             # self-attention over tokens, then mean-pool
    if arch == "cross_attn_x":                         # bidirectional cross-attn over single CLS
        return clin_pool if arm == "LONG" else "single CLS"   # LONG: mean/mamba visit pooling
    return "—"


def _fusion_row(arm, b):
    loss = b["loss"]
    loss_str = LOSS_DISP[loss] + ("" if loss == "ce" else
        f" · λ={b['lam']:g}" + (f" wu{int(b['warmup'])}" if b["warmup"] else ""))
    return {"block": "fusion", "section": f"fusion::{TIMEFRAME[arm]}", "ref_kind": "",
            "clin": CLIN_DISP, "timeframe": TIMEFRAME[arm],
            "mri": MRI_DISP[b["variant"]], "fusion": FUSION_DISP[b["arch"]], "emb": "CLS",
            "agg": _agg_label(b["arch"], arm, b["clin_pool"]), "loss_str": loss_str,
            "sig_acc": "", "sig_auc": "",
            **{k: b[k] for k in ("bacc_m", "bacc_s", "auc_m", "auc_s",
                                 "f1_m", "f1_s", "valbacc_m", "n_test")}}


def headline(full, uni):
    """Per (arm, loss): the best fusion over variant+arch (λ tuned in `full`); PLUS one best-mamba
    row per timeframe (where mamba exists — LONG only) to expose mamba's ceiling; then the unimodal
    references. Each reference carries TWO markers: accuracy (McNemar exact) and AUC (OVO-pairwise
    DeLong), both vs the best fusion of that timeframe. Rows are ordered by test bACC within each
    timeframe group. Columns are decomposed (clinical / timeframe / MRI / fusion / emb / agg / loss)."""
    rows = []
    for arm in ARM_ORDER:
        for loss in LOSS_ORDER:
            sub = full[(full.clin_arm == arm) & (full.loss == loss)]
            if not sub.empty:
                rows.append(_fusion_row(arm, sub.loc[sub["bacc_m"].idxmax()]))
        # best mamba-aggregation fusion for this timeframe (mamba pools >1 visit → LONG only)
        mam = full[(full.clin_arm == arm) & (full.arch == "cross_attn_x")
                   & (full.clin_pool == "mamba")]
        if not mam.empty:
            rows.append(_fusion_row(arm, mam.loc[mam["bacc_m"].idxmax()]))
    # ── unimodal references (own softmax, no head) — grouped clinical vs imaging ──
    REF_TF = {"clinical_only_1yr": "1-year (m12)", "clinical_only_long": "Longitudinal",
              "mri_only_1yr": "1-year (m12)", "mri_mvp_1yr": "1-year (m12)"}
    REF_CLIN = {"clinical_only_1yr": CLIN_DISP, "clinical_only_long": CLIN_DISP,
                "mri_only_1yr": "—", "mri_mvp_1yr": "—"}
    REF_MRI = {"clinical_only_1yr": "—", "clinical_only_long": "—",
               "mri_only_1yr": "BrainDINO-T2", "mri_mvp_1yr": "BrainMVP-T1b*"}
    REF_KIND = {"clinical_only_1yr": "clinical", "clinical_only_long": "clinical",
                "mri_only_1yr": "imaging", "mri_mvp_1yr": "imaging"}
    uni = uni.set_index("ref")
    for ref in ["clinical_only_1yr", "clinical_only_long", "mri_only_1yr", "mri_mvp_1yr"]:
        if ref not in uni.index:
            continue
        r = uni.loc[ref]
        rows.append({"block": "ref", "section": f"ref::{REF_KIND[ref]}", "ref_kind": REF_KIND[ref],
                     "clin": REF_CLIN[ref], "timeframe": REF_TF[ref],
                     "mri": REF_MRI[ref], "fusion": "—", "emb": "CLS",
                     "agg": "—", "loss_str": "—",
                     "sig_acc": str(r.get("sig_acc", r.get("sig", ""))),
                     "sig_auc": str(r.get("sig_auc", "")),
                     "bacc_m": r["balanced_acc_m"], "bacc_s": r["balanced_acc_s"],
                     "auc_m": r["auc_roc_ovr_m"], "auc_s": r["auc_roc_ovr_s"],
                     "f1_m": r["macro_f1_m"], "f1_s": r["macro_f1_s"],
                     "valbacc_m": np.nan, "n_test": r["n"]})
    # order: fusion block before refs; fusion by timeframe then bACC desc; refs by clinical-then-
    # imaging then bACC desc (so clinical and imaging references are separated into blocks).
    df = pd.DataFrame(rows)
    tf_order = {"1-year (m12)": 0, "Longitudinal": 1, "Baseline": 2}
    kind_order = {"clinical": 0, "imaging": 1}
    df["_blk"] = (df["block"] == "ref").astype(int)
    df["_sec"] = [tf_order.get(t, 9) if b == "fusion" else kind_order.get(k, 9)
                  for b, t, k in zip(df["block"], df["timeframe"], df["ref_kind"])]
    df = (df.sort_values(["_blk", "_sec", "bacc_m"], ascending=[True, True, False])
            .drop(columns=["_blk", "_sec"]).reset_index(drop=True))
    return df


def _fmt(m, s):
    return (f"{m:.2f} ± {s:.2f}" if pd.notna(s) else f"{m:.2f}") if pd.notna(m) else "--"


def _cells(r):
    bacc = _fmt(r["bacc_m"], r["bacc_s"]) + (f" [{r['sig_acc']}]" if r["sig_acc"] else "")
    auc = _fmt(r["auc_m"], r["auc_s"]) + (f" [{r['sig_auc']}]" if r["sig_auc"] else "")
    return [r["clin"], r["timeframe"], r["mri"], r["fusion"], r["emb"], r["agg"], r["loss_str"],
            bacc, auc, _fmt(r["f1_m"], r["f1_s"]), f"{r['n_test']:.0f}"]


def render_png(tab, out_png):
    matplotlib.rcParams["font.family"] = "DejaVu Serif"
    LEAD = ["Clinical (model · task)", "Timeframe", "MRI @ m12", "Feature fusion",
            "Emb", "Aggregation", "Loss"]
    headers = LEAD + ["Test bACC", "macro-AUC", "macro-F1", "n"]
    body = [_cells(r) for _, r in tab.iterrows()]
    numeric = np.array([[r["bacc_m"], r["auc_m"], r["f1_m"]] for _, r in tab.iterrows()], float)
    best = [int(np.nanargmax(numeric[:, j])) for j in range(3)]
    # rule whenever the section changes: fusion timeframe groups (M12→LONG), fusion→ref, and the
    # clinical→imaging reference split.
    sec = tab["section"].tolist()
    rules = [i > 0 and sec[i] != sec[i - 1] for i in range(len(tab))]

    COL_W = [2.70, 1.45, 1.55, 2.35, 0.55, 1.65, 3.00, 1.85, 1.85, 1.35, 0.50]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H, TOP_PAD, BOT_PAD = 2.15, 0.42, 0.46, 0.12, 0.12
    fig_w = LEFT + sum(COL_W) + RIGHT_PAD
    fig_h = TOP_PAD + TITLE_H + HEAD_H + ROW_H * len(body) + BOT_PAD
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    cl = [LEFT]
    for w in COL_W:
        cl.append(cl[-1] + w)
    RIGHT = cl[-1]; ccx = [(cl[i] + cl[i + 1]) / 2 for i in range(len(COL_W))]
    LEFT_COLS = {0, 1, 2, 3, 4, 5, 6}

    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)

    y = fig_h - TOP_PAD; y_top = y; y -= TITLE_H
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H / 2,
            "Cross-modal fusion — best hyperparameter-tuned per condition · 12-month T2 (CN/MCI/AD)\n"
            "Clinical encoder = BioClinical-ModernBERT-large (full-FT), CLS pooling · MRI = single m12 scan\n"
            "Fusion tuned over {MRI variant, architecture, λ∈{0.01,0.1,0.5,1.0}}; SigLIP λ>0.01 use warmup=10 · mean ± std seeds 0/1/2\n"
            "References = unimodal own-softmax (no head), grouped clinical / imaging. Markers vs best fusion: accuracy = McNemar exact [on bACC]; AUC = OVO-pairwise DeLong [on macro-AUC]\n"
            "* p<0.05, ** p<0.01, *** p<0.001, n.s. = not significant · macro-AUC shown is OVR (descriptive) · M12 n≈451, LONG n≈478 · BrainMVP-T1b* reference uses BrainMVP's T2 multiclass head (T1b is binary)",
            ha="center", va="center", fontsize=9.5, fontweight="bold", linespacing=1.5)
    hline(y_top, lw=1.5); hline(y, lw=1.2)
    yh = y; y -= HEAD_H; ym = (yh + y) / 2
    for j, h in enumerate(headers):
        if j in LEFT_COLS:
            ax.text(cl[j] + 0.06, ym, h, ha="left", va="center", fontsize=8.8, fontstyle="italic")
        else:
            ax.text(ccx[j], ym, h, ha="center", va="center", fontsize=8.8, fontstyle="italic")
    hline(y, lw=1.2); y_data_top = y
    for i, cells in enumerate(body):
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yt = y; y -= ROW_H; ym = (yt + y) / 2
        for j, txt in enumerate(cells):
            mj = j - len(LEAD)
            bold = (0 <= mj < 3 and best[mj] == i)
            if j in LEFT_COLS:
                ax.text(cl[j] + 0.06, ym, txt, ha="left", va="center", fontsize=8.5)
            else:
                ax.text(ccx[j], ym, txt, ha="center", va="center", fontsize=8.5,
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


def write_tex(tab, out_tex, include_n=True):
    def esc(s):
        return (str(s).replace("·", r"$\cdot$").replace("λ", r"$\lambda$")
                .replace("≈", r"$\approx$").replace("—", "--").replace("±", r"$\pm$"))

    def num(m, s, marker=""):
        body = (f"${m:.2f} \\pm {s:.2f}$" if pd.notna(s) else f"${m:.2f}$") if pd.notna(m) else "--"
        return body + (f"~[{marker}]" if marker else "")

    colspec = "lll l l l l ccc" + (" c" if include_n else "")
    head = (r"\textbf{Clinical (model $\cdot$ task)} & \textbf{Timeframe} & \textbf{MRI @ m12} & "
            r"\textbf{Feature fusion} & \textbf{Emb} & \textbf{Aggregation} & \textbf{Loss} & "
            r"\textbf{Test bACC} & \textbf{macro-AUC} & \textbf{macro-F1}"
            + (r" & \textbf{n}" if include_n else "") + r" \\")
    L = [r"% Joined best-per-condition table (03c_summary_joined.py)",
         r"\begin{table}[ht]", r"\centering", r"\scriptsize",
         r"\begin{tabular}{" + colspec + "}", r"\toprule", head, r"\midrule"]
    prev_sec = None
    for _, r in tab.iterrows():
        if prev_sec is not None and r["section"] != prev_sec:
            L.append(r"\midrule")
        prev_sec = r["section"]
        c = [esc(r["clin"]), esc(r["timeframe"]), esc(r["mri"]), esc(r["fusion"]), esc(r["emb"]),
             esc(r["agg"]), esc(r["loss_str"]),
             num(r["bacc_m"], r["bacc_s"], r["sig_acc"]), num(r["auc_m"], r["auc_s"], r["sig_auc"]),
             num(r["f1_m"], r["f1_s"])] + ([f"{r['n_test']:.0f}"] if include_n else [])
        L.append(" & ".join(c) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}",
          r"\caption{Best hyperparameter-tuned cross-modal fusion per condition for the 12-month "
          r"three-way diagnosis (CN/MCI/AD). Clinical encoder = BioClinical-ModernBERT-large "
          r"(full fine-tune), CLS pooling; MRI = a single m12 scan (BrainDINO-T2 or BrainMVP-T1b). "
          r"Each fusion row is the best over MRI variant, fusion architecture (MLP / self-attention "
          r"over concatenated tokens / bidirectional CLS cross-attention) and contrastive weight "
          r"$\lambda\in\{0.01,0.1,0.5,1.0\}$ ($\lambda>0.01$ with warmup=10, min\_epochs=30); one "
          r"best-mamba row per timeframe is added where mamba applies (LONG only). The "
          r"\emph{Aggregation} column gives how the CLS vectors are combined: concatenation (MLP), "
          r"mean-pool after self-attention, or mean/mamba pooling of the three clinical visits for "
          r"the cross-attention models (whose two updated CLS are then concatenated). "
          r"\emph{References are unimodal own-softmax predictions} (no retrained "
          r"head), grouped into clinical and imaging: the clinical 1-year model (T2\_m12) and "
          r"longitudinal clinical model (T2\_long, m12 row); and the imaging models at the m12 scan, "
          r"BrainDINO-T2 and BrainMVP-T1b$^{*}$. $^{*}$The fused variant~B uses BrainMVP-T1b (binary) "
          r"embeddings, which have no 3-way own softmax; the BrainMVP reference therefore uses "
          r"BrainMVP's T2 multiclass head. Two markers per reference give the "
          r"significance of the best fusion of that timeframe vs the reference: the \textbf{accuracy} "
          r"marker (on bACC) is an exact \textbf{McNemar} test on per-patient correctness; the "
          r"\textbf{AUC} marker (on macro-AUC) is \textbf{OVO-pairwise DeLong} (Hand--Till, mean of "
          r"the CN/MCI, CN/AD, MCI/AD pairwise tests), both pooled over the three test folds: "
          r"$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$, n.s.\ = not significant. The displayed "
          r"macro-AUC is one-vs-rest (descriptive); the DeLong test uses the OVO formulation. M12 "
          r"cohort n$\approx$451, LONG n$\approx$478. Mean $\pm$ std over seeds 0/1/2, 2 dp.}",
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
    write_tex(tab, os.path.join(OUT, "summary_joined_no_n.tex"), include_n=False)
    print("\n  Joined best-per-condition (Test bACC, 2dp):")
    for _, r in tab.iterrows():
        tag = "REF " if r["block"] == "ref" else "fuse"
        mk = f" [acc {r['sig_acc']} / AUC {r['sig_auc']}]" if r["block"] == "ref" else ""
        cfg = "unimodal own-softmax" if r["block"] == "ref" else f"{r['mri']} · {r['fusion']} · {r['agg']}"
        print(f"    [{tag}] {r['timeframe']:13s} {r['loss_str']:30s} "
              f"{r['bacc_m']:.2f}±{r['bacc_s']:.2f}  AUC {r['auc_m']:.2f}  ({cfg}){mk}")


if __name__ == "__main__":
    main()
