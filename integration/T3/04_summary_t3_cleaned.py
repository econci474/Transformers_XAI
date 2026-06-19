"""
04_summary_t3_cleaned.py
========================
Cleaned house-style copy of the T3 per-horizon late-fusion summary (clinical@bl + MRI@m12), matching
the T2/T1d cleaned tables (DejaVu-Serif, bordered, ruled, italic headers, dashed dividers, bold
best-per-block). Long layout grouped by horizon x cohort; the exact clinical / MRI model is named in
the Method column (clinical-only / MRI-only rows). Metrics = Test bACC (seed-mean +/- std), macro-AUC,
F1 conv. Metrics read from SUMMARY_T3_present_only_per_seed.csv; the per-seed non-conv:conv breakdown
is computed by reusing the 03 loaders (no fusion re-run).

Out: integration/T3/summary_t3_cleaned.{csv,png,pdf,tex}
     integration/T3/summary_t3_filtered.{csv,png,pdf,tex}  (top fusion + unimodal refs per group)
     integration/T3/summary_t3_small.{csv,png,pdf,tex}     (best configuration per window)
Run: python integration/T3/04_summary_t3_cleaned.py
"""
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = os.path.dirname(os.path.abspath(__file__))

HZ = [("3yr", "T3a (≤3y)"), ("5yr", "T3b (≤5y)"), ("7yr", "T3d (≤7y)")]
COH = [("both_present", "both-present"), ("present_only", "all-available")]
METHODS = ["clinical_only", "mri_only", "equal_0.5", "weighted_avg", "weighted_avg_J", "elastic_net"]
# exact model abbreviations (named in the table, per horizon)
CLIN_ABBR = {"3yr": "BioClin-L-ft", "5yr": "BioClin-B-ft", "7yr": "MBERT-L-ft"}
MRI_ABBR = {"3yr": "BrainMVP-ft (plus-orig)", "5yr": "BrainMVP-frozen (none)",
            "7yr": "BrainMVP-ft (plus-orig)"}
FUSE_LABEL = {"equal_0.5": "weighted-avg (equal)", "weighted_avg": "weighted-avg (val bACC)",
              "weighted_avg_J": "weighted-avg (Youden J)", "elastic_net": "elastic-net"}


def method_label(v, hk):
    if v == "clinical_only":
        return CLIN_ABBR[hk]
    if v == "mri_only":
        return MRI_ABBR[hk]
    return FUSE_LABEL[v]


def cohort_breakdown():
    """Per-seed mean (non-conv, conv) counts per (horizon, cohort), reusing 03's loaders."""
    spec = importlib.util.spec_from_file_location("t3fuse",
                                                  os.path.join(HERE, "03_fuse_T3_present_only.py"))
    t3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t3)
    h = t3.h
    bd = {}
    for hk, cfg in t3.HORIZONS.items():
        tmpl = t3._mri_template(cfg)
        if not Path(tmpl.format(seed=0)).exists():
            continue
        acc = {"present_only": [], "both_present": []}
        for seed in t3.SEEDS:
            clin = h.load_clinical(seed, t3.CLIN_ROOT, cfg["model"], "full_ft", cfg["clin_task"])
            mri = h.load_mri(seed, tmpl, "m12")
            f = clin.merge(mri, on="Patient_ID", how="left")
            pres = ~np.isnan(f[h.MP_COLS].to_numpy(float)).any(axis=1)
            y = f["y_clin"].to_numpy(int); tm = (f["split"] == "test").to_numpy()
            yt, prt = y[tm], pres[tm]
            acc["present_only"].append([(yt == 0).sum(), (yt == 1).sum()])
            acc["both_present"].append([(yt[prt] == 0).sum(), (yt[prt] == 1).sum()])
        for ck in acc:
            a = np.array(acc[ck], float).mean(0)
            bd[(hk, ck)] = (a[0], a[1])
    return bd


def aggregate():
    d = pd.read_csv(os.path.join(HERE, "SUMMARY_T3_present_only_per_seed.csv"))
    recs = []
    for hk, hlab in HZ:
        for ck, clab in COH:
            for v in METHODS:
                vd = d[(d.horizon == hk) & (d.cohort == ck) & (d.variant == v)]
                if vd.empty:
                    continue
                recs.append({"hk": hk, "ck": ck, "variant": v, "Horizon": hlab, "Cohort": clab,
                             "Method": method_label(v, hk), "bacc": vd.bacc.mean(),
                             "bstd": vd.bacc.std(), "auc": vd.macro_auc.mean(),
                             "f1c": vd.f1_conv.mean(), "mf1": vd.macro_f1.mean(), "n": vd.n.mean()})
    return pd.DataFrame(recs)


def render(df, subtitle2, out_path):
    headers = ["Horizon", "Cohort", "Method", "Test bACC", "macro-AUC", "F1 conv", "n"]
    LEFT_COLS = {0, 1, 2}
    best = {}
    for (hk, ck), g in df.groupby(["hk", "ck"], sort=False):
        for key in ("bacc", "auc", "f1c"):
            best[(hk, ck, key)] = g[key].idxmax()

    body, rules, bolds, show0 = [], [], [], []
    prev_hk, prev_ck = None, None
    for idx, r in df.iterrows():
        sh = r.hk != prev_hk
        sc = (r.hk, r.ck) != (prev_hk, prev_ck)
        body.append([r.Horizon if sh else "", r.Cohort if sc else "", r.Method,
                     f"{r.bacc:.3f} ± {r.bstd:.3f}", f"{r.auc:.3f}", f"{r.f1c:.3f}", f"{r.n:.0f}"])
        bolds.append({3: best[(r.hk, r.ck, "bacc")] == idx, 4: best[(r.hk, r.ck, "auc")] == idx,
                      5: best[(r.hk, r.ck, "f1c")] == idx})
        rules.append("horizon" if (prev_hk is not None and r.hk != prev_hk)
                     else ("cohort" if (prev_ck is not None and r.ck != prev_ck) else None))
        show0.append(sh)
        prev_hk, prev_ck = r.hk, r.ck

    COL_W = [1.35, 1.55, 2.95, 1.62, 1.15, 1.05, 0.55]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H = 1.05, 0.40, 0.42
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
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H * 0.68, "T3 Late Fusion",
            ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H * 0.26,
            "mean ± std across seeds 0/1/2  ·  both-present (complete-case) vs "
            "all-available (full clinical cohort)\n" + subtitle2,
            ha="center", va="center", fontsize=8.3, fontstyle="italic", linespacing=1.45)
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
        if rules[i] == "horizon":
            hline(y, lw=1.0)
        elif rules[i] == "cohort":
            hline(y, lw=0.5, ls=(0, (3, 3)))
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                w = "bold" if (j == 0 and show0[i]) else "normal"
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center",
                        fontsize=9.0, fontweight=w)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bolds[i].get(j, False) else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def write_t3_latex(df, out_path, subtitle2, label="tab:t3_late_fusion", caption_note=""):
    """Portrait booktabs LaTeX of the T3 late-fusion table, organised as Window | Clinical | MRI |
    Aggregation | Missingness + Test bACC / macro-AUC / macro-F1 / n. Window `\\midrule`, cohort
    `\\cmidrule`; bold-best per (window, missingness) group; bottom caption with N-breakdown + glossary."""
    def esc(s):
        return (str(s).replace("≤", r"$\leq$").replace("±", r"$\pm$").replace("—", "--")
                .replace("·", r"$\cdot$").replace("&", r"\&").replace("%", r"\%")
                .replace("_", r"\_").replace("#", r"\#"))
    best = {}
    for (hk, ck), g in df.groupby(["hk", "ck"], sort=False):
        for key in ("bacc", "auc", "mf1"):
            best[(hk, ck, key)] = g[key].idxmax()

    def num(x, bold):
        return f"$\\mathbf{{{x}}}$" if bold else f"${x}$"

    L = [r"% T3 Late Fusion -- per-window clinical@bl + MRI@m12 late fusion.",
         r"% Requires: \usepackage{booktabs}",
         r"\begin{table}[ht]", r"\centering", r"\footnotesize",
         r"\begin{tabular}{lllll ccc c}", r"\toprule",
         r"\multicolumn{9}{c}{\textbf{T3 Late Fusion}} \\", r"\midrule",
         r"\textit{Window} & \textit{Clinical} & \textit{MRI} & \textit{Aggregation} & "
         r"\textit{Missingness} & \textit{Test bACC} & \textit{macro-AUC} & \textit{macro-F1} & "
         r"\textit{n} \\", r"\midrule"]
    prev_hk, prev_ck = None, None
    for idx, r in df.iterrows():
        new_hz = r.hk != prev_hk
        new_co = (r.hk, r.ck) != (prev_hk, prev_ck)
        if prev_hk is not None and new_hz:
            L.append(r"\midrule")
        elif prev_ck is not None and new_co:
            L.append(r"\cmidrule(lr){2-9}")
        clin, mri, agg, miss = _decompose_t3(r)
        win = esc(r.Horizon) if new_hz else ""
        miss_d = esc(miss) if new_co else ""
        bacc = num(f"{r.bacc:.3f} \\pm {r.bstd:.3f}", best[(r.hk, r.ck, 'bacc')] == idx)
        auc = num(f"{r.auc:.3f}", best[(r.hk, r.ck, 'auc')] == idx)
        mf1 = num(f"{r.mf1:.3f}", best[(r.hk, r.ck, 'mf1')] == idx)
        L.append(f"{win} & {esc(clin)} & {esc(mri)} & {esc(agg)} & {miss_d} & {bacc} & {auc} & {mf1} "
                 f"& {r.n:.0f} " + r"\\")
        prev_hk, prev_ck = r.hk, r.ck
    cap = (r"T3 per-window late fusion (clinical@bl $\oplus$ MRI@m12), mean $\pm$ std across seeds 0/1/2. "
           r"\textit{Window} = conversion horizon; \textit{Clinical} / \textit{MRI} = the encoders fused "
           r"(-- = single-modality reference); \textit{Aggregation} = the late-fusion combiner; "
           r"\textit{Missingness} = complete-case (both modalities present) vs all-available (clinical "
           r"proceeds where MRI absent). \textit{macro-F1} = mean of the two per-class F1. Per-window "
           r"clinical / MRI models: T3a ($\leq$3y) BioClin-L-ft / BrainMVP-ft (plus-orig); T3b ($\leq$5y) "
           r"BioClin-B-ft / BrainMVP-frozen (none); T3d ($\leq$7y) MBERT-L-ft / BrainMVP-ft (plus-orig). "
           r"weighted-avg (equal) = fixed equal weights; (val bACC) / (Youden J) = tuned on validation; "
           r"elastic-net = stacked logistic. " + esc(subtitle2) + ". Bold = best per window $\times$ "
           r"missingness." + caption_note)
    L += [r"\bottomrule", r"\end{tabular}", r"\caption{" + cap + r"}", r"\label{" + label + r"}",
          r"\end{table}"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  TEX: {out_path}")


MISS_LABEL = {"both_present": "complete-case", "present_only": "all-available"}


def _decompose_t3(r):
    """variant + horizon -> (Clinical, MRI, Aggregator, Missingness)."""
    clin = "—" if r.variant == "mri_only" else CLIN_ABBR[r.hk]
    mri = "—" if r.variant == "clinical_only" else MRI_ABBR[r.hk]
    agg = "—" if r.variant in ("clinical_only", "mri_only") else FUSE_LABEL[r.variant]
    return clin, mri, agg, MISS_LABEL[r.ck]


def window_refs_and_best(df):
    """Per window: the clinical-only & MRI-only REFERENCE models + the best INTEGRATION (fusion) method,
    taken from the cohort that holds the best fusion ('best cohort' per window). The MRI-only reference
    is absent for an all-available cohort (MRI not present for those subjects). Splits each row into
    Clinical | MRI | Aggregator | Missingness."""
    FUS = ["equal_0.5", "weighted_avg", "weighted_avg_J", "elastic_net"]
    keep = []
    for hk, _ in HZ:
        g = df[df.hk == hk]
        fus_all = g[g.variant.isin(FUS)]
        ck = fus_all.loc[fus_all.bacc.idxmax(), "ck"]        # cohort holding the best fusion
        gc = g[g.ck == ck]
        bf = gc[gc.variant.isin(FUS)]
        keep.append(pd.concat([gc[gc.variant == "clinical_only"],
                               gc[gc.variant == "mri_only"],
                               bf.loc[[bf.bacc.idxmax()]]]))
    b = pd.concat(keep).reset_index(drop=True)
    dec = [_decompose_t3(r) for _, r in b.iterrows()]
    b["Clinical"] = [d[0] for d in dec]; b["MRI"] = [d[1] for d in dec]
    b["Aggregator"] = [d[2] for d in dec]; b["Missingness"] = [d[3] for d in dec]
    return b


def render_best(df, subtitle2, out_path, show_n=True):
    headers = ["Window", "Clinical", "MRI", "Aggregation", "Missingness",
               "Test bACC", "macro-AUC", "macro-F1"] + (["n"] if show_n else [])
    LEFT_COLS = {0, 1, 2, 3, 4}
    METRIC_J0 = 5
    best = {}
    for hk, g in df.groupby("hk", sort=False):
        for j, key in enumerate(("bacc", "auc", "mf1")):
            best[(hk, j)] = g[key].idxmax()

    body, rules, bolds, show_w = [], [], [], []
    prev_hk = None
    for idx, r in df.iterrows():
        sw = r.hk != prev_hk
        cells = [r.Horizon if sw else "", r.Clinical, r.MRI, r.Aggregator,
                 r.Missingness if sw else "",
                 f"{r.bacc:.3f} ± {r.bstd:.3f}", f"{r.auc:.3f}", f"{r.mf1:.3f}"]
        if show_n:
            cells.append(f"{r.n:.0f}")
        body.append(cells)
        bolds.append({METRIC_J0 + j: best[(r.hk, j)] == idx for j in range(3)})
        rules.append(prev_hk is not None and sw)
        show_w.append(sw)
        prev_hk = r.hk

    COL_W = [1.45, 1.65, 2.45, 2.15, 1.75, 1.62, 1.15, 1.15] + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB = ("clinical-only & MRI-only references + the best integration (fusion) method per window "
           "(best-scoring missingness cohort)\n" + subtitle2)
    SUB_H = 0.30 * (SUB.count("\n") + 1)
    TITLE_H, HEAD_H, ROW_H = 0.62 + SUB_H, 0.40, 0.44
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
    ax.text((LEFT + RIGHT) / 2, y + TITLE_H - 0.26,
            "T3 Late Fusion — references vs best integration per window",
            ha="center", va="center", fontsize=13, fontweight="bold")
    ax.text((LEFT + RIGHT) / 2, y + SUB_H / 2 + 0.02, SUB, ha="center", va="center",
            fontsize=8.3, fontstyle="italic", linespacing=1.45)
    hline(y_title_top, lw=1.5); hline(y, lw=1.2)

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
            hline(y, lw=0.7)
        yt = y; y -= ROW_H
        ym = (yt + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                w = "bold" if (j == 0 and show_w[i]) else "normal"
                ax.text(col_left[j] + 0.06, ym, cells[j], ha="left", va="center",
                        fontsize=9.0, fontweight=w)
            else:
                ax.text(col_cx[j], ym, cells[j], ha="center", va="center", fontsize=9.0,
                        fontweight="bold" if bolds[i].get(j, False) else "normal")
    BOTTOM = y
    for x in col_left[1:-1]:
        ax.plot([x, x], [BOTTOM, y_data_top], color="black", linewidth=0.7,
                linestyle=(0, (3, 3)), zorder=2)
    hline(BOTTOM, lw=1.5)
    ax.add_patch(plt.Rectangle((LEFT, BOTTOM), RIGHT - LEFT, y_title_top - BOTTOM,
                               facecolor="none", edgecolor="black", linewidth=1.5, zorder=5))
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def write_t3_best_latex(df, out_path, subtitle2, label="tab:t3_best_per_window"):
    """Portrait booktabs LaTeX: per window, the clinical-only & MRI-only reference models + the best
    integration method. Window | Clinical | MRI | Aggregation | Missingness + Test bACC / macro-AUC /
    macro-F1 / n. `\\midrule` between windows; bold-best per window."""
    def esc(s):
        return (str(s).replace("≤", r"$\leq$").replace("±", r"$\pm$").replace("·", r"$\cdot$")
                .replace("—", "--").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
                .replace("#", r"\#"))
    best = {}
    for hk, g in df.groupby("hk", sort=False):
        for key in ("bacc", "auc", "mf1"):
            best[(hk, key)] = g[key].idxmax()

    def num(x, bold):
        return f"$\\mathbf{{{x}}}$" if bold else f"${x}$"

    L = [r"% T3 Late Fusion -- references vs best integration per conversion window.",
         r"% Requires: \usepackage{booktabs}",
         r"\begin{table}[ht]", r"\centering", r"\small",
         r"\begin{tabular}{lllll ccc c}", r"\toprule",
         r"\multicolumn{9}{c}{\textbf{T3 Late Fusion --- references vs best integration per window}} \\",
         r"\midrule",
         r"\textit{Window} & \textit{Clinical} & \textit{MRI} & \textit{Aggregation} & "
         r"\textit{Missingness} & \textit{Test bACC} & \textit{macro-AUC} & \textit{macro-F1} & "
         r"\textit{n} \\", r"\midrule"]
    prev_hk = None
    for idx, r in df.iterrows():
        sw = r.hk != prev_hk
        if prev_hk is not None and sw:
            L.append(r"\midrule")
        win = esc(r.Horizon) if sw else ""
        miss = esc(r.Missingness) if sw else ""
        bacc = num(f"{r.bacc:.3f} \\pm {r.bstd:.3f}", best[(r.hk, 'bacc')] == idx)
        auc = num(f"{r.auc:.3f}", best[(r.hk, 'auc')] == idx)
        mf1 = num(f"{r.mf1:.3f}", best[(r.hk, 'mf1')] == idx)
        L.append(f"{win} & {esc(r.Clinical)} & {esc(r.MRI)} & {esc(r.Aggregator)} & {miss} & "
                 f"{bacc} & {auc} & {mf1} & {r.n:.0f} " + r"\\")
        prev_hk = r.hk
    cap = (r"T3 per-window late fusion: the clinical-only and MRI-only \textbf{reference models} and the "
           r"\textbf{best integration} (fusion) method per conversion window, drawn from the cohort that "
           r"holds the best fusion. \textit{Aggregation} = the late-fusion combiner; \textit{Missingness} "
           r"= complete-case (both modalities present) vs all-available (clinical proceeds where MRI "
           r"absent; no MRI-only reference exists there). \textit{macro-F1} = mean of the two per-class "
           r"F1. At $\leq$7y the clinical-only reference still beats fusion. mean $\pm$ std across seeds "
           r"0/1/2. " + esc(subtitle2) + ". Bold = best per window.")
    L += [r"\bottomrule", r"\end{tabular}", r"\caption{" + cap + r"}", r"\label{" + label + r"}",
          r"\end{table}"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  TEX: {out_path}")


def main():
    df = aggregate()
    bd = cohort_breakdown()
    parts = []
    for hk, hlab in HZ:
        if (hk, "present_only") in bd:
            nc, cv = bd[(hk, "present_only")]
            parts.append(f"{hlab} {nc:.0f}:{cv:.0f}")
    subtitle2 = ("n = N(test) per seed.  non-conv : conv per horizon (all-available) — "
                 + "  ·  ".join(parts))

    out_cols = ["Horizon", "Cohort", "Method", "bacc", "bstd", "auc", "mf1", "f1c", "n"]
    ren = {"bacc": "Test_bACC_mean", "bstd": "Test_bACC_std", "auc": "macro_AUC",
           "mf1": "macro_F1", "f1c": "F1_conv"}
    df[out_cols].rename(columns=ren).to_csv(os.path.join(HERE, "summary_t3_cleaned.csv"), index=False)
    render(df, subtitle2, os.path.join(HERE, "summary_t3_cleaned.png"))
    write_t3_latex(df, os.path.join(HERE, "summary_t3_cleaned.tex"), subtitle2, "tab:t3_late_fusion")
    print(f"  CSV: {os.path.join(HERE, 'summary_t3_cleaned.csv')}  ({len(df)} rows)")

    # filtered: per (horizon, cohort) keep the unimodal reference models + the top fusion (max bACC)
    REFS = {"clinical_only", "mri_only"}
    keep = []
    for (hk, ck), g in df.groupby(["hk", "ck"], sort=False):
        keep.append(g[g.variant.isin(REFS)])
        fus = g[~g.variant.isin(REFS)]
        if not fus.empty:
            keep.append(fus.loc[[fus.bacc.idxmax()]])
    filt = pd.concat(keep).sort_index().reset_index(drop=True)
    sub_filt = "top fusion combination + unimodal reference models  ·  " + subtitle2
    filt[out_cols].rename(columns=ren).to_csv(os.path.join(HERE, "summary_t3_filtered.csv"), index=False)
    render(filt, sub_filt, os.path.join(HERE, "summary_t3_filtered.png"))
    write_t3_latex(filt, os.path.join(HERE, "summary_t3_filtered.tex"), subtitle2,
                   "tab:t3_late_fusion_filtered",
                   caption_note=(r" Filtered to the top fusion combination (weighted-avg equal) and the "
                                 r"unimodal reference models; full results in "
                                 r"Table~\ref{tab:t3_late_fusion}."))
    print(f"  filtered ({len(filt)} rows)")

    # summary_t3_small: per window, the unimodal reference models + the best integration method
    # (Window | Clinical | MRI | Aggregation | Missingness + Test bACC / macro-AUC / macro-F1 / n)
    best = window_refs_and_best(df)
    small_cols = ["Horizon", "Clinical", "MRI", "Aggregator", "Missingness", "bacc", "bstd", "auc",
                  "mf1", "n"]
    small_ren = {**ren, "Horizon": "Window", "Aggregator": "Aggregation"}
    best[small_cols].rename(columns=small_ren).to_csv(
        os.path.join(HERE, "summary_t3_small.csv"), index=False)
    render_best(best, subtitle2, os.path.join(HERE, "summary_t3_small.png"), show_n=True)
    render_best(best, subtitle2, os.path.join(HERE, "summary_t3_small_no_n.png"), show_n=False)
    write_t3_best_latex(best, os.path.join(HERE, "summary_t3_small.tex"), subtitle2)
    print(f"  small ({len(best)} rows)")
    print("  breakdown:", subtitle2)


if __name__ == "__main__":
    main()
