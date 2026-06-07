"""
04_summary_t3_cleaned.py
========================
Cleaned house-style copy of the T3 per-horizon late-fusion summary (clinical@bl + MRI@m12), matching
the T2/T1d cleaned tables (DejaVu-Serif, bordered, ruled, italic headers, dashed dividers, bold
best-per-block). Long layout grouped by horizon x cohort; the exact clinical / MRI model is named in
the Method column (clinical-only / MRI-only rows). Metrics = Test bACC (seed-mean +/- std), macro-AUC,
F1 conv. Metrics read from SUMMARY_T3_present_only_per_seed.csv; the per-seed non-conv:conv breakdown
is computed by reusing the 03 loaders (no fusion re-run).

Out: integration/T3/summary_t3_cleaned.{csv,png,pdf}
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
                recs.append({"hk": hk, "ck": ck, "Horizon": hlab, "Cohort": clab,
                             "Method": method_label(v, hk), "bacc": vd.bacc.mean(),
                             "bstd": vd.bacc.std(), "auc": vd.macro_auc.mean(),
                             "f1c": vd.f1_conv.mean(), "n": vd.n.mean()})
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
    df[["Horizon", "Cohort", "Method", "bacc", "bstd", "auc", "f1c", "n"]].rename(
        columns={"bacc": "Test_bACC_mean", "bstd": "Test_bACC_std", "auc": "macro_AUC",
                 "f1c": "F1_conv"}).to_csv(os.path.join(HERE, "summary_t3_cleaned.csv"), index=False)
    render(df, subtitle2, os.path.join(HERE, "summary_t3_cleaned.png"))
    print(f"  CSV: {os.path.join(HERE, 'summary_t3_cleaned.csv')}  ({len(df)} rows)")
    print("  breakdown:", subtitle2)


if __name__ == "__main__":
    main()
