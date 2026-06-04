"""
04_summary_top_methods.py
=========================
Cross-timeframe T2 summary: gather every `fusion_metrics.csv` under outputs/{baseline_only, m12_only,
up_to_m12}, unify the (differing) schemas, and render the TOP 5 methods by balanced accuracy.

Rules (locked with user):
  - include reference rows AND fusion methods together;
  - DEDUP rows with identical metrics (collapses clinical_only [CL_m12] reported in many subfolders);
  - EXCLUDE rows scored on n < 30 test patients (drops artefacts like baseline_only weighted-avg =
    bACC 1.000 on n=2);
  - rank by bACC mean (over seeds 0/1/2), take top 5.

CAVEAT (stated on the figure): the cohort/n differs across timeframes (baseline base = CL_bl @bl;
m12 / up_to_m12 base = CL_m12 @m12; complete-case subsets are smaller), so bACC is NOT on an identical
cohort across rows -- read every number together with n.

Outputs (at outputs/):
  summary_top5_bACC.png            top-5 table (best-per-column in bold) + criteria/caveat footnote
  summary_methods_all_ranked.csv   every method (n>=30 ranked first, then n<30 flagged), with provenance

Run:  python integration/T2_classification/04_summary_top_methods.py
"""
from __future__ import annotations

import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "outputs"
METRIC_MEANS = ["bacc_mean", "macro_f1_mean", "macro_auc_mean",
                "f1_CN_mean", "f1_MCI_mean", "f1_AD_mean"]
MIN_N = 30
LABEL_MAP = {  # tidy the raw id-column tokens
    "clinical_only": "clinical-only", "mri_only": "MRI-only", "MRI_only": "MRI-only",
    "temporal_CL_T2_EN": "clinical-temporal-EN", "clinical_temporal_EN": "clinical-temporal-EN",
    "elastic_net": "elastic-net", "weighted_avg": "weighted-avg", "weighted_avg_J": "weighted-avg-J",
    "detector_structured": "detector-structured", "detector_EN": "detector-EN",
    "detector_augmented": "detector-augmented", "lr": "logreg",
}


def pretty_method(row):
    """Clean display name: no timeframe prefix, no [timeframe/family] bracket; spell out the
    clinical source (CL T2 / CL T1/T1b) and put the WEIGHTING in a trailing [square bracket]."""
    m = str(row["method"]); tf = row["timeframe"]
    core = re.sub(r"^(M12X|M12_DET|M12|BL)\s*/\s*", "", m)          # drop the framing prefix
    if "MRI-only" in core:
        return core                                                # already names the BrainMVP model
    if "clinical-only" in core:
        base = {"m12_only": "@m12", "baseline_only": "@bl"}.get(tf, "")
        return f"CL T2 {base}".strip()
    if "clinical-temporal-EN" in core:
        return "CL T2 @bl+m06+m12  [learned weights: elastic-net]"
    mm = re.match(r"^([AB])_(LR|EN|meanP)\s*/\s*(CL\+MRI|CL)$", core)
    if mm:
        form, comb, src = mm.groups()
        clin = "CL T2 @bl+m06+m12" if form == "A" else "CL T1/T1b @bl+m06+m12"
        mri = " + MRI(BrainMVP T1/T1b stoch @m12)" if src == "CL+MRI" else ""
        wt = {"meanP": "equal weighting", "LR": "learned weights: logreg",
              "EN": "learned weights: elastic-net"}[comb]
        return f"{clin}{mri}  [{wt}]"
    return core                                                    # detectors / weighted-avg / etc.


def _label(row, idcols):
    parts = []
    for c in idcols:
        v = str(row[c])
        if v in ("-", "nan", "", "None") or pd.isna(row[c]):
            continue
        parts.append(LABEL_MAP.get(v, v))
    return " / ".join(parts) if parts else "(unnamed)"


def collect():
    files = glob.glob(str(OUT / "**" / "fusion_metrics.csv"), recursive=True)
    # shallower paths first so references are attributed to the top-level timeframe folder
    files.sort(key=lambda p: (len(Path(p).relative_to(OUT).parts), p))
    rows = []
    for f in files:
        rel = Path(f).relative_to(OUT).parts
        timeframe = rel[0]
        family = "/".join(rel[1:-1])
        where = timeframe + (f" / {family}" if family else "")
        d = pd.read_csv(f)
        if "bacc_mean" not in d.columns or "n_test_mean" not in d.columns:
            continue
        idcols = [c for c in d.columns if not any(
            c.startswith(p) for p in ("bacc", "macro", "f1_", "n_test"))]
        for _, r in d.iterrows():
            rec = {"timeframe": timeframe, "family": family, "where": where,
                   "method": _label(r, idcols), "source_file": str(Path(f).relative_to(OUT))}
            for m in METRIC_MEANS + ["bacc_std", "n_test_mean"]:
                rec[m] = float(r[m]) if m in d.columns and pd.notna(r[m]) else np.nan
            rows.append(rec)
    return pd.DataFrame(rows)


def main() -> int:
    df = collect()
    if df.empty:
        print("[error] no fusion_metrics.csv found under", OUT); return 1
    # always name the MRI model: the m12_only/flat MRI-only refs use the BrainMVP T2 3-class head
    mri_unnamed = df["method"].str.contains("MRI-only", case=False) & ~df["method"].str.contains("BrainMVP")
    df.loc[mri_unnamed, "method"] = df.loc[mri_unnamed, "method"] + " [BrainMVP T2 3-class, stoch @m12]"
    df["display"] = df.apply(pretty_method, axis=1)

    # dedup identical-metric rows (same underlying result reported in multiple folders)
    keytuple = lambda r: tuple(round(r[m], 4) if pd.notna(r[m]) else np.nan
                               for m in METRIC_MEANS + ["n_test_mean"])
    df["_key"] = df.apply(keytuple, axis=1)
    df["_dups"] = df.groupby("_key")["_key"].transform("size")
    # keep the SHALLOWEST source (attribute a duplicated reference to its top-level timeframe folder)
    df["_depth"] = df["source_file"].map(lambda s: len(Path(s).parts))
    df = df.sort_values(["_depth", "source_file"]).drop_duplicates("_key", keep="first")

    df["n_ge_30"] = df["n_test_mean"] >= MIN_N
    df = df.sort_values("bacc_mean", ascending=False).reset_index(drop=True)

    # full ranked CSV (n>=30 first, then the excluded low-n rows)
    inc = df[df.n_ge_30].copy(); inc.insert(0, "rank", range(1, len(inc) + 1))
    exc = df[~df.n_ge_30].copy(); exc.insert(0, "rank", "excluded_n<30")
    cols = (["rank", "timeframe", "family", "method"] + METRIC_MEANS
            + ["bacc_std", "n_test_mean", "_dups", "source_file"])
    full = pd.concat([inc, exc], ignore_index=True)[cols].rename(columns={"_dups": "n_folders_reporting"})
    full.to_csv(OUT / "summary_methods_all_ranked.csv", index=False)

    top = inc.head(5).copy()
    # pin the MRI-only @m12 reference (BrainMVP T2 3-class) even though it ranks low
    mri = inc[inc["method"].str.contains("MRI-only", case=False) & (inc["timeframe"] == "m12_only")]
    if len(mri) and mri.iloc[0]["method"] not in set(top["method"]):
        top = pd.concat([top, mri.head(1)], ignore_index=True)
    render_png(top, OUT / "summary_top5_bACC.png", "TOP 5 by Test balanced accuracy (+ MRI-only ref)")
    render_png(inc, OUT / "summary_methods_all_ranked.png",
               f"ALL methods (n≥{MIN_N}), ranked by Test balanced accuracy  [{len(inc)} methods]")

    print("=" * 80)
    print("  T2 cross-timeframe summary — TOP 5 by bACC (n>=30, refs+fusion, dedup)")
    print("=" * 80)
    show = ["rank", "method", "where", "bacc_mean", "macro_f1_mean", "macro_auc_mean", "n_test_mean"]
    print(top.assign(where=top["where"])[show].to_string(index=False))
    print(f"\n  wrote -> {OUT/'summary_top5_bACC.png'}\n         {OUT/'summary_methods_all_ranked.csv'}")
    return 0


# --------------------------------------------------------------------------- #
def render_png(top, out_path, subtitle):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [WARN] matplotlib missing -- skip PNG."); return
    COLS = [("bacc_mean", "Test bACC", "bacc_std"), ("macro_f1_mean", "macro-F1", None),
            ("macro_auc_mean", "macro-AUC", None), ("f1_CN_mean", "F1 CN", None),
            ("f1_MCI_mean", "F1 MCI", None), ("f1_AD_mean", "F1 AD", None)]
    headers = ["#", "Method"] + [c[1] for c in COLS] + ["n"]
    body, numeric = [], []
    for i, (_, r) in enumerate(top.iterrows(), 1):
        rk = str(int(r["rank"])) if "rank" in r.index and pd.notna(r["rank"]) else str(i)
        cells = [rk, r["display"]]
        nums = []
        for key, _, stdk in COLS:
            v = r[key]
            s = "-" if pd.isna(v) else (f"{v:.3f}" + (f" ± {r[stdk]:.3f}"
                                        if stdk and pd.notna(r.get(stdk)) else ""))
            cells.append(s); nums.append(v if pd.notna(v) else np.nan)
        cells.append(str(int(round(r["n_test_mean"]))))
        body.append(cells); numeric.append(nums)
    numeric = np.array(numeric, float)

    n_rows, n_cols = len(body), len(headers)
    col_chars = [max(len(str(x)) for x in [headers[j]] + [b[j] for b in body]) + 2
                 for j in range(n_cols)]
    total = sum(col_chars)
    fig, ax = plt.subplots(figsize=(total * 0.102, 0.7 + 0.4 * (n_rows + 2)))
    ax.axis("off")
    tab = ax.table(cellText=body, colLabels=headers, loc="upper center", cellLoc="center",
                   colLoc="center", colWidths=[c / total for c in col_chars])
    tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.0, 1.5)
    # bold best per metric column (cols 2..2+len(COLS)-1)
    for j, (key, _, _) in enumerate(COLS):
        col = numeric[:, j]
        if not np.all(np.isnan(col)):
            tab[int(np.nanargmax(col)) + 1, 2 + j].set_text_props(weight="bold")
    for j in range(n_cols):
        tab[0, j].set_facecolor("#ECECEC"); tab[0, j].set_text_props(weight="bold")
    for i in range(n_rows + 1):
        tab[i, 1].set_text_props(ha="left")
    plt.title(f"T2 late fusion — best methods across timeframes ({subtitle})\n"
              "baseline_only · m12_only · up_to_m12   (mean ± std across seeds 0/1/2, n = mean #TEST)",
              pad=8, fontsize=11)
    foot = [
        "ALL metrics are on the TEST fold (meta-learners fit per-seed on VAL); 'Test bACC' = test balanced accuracy.",
        "Ranked by Test bACC mean over seeds 0/1/2. References (clinical-only / MRI-only / clinical-temporal-EN)",
        "and learned fusion methods are pooled; rows with identical metrics across folders are de-duplicated;",
        "rows scored on n < 30 TEST patients are EXCLUDED (e.g. baseline_only weighted-avg bACC=1.000 on n=2).",
        "MRI model throughout = BrainMVP (uniformer, full_ft, stochastic aug). MRI-only @m12 pinned as a reference.",
        "Class-wise rows (A_/B_*, from up_to_m12 class_wise_detector_fusion): per-class combiner over a class's",
        "  probabilities across bl/m06/m12 (+MRI), renormalised. A = clinical source is the T2 3-class marginal;",
        "  B = clinical source is the T1(->CN)/T1b(->AD) binary detectors. _LR/_EN = learned logreg / elastic-net",
        "  combiner; _meanP = PRESENT-ONLY per-class MEAN (no imputation, no weights). CL = clinical only; CL+MRI adds MRI.",
        "CAVEAT: cohort/n differs by timeframe — baseline base = CL_bl @bl; m12 / up_to_m12 base = CL_m12 @m12;",
        "complete-case subsets are smaller. Test bACC is NOT on an identical cohort across rows — read with n.",
        "Full ranked list (incl. excluded n<30) in summary_methods_all_ranked.csv.",
    ]
    ax.text(0.0, -0.04, "\n".join(foot), transform=ax.transAxes, ha="left", va="top",
            fontsize=7.6, family="monospace", linespacing=1.4)
    fig.savefig(out_path, dpi=180, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  PNG: {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
