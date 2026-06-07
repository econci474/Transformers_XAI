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

import argparse
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


def cleaned_name(row):
    """Compact thesis name (the MRI foundation model is shown in its own column, NOT in the name).
    The CL T2 3-class (formulation A) model is BioClinical-ModernBERT-large/full_ft -> 'BioClin-L-ft';
    formulation B uses the dedicated T1/T1b binary DETECTORS -> 'CL T1/T1b'."""
    m = str(row["method"]); tf = row["timeframe"]; fam = str(row.get("family", ""))
    mfr = re.match(r"^(M12X|M12_DET|M12|BL)\s*/\s*", m)
    framing = mfr.group(1) if mfr else ""
    core = m[mfr.end():] if mfr else m
    # clinical timepoint for the cross-sectional rows (CLbl_Mm12 = @bl, CLm12_Mm12/M12X = @m12)
    if framing == "BL" or "CLbl" in fam:
        cl_tp = "@bl"
    elif framing in ("M12X", "M12", "M12_DET") or "CLm12" in fam:
        cl_tp = "@m12"
    else:
        cl_tp = {"baseline_only": "@bl", "m12_only": "@m12"}.get(tf, "@m12")

    if "MRI-only" in core:
        return "MRI-only [T2 3-class @m12]"
    # cross-sectional m12 fusions (CLbl_Mm12 / CLm12_Mm12): CL T2 @tp fused with MRI T2 @m12
    mm2 = re.match(r"^(weighted-avg(?:-J)?|elastic-net|logreg)\s*/\s*(\w+)$", core)
    if mm2 and ("CLbl" in fam or "CLm12" in fam):
        meth, miss = mm2.groups()                      # missing-mode (complete_case) shown via the n column
        return f"BioClin-L-ft T2 {cl_tp} + MRI T2 @m12  [{meth}]"
    if "clinical-only" in core:
        return f"BioClin-L-ft T2 {cl_tp}"
    if "clinical-temporal-EN" in core:
        return "BioClin-L-ft T2 @bl+m06+m12  [EN]"
    # up_to_m12 class-wise detector fusion (A = T2 marginal, B = T1/T1b detectors)
    mm = re.match(r"^([AB])_(LR|EN|meanP)\s*/\s*(CL\+MRI|CL\+T1b|CL)$", core)
    if mm:
        form, comb, src = mm.groups()
        clin = "BioClin-L-ft T2 @bl+m06+m12" if form == "A" else "CL T1/T1b @bl+m06+m12"
        mri = {"CL": "", "CL+MRI": " + MRI T1/T1b @m12", "CL+T1b": " + MRI T1b @m12"}[src]
        wt = {"meanP": "", "LR": "  [logreg]", "EN": "  [EN]"}[comb]
        return f"{clin}{mri}{wt}"
    # up_to_m12 remainder-product detectors
    mm3 = re.match(r"^detector-(structured|EN)\s*/\s*(CL\+MRI|CL)$", core)
    if mm3:
        meth, src = mm3.groups()
        mri = " + MRI T1/T1b @m12" if src == "CL+MRI" else ""
        return f"CL T1/T1b detectors @bl+m06+m12{mri}  [remainder-{meth}]"
    # up_to_m12 flat / hierarchical temporal (clinical T2 + MRI T2 3-class)
    if fam in ("flat", "hierarchical"):
        return f"BioClin-L-ft T2 @bl+m06+m12 + MRI T2 @m12  [{fam}: {core}]"
    return core


def uses_mri(method):
    """Whether a method actually fuses MRI (controls the MRI column: model name vs '—')."""
    s = str(method)
    if "MRI-only" in s:
        return True
    if "clinical-only" in s or "clinical-temporal-EN" in s:
        return False
    if re.search(r"/\s*CL\+(MRI|T1b)\b", s):     # class-wise / remainder CL+MRI / CL+T1b
        return True
    if re.search(r"/\s*CL\b\s*$", s):            # CL-only source -> no MRI
        return False
    # cross-sectional CLbl/CLm12 + flat/hierarchical fusions all consume MRI
    return bool(re.search(r"(weighted-avg|elastic-net|logreg|avg)", s))


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
        # outputs are namespaced by MRI foundation model: {timeframe}/{BrainMVP|BrainDINO}/{family}/
        if len(rel) > 1 and rel[1] in ("BrainMVP", "BrainDINO"):
            mri_model, fam_parts = rel[1], rel[2:-1]
        else:
            mri_model, fam_parts = "", rel[1:-1]
        family = "/".join(fam_parts)
        where = timeframe + (f" / {mri_model}" if mri_model else "") + (f" / {family}" if family else "")
        d = pd.read_csv(f)
        if "bacc_mean" not in d.columns or "n_test_mean" not in d.columns:
            continue
        # exclude metric columns AND the val_* tuning columns (which otherwise leak into method labels)
        idcols = [c for c in d.columns if not any(
            c.startswith(p) for p in ("bacc", "macro", "f1_", "n_test", "val"))]
        for _, r in d.iterrows():
            rec = {"timeframe": timeframe, "mri_model": mri_model, "family": family, "where": where,
                   "method": _label(r, idcols), "source_file": str(Path(f).relative_to(OUT))}
            for m in METRIC_MEANS + ["bacc_std", "n_test_mean"]:
                rec[m] = float(r[m]) if m in d.columns and pd.notna(r[m]) else np.nan
            rows.append(rec)
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cleaned", action="store_true",
                    help="Write a cleaned thesis table summary_t2_integration_cleaned.{csv,png}: "
                         "house style (serif/bordered/ruled), NO footnote, CL model named in the "
                         "title with abbrev BioClin-L-ft in the rows; pins method-6 (T1/T1b) + the "
                         "new T1b-only temporal variant + the MRI-only ref. Leaves the existing "
                         "summary_top5_bACC / summary_methods_all_ranked outputs untouched.")
    args = ap.parse_args()

    df = collect()
    if df.empty:
        print("[error] no fusion_metrics.csv found under", OUT); return 1
    # name a bare MRI-only ref (m12_only CLbl/CLm12 'mri_only', no model token) as the T2 3-class head;
    # the namespaced rows already carry 'BrainMVP'/'BrainDINO' in their source label and mri_model column.
    mri_unnamed = (df["method"].str.contains("MRI-only", case=False)
                   & ~df["method"].str.contains("BrainMVP")
                   & ~df["method"].str.contains("BrainDINO"))
    df.loc[mri_unnamed, "method"] = df.loc[mri_unnamed, "method"] + " [T2 3-class @m12]"
    df["display"] = df.apply(pretty_method, axis=1)

    # dedup identical-metric rows (same underlying result reported in multiple folders); the MRI model
    # is part of the key so BrainMVP/BrainDINO rows never collapse into each other.
    keytuple = lambda r: (r["mri_model"],) + tuple(round(r[m], 4) if pd.notna(r[m]) else np.nan
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

    if args.cleaned:
        return render_cleaned(inc)

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
TOP_N_CLEANED = 10


def render_cleaned(inc) -> int:
    """Cleaned cross-model thesis table: top-N by bACC unioned with the rows the user wants visible,
    for BOTH MRI models (BrainMVP, BrainDINO): method-6 (A_meanP/CL+MRI), the T1b-only variant
    (A_meanP/CL+T1b), the CLbl_Mm12 best non-weighted row, and one MRI-only ref."""
    inc = inc.copy()
    inc["display"] = inc.apply(cleaned_name, axis=1)
    inc["mri_disp"] = [r["mri_model"] if (r["mri_model"] and uses_mri(r["method"])) else "—"
                       for _, r in inc.iterrows()]

    pins = [inc.head(TOP_N_CLEANED)]
    for model in ("BrainMVP", "BrainDINO"):
        mm = inc[inc.mri_model == model]
        pins.append(mm[(mm.timeframe == "up_to_m12") & (mm.method == "A_meanP / CL+MRI")].head(1))
        pins.append(mm[(mm.timeframe == "up_to_m12") & (mm.method == "A_meanP / CL+T1b")].head(1))
        pins.append(mm[mm.family.str.contains("CLbl_Mm12", na=False)
                       & ~mm.method.str.contains("weighted-avg", na=False)].head(1))
        pins.append(mm[mm.method.str.contains("MRI-only", case=False, na=False)].head(1))
    sel = (pd.concat(pins)
           .drop_duplicates(subset=["mri_model", "timeframe", "family", "method"])
           .sort_values("bacc_mean", ascending=False)
           # collapse model-agnostic (non-MRI, "—") rows that recur across both model namespaces
           .drop_duplicates(subset=["mri_disp", "display"], keep="first")
           .reset_index(drop=True))

    out_cols = (["rank", "display", "mri_disp", "timeframe", "family"]
                + METRIC_MEANS + ["bacc_std", "n_test_mean"])
    sel[out_cols].rename(columns={"display": "method", "mri_disp": "MRI"}).to_csv(
        OUT / "summary_t2_integration_cleaned.csv", index=False)
    render_house(sel, OUT / "summary_t2_integration_cleaned.png")

    print("=" * 90)
    print("  T2 cross-model summary — CLEANED (house style, BioClin-L-ft; BrainMVP vs BrainDINO)")
    print("=" * 90)
    print(sel[["rank", "mri_disp", "display", "bacc_mean", "bacc_std", "macro_auc_mean",
               "n_test_mean"]].to_string(index=False))
    print(f"\n  wrote -> {OUT/'summary_t2_integration_cleaned.png'}"
          f"\n         {OUT/'summary_t2_integration_cleaned.csv'}")
    return 0


def render_house(sel, out_path):
    """Manual house-style table (DejaVu Serif, bordered, ruled, italic headers, dashed dividers,
    bold best-per-column) — mirrors clinical_pipeline/01e_render_t4_stratification.py."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = "DejaVu Serif"

    COLS = [("bacc_mean", "Test bACC", "bacc_std"), ("macro_f1_mean", "macro-F1", None),
            ("macro_auc_mean", "macro-AUC", None), ("f1_CN_mean", "F1 CN", None),
            ("f1_MCI_mean", "F1 MCI", None), ("f1_AD_mean", "F1 AD", None)]
    headers = ["#", "Method", "MRI"] + [c[1] for c in COLS] + ["n"]
    N_LEAD = 3                                          # leading non-metric cols: #, Method, MRI

    body, numeric, is_ref = [], [], []
    for _, r in sel.iterrows():
        cells = [str(int(r["rank"])), r["display"], str(r.get("mri_disp", "—"))]
        nums = []
        for key, _, stdk in COLS:
            v = r[key]
            if pd.isna(v):
                cells.append("—"); nums.append(np.nan)
            else:
                s = f"{v:.3f}" + (f" ± {r[stdk]:.3f}" if stdk and pd.notna(r.get(stdk)) else "")
                cells.append(s); nums.append(float(v))
        cells.append(str(int(round(r["n_test_mean"]))))
        body.append(cells); numeric.append(nums)
        is_ref.append(str(r["display"]).startswith("MRI-only"))
    numeric = np.array(numeric, float)
    best = [int(np.nanargmax(numeric[:, j])) if not np.all(np.isnan(numeric[:, j])) else -1
            for j in range(len(COLS))]

    # column widths (inches): #, Method (wide, left), MRI, 6 metric cols, n
    COL_W = [0.45, 4.95, 1.10, 1.45, 1.10, 1.20, 0.90, 0.96, 0.90, 0.50]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H = 1.06, 0.40, 0.44
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
    ax.text((LEFT + RIGHT) / 2, (y_title_top + y) / 2,
            "T2 Late Fusion: Best Methods Across Timeframes & MRI Models\n"
            "mean ± std across seeds 0/1/2 · ranked by Test balanced accuracy\n"
            "clinical T2 = BioClinical-ModernBERT-large (full fine-tune, \"BioClin-L-ft\")  ·  "
            "MRI = BrainMVP or BrainDINO (frozen)",
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    hline(y_title_top, lw=1.5)
    hline(y, lw=1.2)

    # header row (italic): #, Method left-aligned; the rest centered
    y_head_top = y
    y -= HEAD_H
    ymid = (y_head_top + y) / 2
    ax.text(col_left[0] + 0.06, ymid, headers[0], ha="left", va="center",
            fontsize=9.5, fontstyle="italic")
    ax.text(col_left[1] + 0.06, ymid, headers[1], ha="left", va="center",
            fontsize=9.5, fontstyle="italic")
    for j in range(2, len(headers)):
        ax.text(col_cx[j], ymid, headers[j], ha="center", va="center",
                fontsize=9.5, fontstyle="italic")
    hline(y, lw=1.2)

    y_data_top = y
    for i, cells in enumerate(body):
        if is_ref[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))      # dashed rule before the MRI-only reference
        yr_top = y
        y -= ROW_H
        ymid = (yr_top + y) / 2
        ax.text(col_left[0] + 0.06, ymid, cells[0], ha="left", va="center", fontsize=9.0)
        ax.text(col_left[1] + 0.06, ymid, cells[1], ha="left", va="center", fontsize=9.0)
        for j in range(2, len(headers)):
            metric_j = j - N_LEAD                  # metric columns start after #, Method, MRI
            bold = (0 <= metric_j < len(COLS) and best[metric_j] == i)
            ax.text(col_cx[j], ymid, cells[j], ha="center", va="center", fontsize=9.0,
                    fontweight="bold" if bold else "normal")
    BOTTOM = y

    # dashed vertical dividers between every column
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
