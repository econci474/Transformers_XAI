"""
03_summary_T1e_cleaned.py
=========================
Cleaned house-style copy of the T1e (sCN vs pCN) late-fusion LEADERBOARD
(`outputs/baseline/fusion_table_T1e_summary.png`), matching the T1d / T2 cleaned tables
(DejaVu-Serif, bordered, ruled, italic headers, dashed dividers, bold best-per-column, clean
bold title + italic subtitle). No fusion re-fit — reads the existing per-combo CSVs.

Differences vs the raw leaderboard (per user):
  - "Tuned params" column removed from the main table -> moved to a separate `_hp` table.
  - AUC and macro-F1 swapped (-> Test bACC, AUC, macro-F1 ... matching T1d/T2 ordering).
  - House multi-line bold title + italic N(test)= subtitle with the sCN/pCN class breakdown.
  - "avg" -> "mean"; same model / aggregation nomenclature as the T1d/T2 cleaned tables.
  - "Rank" column dropped; references grouped below dotted rules (house style).

The leaderboard = top-6 meta_prs fusion rows + 2 CL-only refs + 1 SNP-only ref (9 rows).

Out: integration/T1e_classification/outputs/baseline/summary_t1e_cleaned.{csv,png,pdf} (+ _no_n)
     integration/T1e_classification/outputs/baseline/summary_t1e_hp.{csv,png,pdf}
Run: python integration/T1e_classification/03_summary_T1e_cleaned.py
"""
import importlib.util
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Serif"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs", "baseline")

# reuse the fusion harness's loader + tuned-param formatter (module import has no side effects)
spec = importlib.util.spec_from_file_location("t1e_fuse", os.path.join(HERE, "01_fuse_T1e.py"))
t1e = importlib.util.module_from_spec(spec); spec.loader.exec_module(t1e)

SEEDS = [0, 1, 2]
FUSION_METHODS = ["equal_avg", "class_avg", "best_alpha", "best_alpha_J", "stacked_en"]

CL_ABBR = {"bioclin_mb_base_frozen": "BioClin-B-frozen", "modernbert_base_full_ft": "MBERT-B-ft"}
SNP_ABBR = {"meta_prs": "meta-PRS (EN)", "kosteridis_prs": "Kosteridis PRS",
            "meta_prs_filtered": "meta-PRS (filtered)",
            "kosteridis_shared_AD_CV": "Kosteridis (shared AD-CV)",
            "prscs_bellenguez": "Bellenguez (PRS-CS)",
            "prscs_derojas":    "DeRojas (PRS-CS)",
            "prscs_kunkle":     "Kunkle (PRS-CS)",
            "prs_all_dedup":     "all-dedup PRS",
            "prs_all_dedup_ivw": "all-dedup PRS (IVW)",
            "bmfm_none_mlp": "BMFM-SNP (none·MLP)", "ntv2_permod_mlp": "NTv2 (per-mod·MLP)"}
# model-type ordering for the expanded table: PRS block then FM block (dotted divider between)
SNP_TYPE = {"meta_prs": "PRS", "kosteridis_prs": "PRS",
            "meta_prs_filtered": "PRS", "kosteridis_shared_AD_CV": "PRS",
            "prscs_bellenguez": "PRS", "prscs_derojas": "PRS", "prscs_kunkle": "PRS",
            "prs_all_dedup": "PRS", "prs_all_dedup_ivw": "PRS",
            "bmfm_none_mlp": "FM", "ntv2_permod_mlp": "FM"}
SNP_TYPE_ORDER = {"PRS": 0, "FM": 1}
AGG_LABEL = {"equal_avg": "mean", "best_alpha": "weighted-avg (val bACC)",
             "best_alpha_J": "weighted-avg (Youden J)", "class_avg": "class-wise (val bACC)",
             "stacked_en": "elastic-net", "clinical_only": "—", "snp_only": "—"}

CLEAN_METRICS = [("bacc", "Test bACC"), ("auc", "AUC"), ("macro_f1", "macro-F1"),
                 ("f1_sCN", "F1 sCN"), ("f1_pCN", "F1 pCN"),
                 ("R2_lia_null", "R^2 liab|null"),
                 ("brier", "Brier")]

# K_p for Lee 2012 liability conversion (population AD prevalence; matches
# the corrected-GTEx master leaderboard's val_R2_lia_null column).
K_POP = 0.283


def _cox_snell_R2(y, p):
    """Cox-Snell R² = 1 - exp(2(LL_null - LL_fit)/n) with LL_null at sample prevalence."""
    from scipy import stats as _sst   # noqa: F401  (kept local to function scope)
    eps = 1e-10
    p = np.clip(p, eps, 1 - eps)
    n = len(y); ybar = float(np.mean(y))
    if n < 5 or ybar in (0.0, 1.0):
        return float("nan")
    LL_null = np.sum(y * np.log(ybar) + (1 - y) * np.log(1 - ybar))
    LL_fit  = np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
    return float(1 - np.exp(2 * (LL_null - LL_fit) / n))


def _lee_liability_R2(R2cs, P, K=K_POP):
    """Lee 2012 eq.8 liability-scale conversion from observed-scale R²_CS."""
    from scipy import stats as _sst
    if np.isnan(R2cs) or P in (0.0, 1.0):
        return float("nan")
    t = _sst.norm.ppf(1 - K); phi_t = _sst.norm.pdf(t)
    if phi_t <= 0 or P * (1 - P) <= 0:
        return float("nan")
    return float(R2cs * (K * (1 - K)) ** 2 / (P * (1 - P) * phi_t ** 2))


_R2_CACHE = None


def _brier_score(y, p):
    """Brier score = mean((p - y)^2). Lower = better. Range [0, 1]."""
    eps = 1e-12
    p = np.clip(p, eps, 1 - eps)
    return float(np.mean((p - y) ** 2))


def get_r2_liab_map():
    """Lazy-compute and cache R²_liab over null AND Brier per (clin_v, snp_v, method).
    Keys:
      ('clin_v', 'snp_v', 'method')         for fusion methods
      ('clin_v', None,    'clinical_only')  for clinical-only references
      (None,     'snp_v', 'snp_only')        for SNP-only references
    Returns dict[key] = {'R2': (mean, std), 'brier': (mean, std)}."""
    global _R2_CACHE
    if _R2_CACHE is not None:
        return _R2_CACHE
    cache = {}
    fused = pd.read_csv(os.path.join(OUT_DIR, "fused_predictions.csv"))
    def _agg(r2vals, brvals, key):
        cache[key] = {
            "R2": (float(np.nanmean(r2vals)),
                    float(np.nanstd(r2vals, ddof=1)) if len(r2vals) > 1 else 0.0),
            "brier": (float(np.nanmean(brvals)),
                       float(np.nanstd(brvals, ddof=1)) if len(brvals) > 1 else 0.0),
        }
    # Fusion methods
    for (cv, sv, mt), grp in fused.groupby(["clin_variant", "snp_variant", "method"]):
        r2vals, brvals = [], []
        for _, sgrp in grp.groupby("seed"):
            y = sgrp["y_true"].astype(int).values
            p = sgrp["prob_pCN"].astype(float).values
            r2cs = _cox_snell_R2(y, p)
            r2vals.append(_lee_liability_R2(r2cs, float(np.mean(y))))
            brvals.append(_brier_score(y, p))
        if r2vals:
            _agg(r2vals, brvals, (cv, sv, mt))
    # Clinical-only: one per CL variant
    for cl in t1e.CLIN_VARIANTS:
        r2vals, brvals = [], []
        for seed in SEEDS:
            cl_df = t1e.load_clinical(seed, cl["model"], cl["strat"])
            cl_te = cl_df[cl_df["split"] == "test"]
            y = cl_te["y_clin"].astype(int).values
            p = cl_te["cp1"].astype(float).values
            r2cs = _cox_snell_R2(y, p)
            r2vals.append(_lee_liability_R2(r2cs, float(np.mean(y))))
            brvals.append(_brier_score(y, p))
        if r2vals:
            _agg(r2vals, brvals, (cl["key"], None, "clinical_only"))
    # SNP-only: one per SNP variant
    for sn in t1e.SNP_VARIANTS:
        r2vals, brvals = [], []
        for seed in SEEDS:
            try:
                sn_df = t1e.load_snp(seed, sn["parquet_tmpl"])
            except FileNotFoundError:
                continue
            sn_te = sn_df[sn_df["split"] == "test"]
            if sn_te.empty:
                continue
            y = sn_te["y_snp"].astype(int).values
            p = sn_te["sp1"].astype(float).values
            r2cs = _cox_snell_R2(y, p)
            r2vals.append(_lee_liability_R2(r2cs, float(np.mean(y))))
            brvals.append(_brier_score(y, p))
        if r2vals:
            _agg(r2vals, brvals, (None, sn["key"], "snp_only"))
    _R2_CACHE = cache
    print(f"[R2/Brier] cached for {len(cache)} (CL × SNP × method) keys")
    return cache


# --------------------------------------------------------------------------- #
# Row assembly
# --------------------------------------------------------------------------- #
def class_breakdown():
    """Per-seed mean (sCN, pCN) test counts, from the clinical loader (label 0=sCN, 1=pCN)."""
    cl = t1e.CLIN_VARIANTS[0]
    counts = []
    for seed in SEEDS:
        d = t1e.load_clinical(seed, cl["model"], cl["strat"])
        t = d[d["split"] == "test"]
        y = t["y_clin"].astype(int).to_numpy()
        counts.append([(y == 0).sum(), (y == 1).sum()])
    a = np.array(counts, float).mean(0)
    return int(round(a[0])), int(round(a[1]))


def tuned_param_map():
    """(clin_variant, snp_variant, method) -> tuned-param string, from per-seed means."""
    per = pd.read_csv(os.path.join(OUT_DIR, "fusion_metrics_per_seed.csv"))
    per = per[per["split"] == "test"]
    g = ["clin_variant", "snp_variant", "method"]
    cols = ["alpha", "threshold", "w0", "w1", "best_C", "best_l1_ratio"]
    pmean = per.groupby(g)[cols].mean().round(3).reset_index()
    out = {}
    for _, r in pmean.iterrows():
        row = {f"{c}_mean": r[c] for c in cols}
        out[(r.clin_variant, r.snp_variant, r.method)] = t1e._tuned_params_str(r.method, pd.Series(row))
    return out


def _type_key(snp_variant):
    return SNP_TYPE_ORDER.get(SNP_TYPE.get(snp_variant, "PRS"), 9)


def assemble(snp_fusion_keys, top_n, snp_ref_keys, group_by_type=False, clin_filter=None):
    """Build the row table. snp_fusion_keys = SNP sources whose fusion combos enter the ranking;
    top_n = keep this many top fusion rows (None -> all); snp_ref_keys = SNP-only reference rows.
    group_by_type: order fusion + SNP-only rows by model type (PRS block, then FM block) with a dotted
    divider at the type transition, instead of a single Test-bACC ranking.
    clin_filter: restrict fusion + clinical-only refs to this CL variant (SNP-only refs are CL-independent
    so they are kept)."""
    fm = pd.read_csv(os.path.join(OUT_DIR, "fusion_metrics.csv"))
    test = fm[fm["split"] == "test"].copy()
    if clin_filter is not None:
        test = test[(test.clin_variant == clin_filter) | (test.method == "snp_only")].copy()

    fus = test[(test.snp_variant.isin(snp_fusion_keys)) & (test.method.isin(FUSION_METHODS))].copy()
    snpref = (test[(test.method == "snp_only") & (test.snp_variant.isin(snp_ref_keys))]
              .drop_duplicates("snp_variant").copy())
    if group_by_type:
        for d in (fus, snpref):
            d["_t"] = d.snp_variant.map(_type_key)
        fus = fus.sort_values(["_t", "bacc_mean"], ascending=[True, False], kind="stable")
        snpref = snpref.sort_values(["_t", "bacc_mean"], ascending=[True, False], kind="stable")
    else:
        fus = fus.sort_values("bacc_mean", ascending=False, kind="stable")
        snpref = snpref.sort_values("bacc_mean", ascending=False, kind="stable")
    if top_n is not None:
        fus = fus.head(top_n)
    clref = (test[test.method == "clinical_only"].drop_duplicates("clin_variant")
             .sort_values("bacc_mean", ascending=False, kind="stable"))

    recs = []
    for blk, sub in (("fusion", fus), ("cl", clref), ("snp", snpref)):
        for _, r in sub.iterrows():
            if r.method == "snp_only":
                method = "—"                                   # SNP source already named in the SNP column
                snp = SNP_ABBR[r.snp_variant]
            else:
                method = CL_ABBR[r.clin_variant]
                snp = "—" if r.method == "clinical_only" else SNP_ABBR[r.snp_variant]
            rec = {"block": blk, "type": SNP_TYPE.get(r.snp_variant, "—"),
                   "clin_variant": r.clin_variant, "snp_variant": r.snp_variant,
                   "raw_method": r.method, "Method": method, "SNP": snp,
                   "Aggregation": AGG_LABEL.get(r.method, r.method), "n": r.n_mean}
            for key, _ in CLEAN_METRICS:
                if key in ("R2_lia_null", "brier"):
                    continue  # filled from cache below
                rec[f"{key}_m"] = r[f"{key}_mean"]; rec[f"{key}_s"] = r[f"{key}_std"]
            # R²_liab over null AND Brier: look up by row type
            r2cache = get_r2_liab_map()
            if r.method == "clinical_only":
                k = (r.clin_variant, None, "clinical_only")
            elif r.method == "snp_only":
                k = (None, r.snp_variant, "snp_only")
            else:
                k = (r.clin_variant, r.snp_variant, r.method)
            entry = r2cache.get(k, {"R2": (float("nan"), float("nan")),
                                       "brier": (float("nan"), float("nan"))})
            rec["R2_lia_null_m"], rec["R2_lia_null_s"] = entry["R2"]
            rec["brier_m"], rec["brier_s"] = entry["brier"]
            recs.append(rec)
    df = pd.DataFrame(recs)
    # dotted rule above each block change, and (if grouping) above each model-type change within
    # the fusion / snp blocks.
    df["rule_above"] = False
    prev_blk, prev_type = None, None
    for i in df.index:
        blk, typ = df.at[i, "block"], df.at[i, "type"]
        change = (prev_blk is not None and blk != prev_blk)
        if group_by_type and blk in ("fusion", "snp") and blk == prev_blk and typ != prev_type:
            change = True
        df.at[i, "rule_above"] = bool(change)
        prev_blk, prev_type = blk, typ
    return df


# --------------------------------------------------------------------------- #
# Renderers (house style; mirrors T1d_classification/02_top5_summary.py::render_house)
# --------------------------------------------------------------------------- #
_CL_LINE = ('clinical = ModernBERT-base ("MBERT-B-ft", full fine-tune) or '
            'BioClinical-ModernBERT-base ("BioClin-B-frozen", frozen)')
TITLE_MAIN = ('T1e Late Fusion (sCN vs pCN)\n'
              'mean ± std across seeds 0/1/2 · ranked by Test balanced accuracy\n'
              + _CL_LINE + '\n'
              'SNP = meta-PRS (elastic-net combined)')
_SNP_LINE_EXPANDED = ('SNP = PRS (meta-PRS EN/filtered, Kosteridis full/shared-AD-CV, '
                      'Bellenguez/DeRojas/Kunkle PRS-CS) or FM (BMFM-SNP none·MLP, NTv2 per-mod·MLP)')
TITLE_EXPANDED = ('T1e Late Fusion (sCN vs pCN)\n'
                  'mean ± std across seeds 0/1/2 · grouped by SNP model type, Test-bACC-sorted within\n'
                  + _CL_LINE + '\n'
                  + _SNP_LINE_EXPANDED)
TITLE_MBERT = ('T1e Late Fusion (sCN vs pCN)\n'
               'mean ± std across seeds 0/1/2 · grouped by SNP model type, Test-bACC-sorted within\n'
               'clinical = ModernBERT-base ("MBERT-B-ft", full fine-tune)\n'
               + _SNP_LINE_EXPANDED)
TITLE_MBERT_FILTERED = ('T1e Late Fusion (sCN vs pCN) — best fusion method per (CL × SNP)\n'
                        'mean ± std across seeds 0/1/2 · grouped by SNP model type, Test-bACC-sorted within\n'
                        'clinical = ModernBERT-base ("MBERT-B-ft", full fine-tune)\n'
                        + _SNP_LINE_EXPANDED)
_SNP_LINE_SMALL = ('SNP = PRS (meta-PRS EN/filtered, Kosteridis, Kunkle PRS-CS, all-dedup PRS IVW) '
                   'or FM (BMFM-SNP none·MLP, NTv2 per-mod·MLP)')
TITLE_MBERT_SMALL = ('T1e Late Fusion (sCN vs pCN) — best fusion method per (CL × SNP)\n'
                     'mean ± std across seeds 0/1/2 · grouped by SNP model type, Test-bACC-sorted within\n'
                     'clinical = ModernBERT-base ("MBERT-B-ft", full fine-tune)\n'
                     + _SNP_LINE_SMALL)


def hline_fn(ax, LEFT, RIGHT):
    def hline(yy, lw=1.0, ls="-"):
        ax.plot([LEFT, RIGHT], [yy, yy], color="black", linewidth=lw, linestyle=ls,
                solid_capstyle="butt", zorder=3)
    return hline


def render_main(df, out_path, subtitle, title=TITLE_MAIN, show_n=True, metrics=None):
    metrics = metrics or CLEAN_METRICS
    headers = ["Clinical", "SNP", "Aggregation"] + [lbl for _, lbl in metrics] + (["n"] if show_n else [])
    N_LEAD = 3
    LEFT_COLS = {0, 2}

    body, numeric, rules = [], [], []
    for _, r in df.iterrows():
        cells = [r["Method"], r["SNP"], r["Aggregation"]]
        nums = []
        for key, _ in metrics:
            v, s = r.get(f"{key}_m"), r.get(f"{key}_s")
            if pd.isna(v):
                cells.append("—"); nums.append(np.nan)
            else:
                cells.append(f"{v:.3f}" + (f" ± {s:.3f}" if pd.notna(s) else "")); nums.append(float(v))
        if show_n:
            cells.append(f"{r['n']:.0f}")
        body.append(cells); numeric.append(nums); rules.append(bool(r["rule_above"]))
    numeric = np.array(numeric, float)
    # For Brier: lower is better → use argmin; everything else: higher is better → argmax.
    _LOWER_IS_BETTER = {"brier"}
    best = []
    for j, (k, _) in enumerate(metrics):
        col = numeric[:, j]
        if np.all(np.isnan(col)):
            best.append(-1)
        elif k in _LOWER_IS_BETTER:
            best.append(int(np.nanargmin(col)))
        else:
            best.append(int(np.nanargmax(col)))

    COL_W = [3.70, 2.45, 3.25] + [1.58] + [1.55] * (len(metrics) - 1) + ([0.55] if show_n else [])
    LEFT, RIGHT_PAD = 0.28, 0.28
    SUB_H = 0.40 if subtitle else 0.0
    TITLE_H, HEAD_H, ROW_H = 1.30 + SUB_H, 0.40, 0.44
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
    hline = hline_fn(ax, LEFT, RIGHT)

    y = fig_h - TOP_PAD
    y_title_top = y
    y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + SUB_H + (TITLE_H - SUB_H) / 2, title,
            ha="center", va="center", fontsize=10.5, fontweight="bold", linespacing=1.5)
    if subtitle:
        ax.text(cx, y + SUB_H / 2, subtitle, ha="center", va="center",
                fontsize=8.5, fontstyle="italic")
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
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y; y -= ROW_H
        ymid = (yr_top + y) / 2
        for j in range(len(headers)):
            metric_j = j - N_LEAD
            bold = (0 <= metric_j < len(metrics) and best[metric_j] == i)
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
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    fig.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  PNG: {out_path}")


def render_hp(df, hp_map, out_path, title="T1e Late Fusion — Tuned Operating Points"):
    headers = ["Clinical", "SNP", "Aggregation", "Tuned params"]
    LEFT_COLS = {0, 2, 3}
    body, rules = [], []
    for _, r in df.iterrows():
        tp = hp_map.get((r.clin_variant, r.snp_variant, r.raw_method), "—")
        body.append([r["Method"], r["SNP"], r["Aggregation"], tp]); rules.append(bool(r["rule_above"]))

    COL_W = [3.70, 1.95, 3.05, 2.50]
    LEFT, RIGHT_PAD = 0.28, 0.28
    TITLE_H, HEAD_H, ROW_H = 0.95, 0.40, 0.44
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
    hline = hline_fn(ax, LEFT, RIGHT)

    y = fig_h - TOP_PAD
    y_title_top = y
    y -= TITLE_H
    cx = (LEFT + RIGHT) / 2
    ax.text(cx, y + TITLE_H * 0.66, title,
            ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(cx, y + TITLE_H * 0.24,
            "operating point chosen on VAL, mean across seeds 0/1/2 (same row order as the leaderboard)\n"
            "best-α = α on VAL bACC (t=0.5) · weighted-avg (Youden J) = (α, t) on VAL Youden's J · "
            "class-wise = (w₀, w₁) on VAL bACC · elastic-net = (C, l1) CV-tuned",
            ha="center", va="center", fontsize=7.6, fontstyle="italic", linespacing=1.4)
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
        if rules[i]:
            hline(y, lw=0.6, ls=(0, (3, 3)))
        yr_top = y; y -= ROW_H
        ymid = (yr_top + y) / 2
        for j in range(len(headers)):
            if j in LEFT_COLS:
                ax.text(col_left[j] + 0.06, ymid, cells[j], ha="left", va="center", fontsize=9.0)
            else:
                ax.text(col_cx[j], ymid, cells[j], ha="center", va="center", fontsize=9.0)
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


def _write_csv(df, stem):
    metric_cols = [f"{k}_m" for k, _ in CLEAN_METRICS] + [f"{k}_s" for k, _ in CLEAN_METRICS]
    df[["Method", "SNP", "Aggregation"] + metric_cols + ["n"]].rename(
        columns={"Method": "Clinical"}).to_csv(
        os.path.join(OUT_DIR, f"{stem}.csv"), index=False, encoding="utf-8")


def _write_hp_csv(df, hp_map, stem):
    hp_df = df[["Method", "SNP", "Aggregation"]].copy()
    hp_df["Tuned params"] = [hp_map.get((r.clin_variant, r.snp_variant, r.raw_method), "—")
                             for _, r in df.iterrows()]
    hp_df.rename(columns={"Method": "Clinical"}).to_csv(
        os.path.join(OUT_DIR, f"{stem}.csv"), index=False, encoding="utf-8")


def write_t1e_latex(df, out_path, drop_perclass=False, label="tab:t1e_late_fusion", caption_note="",
                    snp_desc="SNP = meta-PRS (elastic-net combined)"):
    """LaTeX (sidewaystable, booktabs, multicolumn title + group labels, bottom caption with N(test) +
    abbreviations) — int_t2_simple style. drop_perclass=True drops the F1 sCN / F1 pCN columns;
    intra-block dashed rules (e.g. PRS/FM type splits) become \\midrule."""
    def esc(s):
        return (str(s).replace("—", "--").replace("·", r"$\cdot$").replace("&", r"\&")
                .replace("%", r"\%").replace("_", r"\_").replace("#", r"\#"))
    metric_cols = [("bacc", "Test bACC"), ("auc", "AUC"), ("macro_f1", "macro-F1")]
    if not drop_perclass:
        metric_cols += [("f1_sCN", "F1 sCN"), ("f1_pCN", "F1 pCN")]
    metric_cols += [("R2_lia_null", r"$R^2_{\mathrm{liab}|\mathrm{null}}$"),
                     ("brier", r"Brier")]
    ncol = 3 + len(metric_cols)
    mat = np.array([[r[f"{k}_m"] for k, _ in metric_cols] for _, r in df.iterrows()], float)
    _LOWER_IS_BETTER = {"brier"}
    best = []
    for j, (k, _) in enumerate(metric_cols):
        col = mat[:, j]
        if np.all(np.isnan(col)):
            best.append(-1)
        elif k in _LOWER_IS_BETTER:
            best.append(int(np.nanargmin(col)))
        else:
            best.append(int(np.nanargmax(col)))
    BLOCK = {"fusion": r"Clinical $\oplus$ SNP late fusion", "cl": "Clinical-only references",
             "snp": "SNP-only references"}

    def cell(m, s, bold):
        if pd.isna(m):
            return "--"
        txt = f"{m:.2f} \\pm {s:.2f}" if pd.notna(s) else f"{m:.2f}"
        return f"$\\mathbf{{{txt}}}$" if bold else f"${txt}$"

    head = [r"\textbf{Clinical}", r"\textbf{SNP}", r"\textbf{Aggregation}"] \
        + [r"\textbf{" + lbl + "}" for _, lbl in metric_cols]
    L = [r"% T1e Late Fusion -- clinical + SNP late fusion." + (" Filtered." if drop_perclass else ""),
         r"% Requires: \usepackage{rotating, booktabs, graphicx}",
         r"\begin{sidewaystable}[ht]", r"\centering",
         r"\resizebox{\textwidth}{!}{%", r"\normalsize",
         r"\begin{tabular}{lll " + "c" * len(metric_cols) + "}", r"\toprule",
         r"\multicolumn{" + str(ncol) + r"}{c}{\textbf{T1e Late Fusion (sCN vs pCN)}} \\", r"\midrule",
         " & ".join(head) + r" \\", r"\midrule"]
    prev = None
    for i, (_, r) in enumerate(df.iterrows()):
        b = r["block"]
        if b != prev:
            if prev is not None:
                L.append(r"\midrule")
            L.append(r"\multicolumn{" + str(ncol) + r"}{l}{\textit{" + BLOCK[b] + r"}} \\")
        elif bool(r["rule_above"]):                        # intra-block type split (PRS/FM)
            L.append(r"\midrule")
        cells = [esc(r["Method"]), esc(r["SNP"]), esc(r["Aggregation"])]
        cells += [cell(r[f"{k}_m"], r.get(f"{k}_s"), best[j] == i) for j, (k, _) in enumerate(metric_cols)]
        L.append(" & ".join(cells) + r" \\")
        prev = b
    cap = (r"T1e late fusion (sCN vs pCN), mean $\pm$ std across seeds 0/1/2, ranked by Test balanced "
           r"accuracy; fusion rows on top, single-modality references below. \textit{Clinical} = clinical "
           r"encoder; \textit{SNP} = polygenic-risk source; \textit{Aggregation} = the late-fusion "
           r"combiner. clinical = ModernBERT-base (``MBERT-B-ft'', full fine-tune) or "
           r"BioClinical-ModernBERT-base (``BioClin-B-frozen'', frozen); " + snp_desc + r". "
           r"N(test) = 22 per seed (sCN 18 / pCN 4); CL and SNP both present for all test "
           r"subjects. Bold = best per column." + caption_note)
    L += [r"\bottomrule", r"\end{tabular}%", r"}", r"\caption{" + cap + r"}",
          r"\label{" + label + r"}", r"\end{sidewaystable}"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  TEX: {out_path}")


def main():
    n_scn, n_pcn = class_breakdown()
    subtitle = (f"N(test) = 22 per seed (sCN {n_scn} / pCN {n_pcn})  ·  "
                "CL and SNP both present for all test subjects")
    hp_map = tuned_param_map()

    # --- main table: meta-PRS only (top-6 fusion + CL refs + meta-PRS SNP-only) ---
    main_df = assemble(snp_fusion_keys=["meta_prs"], top_n=6, snp_ref_keys=["meta_prs"])
    _write_csv(main_df, "summary_t1e_cleaned")
    _write_hp_csv(main_df, hp_map, "summary_t1e_hp")
    render_main(main_df, os.path.join(OUT_DIR, "summary_t1e_cleaned.png"), subtitle,
                title=TITLE_MAIN, show_n=True)
    render_main(main_df, os.path.join(OUT_DIR, "summary_t1e_cleaned_no_n.png"), subtitle,
                title=TITLE_MAIN, show_n=False)
    render_hp(main_df, hp_map, os.path.join(OUT_DIR, "summary_t1e_hp.png"))
    write_t1e_latex(main_df, os.path.join(OUT_DIR, "summary_t1e_cleaned.tex"), drop_perclass=True)

    # --- int_t1e_filtered: the headline fusion model + its corresponding unimodal references
    #     (clinical-only of the SAME CL model + the SNP-only ref); drops the F1 sCN / F1 pCN columns. ---
    THREE = [("bacc", "Test bACC"), ("auc", "AUC"), ("macro_f1", "macro-F1")]
    head_row = main_df[main_df.block == "fusion"].sort_values("bacc_m", ascending=False).head(1)
    head_cl = head_row.iloc[0]["Method"]
    clin_ref = main_df[(main_df.block == "cl") & (main_df.Method == head_cl)]
    snp_ref = main_df[main_df.block == "snp"]
    flt = pd.concat([head_row, clin_ref, snp_ref], ignore_index=True)
    flt["rule_above"] = [False] + [flt.block.iloc[i] != flt.block.iloc[i - 1] for i in range(1, len(flt))]
    flt[["Method", "SNP", "Aggregation"] + [f"{k}_m" for k, _ in THREE]
        + [f"{k}_s" for k, _ in THREE] + ["n"]].rename(columns={"Method": "Clinical"}).to_csv(
        os.path.join(OUT_DIR, "int_t1e_filtered.csv"), index=False, encoding="utf-8")
    render_main(flt, os.path.join(OUT_DIR, "int_t1e_filtered.png"), subtitle, title=TITLE_MAIN,
                show_n=True, metrics=THREE)
    render_main(flt, os.path.join(OUT_DIR, "int_t1e_filtered_no_n.png"), subtitle, title=TITLE_MAIN,
                show_n=False, metrics=THREE)
    write_t1e_latex(flt, os.path.join(OUT_DIR, "int_t1e_filtered.tex"), drop_perclass=True,
                    label="tab:t1e_late_fusion_min",
                    caption_note=r" Headline fusion vs its unimodal (clinical-only and SNP-only) "
                                 r"references; full results in Table~\ref{tab:t1e_late_fusion}.")

    # --- expanded table: PRS (meta-PRS EN/filtered, Kosteridis full/shared, PRS-CS Bellenguez/
    #     DeRojas/Kunkle) + FM (BMFM-none, NTv2-per-mod) fusion combos + matching SNP-only refs,
    #     grouped by model type (PRS block then FM block) ---
    snp_all = ["meta_prs", "meta_prs_filtered",
               "kosteridis_prs", "kosteridis_shared_AD_CV",
               "prscs_bellenguez", "prscs_derojas", "prscs_kunkle",
               "prs_all_dedup", "prs_all_dedup_ivw",
               "bmfm_none_mlp", "ntv2_permod_mlp"]
    exp_df = assemble(snp_fusion_keys=snp_all, top_n=None, snp_ref_keys=snp_all, group_by_type=True)
    _write_csv(exp_df, "summary_t1e_expanded")
    _write_hp_csv(exp_df, hp_map, "summary_t1e_expanded_hp")
    render_main(exp_df, os.path.join(OUT_DIR, "summary_t1e_expanded.png"), subtitle,
                title=TITLE_EXPANDED, show_n=True)
    render_main(exp_df, os.path.join(OUT_DIR, "summary_t1e_expanded_no_n.png"), subtitle,
                title=TITLE_EXPANDED, show_n=False)
    render_hp(exp_df, hp_map, os.path.join(OUT_DIR, "summary_t1e_expanded_hp.png"),
              title="T1e Late Fusion (expanded) — Tuned Operating Points")
    write_t1e_latex(exp_df, os.path.join(OUT_DIR, "summary_t1e_expanded.tex"), drop_perclass=True,
                    label="tab:t1e_late_fusion_expanded",
                    snp_desc=(r"SNP = PRS (meta-PRS EN/filtered, Kosteridis full/shared-AD-CV, "
                              r"Bellenguez/DeRojas/Kunkle PRS-CS) or FM (BMFM-SNP none$\cdot$MLP, "
                              r"NTv2 per-mod$\cdot$MLP), grouped by SNP model type (PRS block then FM "
                              r"block, Test-bACC-sorted within)"))

    # --- ModernBERT-B-ft-only expanded variant (drops BioClin-B-frozen fusion + CL ref) ---
    mbert_df = assemble(snp_fusion_keys=snp_all, top_n=None, snp_ref_keys=snp_all,
                        group_by_type=True, clin_filter="modernbert_base_full_ft")
    _write_csv(mbert_df, "summary_t1e_modernbert")
    render_main(mbert_df, os.path.join(OUT_DIR, "summary_t1e_modernbert.png"), subtitle,
                title=TITLE_MBERT, show_n=True)
    render_main(mbert_df, os.path.join(OUT_DIR, "summary_t1e_modernbert_no_n.png"), subtitle,
                title=TITLE_MBERT, show_n=False)

    # --- ModernBERT-B-ft FILTERED: best fusion method per (CL, SNP) combination ---
    fus = mbert_df[mbert_df.block == "fusion"].copy()
    fus_best = (fus.sort_values("bacc_m", ascending=False)
                   .drop_duplicates(["clin_variant", "snp_variant"], keep="first"))
    fus_best["_t"] = fus_best.snp_variant.map(_type_key)
    fus_best = fus_best.sort_values(["_t", "bacc_m"], ascending=[True, False],
                                      kind="stable").drop(columns="_t")
    mbert_filt = pd.concat([fus_best,
                              mbert_df[mbert_df.block == "cl"],
                              mbert_df[mbert_df.block == "snp"]], ignore_index=True)
    # Recompute rule_above for block/type transitions
    mbert_filt["rule_above"] = False
    prev_blk, prev_type = None, None
    for i in mbert_filt.index:
        blk, typ = mbert_filt.at[i, "block"], mbert_filt.at[i, "type"]
        change = (prev_blk is not None and blk != prev_blk)
        if blk in ("fusion", "snp") and blk == prev_blk and typ != prev_type:
            change = True
        mbert_filt.at[i, "rule_above"] = bool(change)
        prev_blk, prev_type = blk, typ
    _write_csv(mbert_filt, "summary_t1e_modernbert_filtered")
    _write_hp_csv(mbert_filt, hp_map, "summary_t1e_modernbert_filtered_hp")
    render_main(mbert_filt, os.path.join(OUT_DIR, "summary_t1e_modernbert_filtered.png"), subtitle,
                title=TITLE_MBERT_FILTERED, show_n=True)
    render_main(mbert_filt, os.path.join(OUT_DIR, "summary_t1e_modernbert_filtered_no_n.png"), subtitle,
                title=TITLE_MBERT_FILTERED, show_n=False)
    render_hp(mbert_filt, hp_map, os.path.join(OUT_DIR, "summary_t1e_modernbert_filtered_hp.png"),
              title="T1e Late Fusion (best per CL × SNP) — Tuned Operating Points")
    write_t1e_latex(mbert_filt,
                    os.path.join(OUT_DIR, "summary_t1e_modernbert_filtered.tex"),
                    drop_perclass=True,
                    label="tab:t1e_late_fusion_mbert_filtered",
                    snp_desc=(r"SNP = PRS (meta-PRS EN/filtered, Kosteridis full/shared-AD-CV, "
                              r"Bellenguez/DeRojas/Kunkle PRS-CS, all-dedup PRS / PRS IVW) or "
                              r"FM (BMFM-SNP none$\cdot$MLP, NTv2 per-mod$\cdot$MLP), grouped by SNP "
                              r"model type (PRS block then FM block); within each block the best "
                              r"fusion method (by Test bACC) per CL $\times$ SNP combination is shown"),
                    caption_note=(" Best fusion method per (Clinical, SNP) combination chosen by Test "
                                  "bACC; full method-wise expansion in Table~\\ref{tab:t1e_late_fusion_expanded}."))

    # --- ModernBERT-B-ft SMALL: curated subset of 7 SNP variants (5 PRS + 2 FM),
    #     best fusion method per (CL, SNP), with full reference rows. ---
    snp_small = ["meta_prs", "meta_prs_filtered",
                 "kosteridis_prs", "prscs_kunkle", "prs_all_dedup_ivw",
                 "bmfm_none_mlp", "ntv2_permod_mlp"]
    small_full = assemble(snp_fusion_keys=snp_small, top_n=None, snp_ref_keys=snp_small,
                          group_by_type=True, clin_filter="modernbert_base_full_ft")
    fus_s = small_full[small_full.block == "fusion"].copy()
    fus_s_best = (fus_s.sort_values("bacc_m", ascending=False)
                       .drop_duplicates(["clin_variant", "snp_variant"], keep="first"))
    fus_s_best["_t"] = fus_s_best.snp_variant.map(_type_key)
    fus_s_best = fus_s_best.sort_values(["_t", "bacc_m"], ascending=[True, False],
                                          kind="stable").drop(columns="_t")
    small_df = pd.concat([fus_s_best,
                            small_full[small_full.block == "cl"],
                            small_full[small_full.block == "snp"]], ignore_index=True)
    small_df["rule_above"] = False
    prev_blk, prev_type = None, None
    for i in small_df.index:
        blk, typ = small_df.at[i, "block"], small_df.at[i, "type"]
        change = (prev_blk is not None and blk != prev_blk)
        if blk in ("fusion", "snp") and blk == prev_blk and typ != prev_type:
            change = True
        small_df.at[i, "rule_above"] = bool(change)
        prev_blk, prev_type = blk, typ
    _write_csv(small_df, "summary_t1e_modernbert_small")
    _write_hp_csv(small_df, hp_map, "summary_t1e_modernbert_small_hp")
    render_main(small_df, os.path.join(OUT_DIR, "summary_t1e_modernbert_small.png"), subtitle,
                title=TITLE_MBERT_SMALL, show_n=True)
    render_main(small_df, os.path.join(OUT_DIR, "summary_t1e_modernbert_small_no_n.png"), subtitle,
                title=TITLE_MBERT_SMALL, show_n=False)
    render_hp(small_df, hp_map, os.path.join(OUT_DIR, "summary_t1e_modernbert_small_hp.png"),
              title="T1e Late Fusion (small) — Tuned Operating Points")
    write_t1e_latex(small_df,
                    os.path.join(OUT_DIR, "summary_t1e_modernbert_small.tex"),
                    drop_perclass=True,
                    label="tab:t1e_late_fusion_mbert_small",
                    snp_desc=(r"SNP = PRS (meta-PRS EN/filtered, Kosteridis, Kunkle PRS-CS, "
                              r"all-dedup PRS IVW) or FM (BMFM-SNP none$\cdot$MLP, NTv2 per-mod$\cdot$MLP), "
                              r"grouped by SNP model type (PRS block then FM block); within each block "
                              r"the best fusion method (by Test bACC) per CL $\times$ SNP combination is shown"),
                    caption_note=(" Curated small subset of 7 SNP variants (5 PRS + 2 FM); best fusion "
                                  "method per (Clinical, SNP) combination chosen by Test bACC. Full "
                                  "method-wise expansion in Table~\\ref{tab:t1e_late_fusion_expanded}."))

    for name, df in (("summary_t1e_cleaned", main_df), ("summary_t1e_expanded", exp_df),
                     ("summary_t1e_modernbert", mbert_df)):
        print("=" * 92)
        print(f"  {name} ({len(df)} rows)  ·  {subtitle}")
        print("=" * 92)
        print(df[["Method", "SNP", "Aggregation", "bacc_m", "auc_m", "macro_f1_m", "n"]]
              .to_string(index=False))


if __name__ == "__main__":
    main()
