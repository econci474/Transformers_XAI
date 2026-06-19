r"""
30c_update_master_with_corrected_fm.py
=======================================
Replace the 7 FM rows of the master leaderboard comparison TSV/PNG with the
corrected-GTEx FM rows from the 21 predictions.tsv files that 30v3 saved
when re-run with --save-predictions on Colab L4 GPU.

Workflow:
  1. Read all 21 predictions.tsv from
     D:/ADNI_SNP_Omni2.5M_20140220/AlphaGenome_gtex_corrected/fm_top_predictions_top_7/
  2. Build a spec list for the 7 FM configs (matching what was ACTUALLY run —
     row 7 BMFM-MLP attn_bias_per_modality is d=0.1 not d=0.3; row 3 BMFM XGB
     is n_estimators=20 not 200, per user 2026-06-08).
  3. For each (spec, split, seed), use _compute_train_fitted_metrics from
     _prs_quantile_plot.py to compute bAcc/AUC/R²_lia_null/R²_lia_cov.
  4. Aggregate across seeds → mean ± std.
  5. Load the original master TSV; drop the FM rows; append the 7 new rows.
  6. Write outputs/corrected_gtex/master_leaderboard_comparison_corrected_gtex
     .{tsv,png}.

Caveats:
  - The original FM rows in master_leaderboard_comparison_train.tsv were
    computed by _prs_quantile_plot.py from an earlier set of predictions
    extracted on Colab A100 with the buggy AG features. The new rows here
    are from the strict-filter corrected-GTEx 25-D layout.
  - The non-FM rows (baseline / LD-C+T / EN / PRS-CS) are NOT recomputed —
    they're copied verbatim from the original TSV.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import numpy as np

HERE = Path(__file__).parent

_FAMILY_ORDER = ("baseline", "LD-C+T", "EN", "PRS-CS", "FM")

# PRS+cov base-names to drop from the filtered leaderboard (2026-06-10 user
# request: declutter the filtered view; the unfiltered master TSV keeps them).
_EXCLUDE_PRS_COV_BASE = {
    "prs_all_dedup",
    "prs_all_dedup_filtered",
    "Kosteridis_MTAG_AD",
    "prs_all_dedup_EN_dosage",
}

# Strip trailing hyperparameter block ` (...)` from FM display names
# (2026-06-10): the master TSV name "X · Y · Z · MLP (w=256, d=0.1, lr=0.001)"
# shows in the filtered view as just "X · Y · Z · MLP".
import re as _re
_FM_HP_TAIL = _re.compile(r"\s*\([^)]*\)$")

# Hot-patch _prs_quantile_plot.py's PRED_DIR before importing it
NEW_PRED_DIR = Path("D:/ADNI_SNP_Omni2.5M_20140220/AlphaGenome_gtex_corrected/"
                     "fm_top_predictions_top_7/diff_attn_v3")
ORIG_MASTER_TSV = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/prs_quantile_plot/"
                        "master_leaderboard_comparison_train.tsv")
OUT_DIR = HERE.parent / "snp_pipeline/outputs/corrected_gtex"


def _import_qp():
    spec = importlib.util.spec_from_file_location("qp", HERE / "_prs_quantile_plot.py")
    qp = importlib.util.module_from_spec(spec)
    # Override PRED_DIR before module-level code runs
    # (PRED_DIR is a module-level constant; we override after exec)
    spec.loader.exec_module(qp)
    qp.PRED_DIR = NEW_PRED_DIR
    return qp


# 7 FM specs — names match what we ACTUALLY ran (not the original d=0.3 etc.)
FM_CONFIGS = [
    # (family_label_for_name, model, head, aggregation, mode, seq_length, hp_dict)
    ("BMFM-SNP", "bmfm_snp", "mlp2", "chrom_hier", "attn_bias_single", "101bp",
     dict(mlp_width=256, mlp_dropout=0.1, lr=0.001)),
    ("BMFM-SNP", "bmfm_snp", "mlp2", "chrom_hier", "none", "101bp",
     dict(mlp_width=256, mlp_dropout=0.1, lr=0.001)),
    ("BMFM-SNP", "bmfm_snp", "xgb", "global_attn", "none", "1001bp",
     dict(xgb_lr=0.05, xgb_max_depth=4, xgb_n_estimators=20)),
    ("NTv2", "ntv2", "mlp2", "chrom_hier", "attn_bias_per_modality", "101bp",
     dict(mlp_width=512, mlp_dropout=0.3, lr=0.0003)),
    ("NTv2", "ntv2", "mlp2", "chrom_hier", "attn_bias_single", "101bp",
     dict(mlp_width=512, mlp_dropout=0.3, lr=0.0003)),
    ("NTv2", "ntv2", "mlp2", "chrom_hier", "none", "101bp",
     dict(mlp_width=512, mlp_dropout=0.3, lr=0.0003)),
    ("BMFM-SNP", "bmfm_snp", "mlp2", "chrom_hier", "attn_bias_per_modality", "101bp",
     dict(mlp_width=256, mlp_dropout=0.1, lr=0.001)),
    # 2026-06-09: 2 NTv2 MLP combo modes added per user request — appeared
    # near the top of the corrected-GTEx GPU sweep table.
    ("NTv2", "ntv2", "mlp2", "chrom_hier", "value_mult", "101bp",
     dict(mlp_width=512, mlp_dropout=0.3, lr=0.0003)),
    ("NTv2", "ntv2", "mlp2", "chrom_hier", "attn_bias_per_modality_plus_value_mult", "101bp",
     dict(mlp_width=512, mlp_dropout=0.3, lr=0.0003)),
]


def _spec_name(backbone, mode, agg, head, hp):
    """Format the row 'name' column to match the master TSV's pattern."""
    head_label = "MLP" if head.startswith("mlp") else "XGB"
    if head_label == "MLP":
        hp_str = f"w={hp['mlp_width']}, d={hp['mlp_dropout']}, lr={hp['lr']}"
    else:
        hp_str = (f"xgb_lr={hp['xgb_lr']}, max_depth={hp['xgb_max_depth']}, "
                  f"n_estimators={hp['xgb_n_estimators']}")
    return f"{backbone} · {mode} · {agg} · {head_label} ({hp_str})"


def _build_fm_spec(qp, cfg_tuple) -> dict:
    """Build a spec dict that _compute_train_fitted_metrics will recognise as FM."""
    backbone, model, head, agg, mode, seq_length, hp = cfg_tuple
    name = _spec_name(backbone, mode, agg, head, hp)
    # Need the canonical compute_fn that returns DataFrame[Patient_ID, prs]
    # for the requested split × seed. Reuse _fm_compute_fn but with our config
    # filter dict matching the predictions.tsv columns.
    cfg_filter = {
        "model": model, "func_integration_mode": mode, "head": head,
        "seq_length": seq_length, "aggregation": agg,
    }
    if head.startswith("mlp"):
        cfg_filter.update({
            "mlp_width": hp["mlp_width"],
            "mlp_dropout": hp["mlp_dropout"],
            "lr": hp["lr"],
        })
    return {
        "family": "FM",
        "name": name,
        "snp_count": 128,
        "compute_fn": qp._fm_compute_fn(cfg_filter, qp._load_fm_predictions_df()),
        "_filter": cfg_filter,
    }


def _parse_mean(s: str) -> float:
    """Pull the mean out of a 'm+/-s' aggregated string. Returns -inf on miss."""
    if not isinstance(s, str) or not s:
        return float("-inf")
    try:
        return float(s.split("+/-")[0])
    except ValueError:
        return float("-inf")


def _sort_within_family_by_balacc(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by _FAMILY_ORDER (stable) then val_bacc descending within each family."""
    fam_rank = {f: i for i, f in enumerate(_FAMILY_ORDER)}
    out = df.copy()
    out["_fam_rank"] = out["family"].map(fam_rank).fillna(len(_FAMILY_ORDER))
    out["_val_bacc_mean"] = out["val_bacc"].map(_parse_mean)
    out = out.sort_values(["_fam_rank", "_val_bacc_mean"],
                            ascending=[True, False],
                            kind="mergesort")
    return out.drop(columns=["_fam_rank", "_val_bacc_mean"]).reset_index(drop=True)


def _fmt_2dp(s: str) -> str:
    """Convert '0.xxx+/-0.xxx' (any precision) → '0.xx ± 0.xx'. Pass empty/garbage through."""
    if not isinstance(s, str) or "+/-" not in s:
        return s if isinstance(s, str) else ""
    try:
        m, sd = s.split("+/-")
        return f"{float(m):.2f} ± {float(sd):.2f}"
    except ValueError:
        return s


def _render_filtered_leaderboard(df: pd.DataFrame, out_dir: Path, K_p: float,
                                    basename: str = "master_leaderboard_filtered") -> None:
    """Filtered, 2dp variant of the master leaderboard. Emits PNG + PDF + TeX.

    Rows kept:
      - covariates-only baseline (family == 'baseline')
      - all PRS+covariates rows (name contains ' + covariates')
      - all FM rows (family == 'FM', raw — FM has no cov-fusion variant)

    Columns kept (val only): family, name, n_snps, val_bacc, val_auc, val_R2_lia_cov.
    Numerics shown at 2 decimal places with unicode '±'.
    Sort: within _FAMILY_ORDER then val_bacc descending."""
    import matplotlib.pyplot as plt

    is_baseline = df["family"] == "baseline"
    is_prs_cov  = df["name"].fillna("").str.contains(" + covariates", regex=False)
    is_fm       = df["family"] == "FM"
    sub = df[is_baseline | is_prs_cov | is_fm].copy()

    # Drop excluded PRS+cov rows (match on base name before " + covariates").
    def _base_name(n):
        return n.split(" + covariates")[0].strip() if isinstance(n, str) else ""
    excluded = sub["name"].map(_base_name).isin(_EXCLUDE_PRS_COV_BASE)
    if excluded.any():
        print(f"  [filter] excluding {int(excluded.sum())} PRS+cov rows: "
              f"{sorted(set(sub.loc[excluded, 'name'].map(_base_name)))}")
        sub = sub.loc[~excluded].copy()

    # Strip trailing hyperparameter block from FM display names.
    fm_mask = sub["family"] == "FM"
    if fm_mask.any():
        sub.loc[fm_mask, "name"] = (sub.loc[fm_mask, "name"]
                                       .str.replace(_FM_HP_TAIL, "", regex=True))

    sub = _sort_within_family_by_balacc(sub)

    keep_cols = ["family", "name", "n_snps", "val_bacc", "val_auc", "val_R2_lia_cov"]
    sub = sub[keep_cols].copy()

    # Fix stale n_snps for meta-PRS rows. The master TSV stamps the
    # number of source-PRS *columns* (17 / 7), not the SNP count. Correct
    # methodology: meta-PRS doesn't apply LD pruning at the meta level (the
    # EN combines per-source PRS values, not SNPs), so n_snps should be the
    # union of unique lead-SNP rsIDs across constituent sources:
    #   meta_prs_EN_combined: 119 unique rsIDs across 17 studies
    #   meta_prs_EN_filtered: 102 unique rsIDs across 7 META_FILTERED_SOURCES
    # (apples-to-apples with the FM rows' 128 unique-SNP convention).
    _META_PRS_NSNPS = {
        "meta_prs_EN_combined": 119,
        "meta_prs_EN_filtered": 102,
    }
    for src_name, n in _META_PRS_NSNPS.items():
        mask = sub["name"].fillna("").str.startswith(src_name)
        sub.loc[mask, "n_snps"] = n
    for c in ("val_bacc", "val_auc", "val_R2_lia_cov"):
        sub[c] = sub[c].map(_fmt_2dp)

    # ── Bolding masks (shared between LaTeX and PNG) ──────────────────────
    # Top per column for BalAcc/AUC; top + any positive for R²_cov.
    def _mean_2dp(s: str) -> float:
        if not isinstance(s, str) or "±" not in s:
            return float("-inf")
        try:
            return float(s.split("±")[0].strip())
        except ValueError:
            return float("-inf")
    bold_rules = [("val_bacc", False), ("val_auc", False),
                  ("val_R2_lia_cov", True)]
    bold_mask = {col: [False] * len(sub) for col, _ in bold_rules}
    for col, bold_positive in bold_rules:
        vals = [_mean_2dp(v) for v in sub[col].tolist()]
        if all(v == float("-inf") for v in vals):
            continue
        top = round(max(vals), 2)
        for i, v in enumerate(vals):
            if v == float("-inf"):
                continue
            if round(v, 2) == top or (bold_positive and v > 0.0):
                bold_mask[col][i] = True

    # ── TeX ───────────────────────────────────────────────────────────────
    # Use sentinels for math-header AND inline bolding so escape=True passes
    # them through unchanged; swap for real LaTeX after rendering.
    _R2_SENTINEL = "R2COVHEADER"
    _BOLD_OPEN = "BOLDSTART"
    _BOLD_CLOSE = "BOLDEND"
    tex_sub = sub.copy()
    for col, _ in bold_rules:
        for i, is_bold in enumerate(bold_mask[col]):
            if is_bold:
                tex_sub.iat[i, tex_sub.columns.get_loc(col)] = (
                    f"{_BOLD_OPEN}{sub[col].iloc[i]}{_BOLD_CLOSE}")
    tex_df = tex_sub.rename(columns={
        "family":          "Family",
        "name":            "Method",
        "n_snps":          "n SNPs",
        "val_bacc":        "BalAcc",
        "val_auc":         "AUC",
        "val_R2_lia_cov":  _R2_SENTINEL,
    })
    tex_body = tex_df.to_latex(index=False, escape=True, na_rep="",
                                column_format="llrrrr")
    tex_body = tex_body.replace(_R2_SENTINEL, r"$R^2_{\mathrm{liab}|\mathrm{cov}}$")
    tex_body = tex_body.replace(_BOLD_OPEN, r"\textbf{").replace(_BOLD_CLOSE, "}")
    # Insert a multicolumn title row at the top of the tabular (above the
    # column headers, below \toprule).
    n_cols = len(tex_df.columns)
    title_row = (f"\\multicolumn{{{n_cols}}}{{c}}{{\\textbf{{T1e: Genetic "
                  "model performance on Validation set}} \\\\\n\\midrule\n")
    tex_body = tex_body.replace("\\toprule\n", "\\toprule\n" + title_row, 1)
    caption = (
        "Top performing PRS and DNA foundation models for T1e task on "
        "validation fold. Mean $\\pm$ std over 3 seeds. "
        f"$R^2$ on the liability scale via Lee 2012 ($K_p={K_p:.3f}$), "
        "incremental over baseline covariates (age, sex, APOE4, APOE2); "
        "the covariates-only row is $0$ by construction. "
        "Sorted within family by BalAcc descending; bold marks the top "
        "value per column (and any positive value in the $R^2$ column). "
        "For meta-PRS rows, n SNPs is the union of unique lead-SNP rsIDs "
        "across constituent sources (combined: 119 unique SNPs over 17 "
        "studies; filtered: 102 unique SNPs over 7 studies). "
        "Family abbreviations: LD-C+T = LD clumping with p-value thresholding "
        "(raw per-source effect sizes); EN = Elastic-Net meta-PRS (logistic "
        "regression over per-source PRS values); PRS-CS = PRS Continuous "
        "Shrinkage (Bayesian posterior effect sizes; Ge et al. 2019); "
        "FM = Foundation-Model embeddings (BMFM-SNP or NTv2 DNA language "
        "models with GTEx AlphaGenome functional features)."
    )
    # Caption + label at the bottom (after the tabular).
    # Wrap tabular in resizebox so it fits page width.
    tex_doc = (
        "% Auto-generated by 30c_update_master_with_corrected_fm.py.\n"
        "% Requires the booktabs and graphicx packages.\n"
        "\\begin{table}[ht]\n"
        "\\centering\n"
        "\\resizebox{\\textwidth}{!}{%\n"
        f"{tex_body}"
        "}\n"
        f"\\caption{{{caption}}}\n"
        "\\label{tab:prs_master_filtered}\n"
        "\\end{table}\n"
    )
    tex_p = out_dir / f"{basename}.tex"
    tex_p.write_text(tex_doc, encoding="utf-8")
    print(f"  wrote {tex_p}")

    # ── PNG + PDF (one Figure, savefig twice) ─────────────────────────────
    n_rows = len(sub)
    fig_h = max(3.0, 0.35 * n_rows + 2.0)
    fig_w = 18  # wide enough for the longest FM method label
    fig, ax = plt.subplots(figsize=(fig_w, fig_h)); ax.axis("off")
    title = (
        "PRS master comparison — filtered (PRS+cov + FM + covariates-only reference)\n"
        f"VAL fold, mean ± std over 3 seeds; R²_liab over cov "
        f"(incremental over age/sex/APOE4/APOE2, Lee 2012, K_p={K_p:.3f}).\n"
        "Sorted within family by BalAcc descending; bold = top per column (R²: any positive)."
    )
    ax.set_title(title, fontsize=11, fontweight="bold", pad=14, loc="left")

    display_cols = ["Family", "Method", "n SNPs", "BalAcc", "AUC", "R²_liab\nover cov"]
    widths = [0.06, 0.58, 0.05, 0.105, 0.105, 0.09]
    s = sum(widths); widths = [w / s for w in widths]
    tbl = ax.table(cellText=sub.values.tolist(), colLabels=display_cols,
                    loc="center", cellLoc="center", colWidths=widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.35)
    HDR_HEIGHT = 0.085
    for j in range(len(display_cols)):
        c = tbl[(0, j)]
        c.set_facecolor("white")
        c.set_text_props(color="black", weight="bold", fontsize=10)
        c.set_height(HDR_HEIGHT); c.set_edgecolor("black")

    # Bold rules: top per column for BalAcc/AUC; top + any positive for R²_cov.
    def _mean_2dp(s: str) -> float:
        if not isinstance(s, str) or "±" not in s:
            return float("-inf")
        try:
            return float(s.split("±")[0].strip())
        except ValueError:
            return float("-inf")
    col_idx = {c: i for i, c in enumerate(sub.columns)}
    metric_cols = [("val_bacc", False), ("val_auc", False), ("val_R2_lia_cov", True)]
    for col, bold_positive in metric_cols:
        vals = [_mean_2dp(v) for v in sub[col].tolist()]
        if all(v == float("-inf") for v in vals): continue
        top = round(max(vals), 2)
        for i, v in enumerate(vals):
            if v == float("-inf"): continue
            is_top = round(v, 2) == top
            is_positive = bold_positive and v > 0.0
            if is_top or is_positive:
                tbl[(i + 1, col_idx[col])].set_text_props(weight="bold")

    # Family-group separators: thin dashed line between adjacent families.
    fam_seq = sub["family"].tolist()
    for i in range(len(fam_seq) - 1):
        if fam_seq[i] != fam_seq[i + 1]:
            for j in range(len(sub.columns)):
                tbl[(i + 1, j)].visible_edges = "BT"

    png_p = out_dir / f"{basename}.png"
    pdf_p = out_dir / f"{basename}.pdf"
    fig.savefig(png_p, dpi=200, bbox_inches="tight")
    fig.savefig(pdf_p, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png_p}")
    print(f"  wrote {pdf_p}")


def main():
    qp = _import_qp()
    # Manually re-load predictions df with our PRED_DIR override
    fm_df = qp._load_fm_predictions_df()
    print(f"[load] {len(fm_df)} prediction rows from {NEW_PRED_DIR}")
    if fm_df.empty:
        sys.exit(f"[ERROR] no predictions.tsv files under {NEW_PRED_DIR}")
    # Print distinct (model, mode, head, seq_length, aggregation, mlp_width,
    # mlp_dropout, lr) combos we see in the predictions
    cfg_cols = [c for c in ["model","func_integration_mode","head","seq_length",
                              "aggregation","mlp_width","mlp_dropout","lr"]
                if c in fm_df.columns]
    distinct = fm_df[cfg_cols].drop_duplicates()
    print(f"[load] {len(distinct)} distinct configs in predictions:")
    for _, r in distinct.iterrows():
        print(f"  {dict(r)}")

    # Load subject covariates for R²_lia_cov computation
    cov_spec = importlib.util.spec_from_file_location(
        "lib", HERE / "_prs_strict_qc_lib.py")
    lib = importlib.util.module_from_spec(cov_spec)
    cov_spec.loader.exec_module(lib)
    cov = lib.load_subject_covariates()
    print(f"[cov] loaded {len(cov)} subject covariate rows")

    K_p = qp.DEFAULT_K_POP
    print(f"[lee] K_p (population AD prevalence) = {K_p:.3f}")

    # Build specs and compute per-(spec, split, seed) metrics
    rows = []
    for cfg_tuple in FM_CONFIGS:
        spec = _build_fm_spec(qp, cfg_tuple)
        print(f"\n[spec] {spec['name']}")
        per_split = {}
        for split in ("val", "test"):
            seed_runs = []
            for seed in (0, 1, 2):
                m = qp._compute_train_fitted_metrics(spec, split, seed, K_p, cov)
                if m is not None:
                    seed_runs.append(m)
                    print(f"  {split} seed={seed}: bAcc={m['bacc']:.3f}  "
                          f"AUC={m['auc']:.3f}  "
                          f"R²_null={m['R2_lia_null']:+.4f}  "
                          f"R²_cov={m['R2_lia_cov']:+.4f}  "
                          f"N={m['N']}")
                else:
                    print(f"  {split} seed={seed}: [empty]")
            per_split[split] = qp._aggregate_comparison_metrics(seed_runs)
        if not per_split["val"] and not per_split["test"]:
            print("  [skip] no rows")
            continue
        row = {"family": "FM", "name": spec["name"], "n_snps": 128}
        for sp in ("val", "test"):
            agg = per_split[sp]
            if not agg:
                for k in ("N","cases","controls","bacc","auc",
                          "R2_lia_null","R2_lia_cov"):
                    row[f"{sp}_{k}"] = ""
                continue
            row[f"{sp}_N"]        = f"{agg['N_mean']:.0f}"
            row[f"{sp}_cases"]    = f"{agg['n_pos_mean']:.0f}"
            row[f"{sp}_controls"] = f"{agg['n_neg_mean']:.0f}"
            row[f"{sp}_bacc"]     = f"{agg['bacc_mean']:.3f}+/-{agg['bacc_std']:.3f}"
            row[f"{sp}_auc"]      = f"{agg['auc_mean']:.3f}+/-{agg['auc_std']:.3f}"
            row[f"{sp}_R2_lia_null"] = (f"{agg['R2_lia_null_mean']:.4f}+/-"
                                        f"{agg['R2_lia_null_std']:.4f}")
            row[f"{sp}_R2_lia_cov"]  = (f"{agg['R2_lia_cov_mean']:.4f}+/-"
                                        f"{agg['R2_lia_cov_std']:.4f}")
        rows.append(row)

    new_fm_df = pd.DataFrame(rows)
    print(f"\n[fm-new] built {len(new_fm_df)} FM rows")

    # Load original master TSV, drop FM rows, append new ones
    orig = pd.read_csv(ORIG_MASTER_TSV, sep="\t", dtype=str)
    non_fm = orig[orig["family"] != "FM"].copy()
    print(f"[orig]  loaded {len(orig)} rows; {len(non_fm)} non-FM kept verbatim")
    out = pd.concat([non_fm, new_fm_df], ignore_index=True)

    # Sort within each family by val_bacc descending (2026-06-09).
    out = _sort_within_family_by_balacc(out)
    print(f"[sort]  rows reordered within family by val_bacc")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_tsv = OUT_DIR / "master_leaderboard_comparison_corrected_gtex.tsv"
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"[out]   {out_tsv}")

    # Render PNG using the existing renderer
    title = (
        "PRS master comparison — TRAIN-FITTED held-out predictive accuracy "
        "(val + test)\n"
        "FM rows: corrected-GTEx AlphaGenome features (strict allow-list, "
        "25-D, 6-modality summaries) — rerun on Colab L4 GPU 2026-06-08.\n"
        "Non-FM rows verbatim from D:/ADNI_SNP_Omni2.5M_20140220/outputs/"
        "prs_quantile_plot/master_leaderboard_comparison_train.tsv.\n"
        f"mean +/- std over 3 seeds; R²_liab via Lee 2012 (K_p={K_p:.3f}); "
        "bold = top per column."
    )
    qp._render_master_comparison_png(
        out, OUT_DIR, K_p,
        output_filename="master_leaderboard_comparison_corrected_gtex.png",
        title_override=title)
    print(f"[out]   {OUT_DIR / 'master_leaderboard_comparison_corrected_gtex.png'}")

    # 2026-06-09: filtered variant — covariates-only + PRS+cov + FM rows,
    # 2dp, PNG + PDF + LaTeX.
    print(f"\n[filt]  rendering filtered leaderboard (PNG + PDF + TeX)")
    _render_filtered_leaderboard(out, OUT_DIR, K_p,
                                   basename="master_leaderboard_filtered")


if __name__ == "__main__":
    main()
