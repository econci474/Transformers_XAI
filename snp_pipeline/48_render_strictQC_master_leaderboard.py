"""
48_render_strictQC_master_leaderboard.py
========================================
Build the single wide cross-config + cross-task master leaderboard.

Layout (wide, landscape):
  rows = source x covar_mode
  per (LD config), grouped columns:
    n_snps | classification: val_bacc±std | cox: val_cindex±std | aao: val_r2±std
  → 4 configs × 4 cells = 16 metric columns + 2 ID columns.

Outputs at: outputs/strict_qc_prs/master_leaderboard.{tsv,png}.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

OUT_BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220/outputs/strict_qc_prs")
LD_CONFIGS = ["ld_1000kb_r2_0.8", "ld_500kb_r2_0.2",
               "ld_250kb_r2_0.1", "ld_50_5_r2_0.5"]

# When --beta-source prscs, the OUT_BASE is rerooted (see scripts 45/46/47).

# (task subdir, short label, kind, cols)
# kind="pm": cell = mean +/- std; cols = (mean_col, std_col)
# kind="ci": cell = mean [lo-hi]; cols = (mean_col, lo_col, hi_col)
def _build_tasks(split: str) -> list:
    """Return TASKS list with column names targeting the requested split.
    Naming convention matches scripts 45/46/47 emit shape:
      bacc / cindex / R2     -> {split}_<metric>_mean/_std
      OR (cls)               -> or_per_sd_{split}_<lo/hi>_mean
      HR (cox)               -> hr_per_sd_{split}_mean / hr_<lo/hi>_{split}_mean
    """
    return [
        ("classification", "bacc",   "pm", (f"{split}_bacc_mean",
                                              f"{split}_bacc_std")),
        ("classification", "OR",     "ci", (f"or_per_sd_{split}_mean",
                                              f"or_per_sd_{split}_lo_mean",
                                              f"or_per_sd_{split}_hi_mean")),
        ("cox",            "cindex", "pm", (f"{split}_cindex_mean",
                                              f"{split}_cindex_std")),
        ("cox",            "HR",     "ci", (f"hr_per_sd_{split}_mean",
                                              f"hr_lo_{split}_mean",
                                              f"hr_hi_{split}_mean")),
        ("aao_regression", "R2",     "pm", (f"{split}_r2_mean",
                                              f"{split}_r2_std")),
    ]


def _fmt_pm(m, s):
    if pd.isna(m): return ""
    if pd.isna(s): return f"{m:.3f}"
    return f"{m:.3f}+/-{s:.3f}"


def _fmt_ci(m, lo, hi):
    if pd.isna(m): return ""
    if pd.isna(lo) or pd.isna(hi): return f"{m:.2f}"
    return f"{m:.2f} [{lo:.2f}-{hi:.2f}]"


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"])
    ap.add_argument("--split", default="val", choices=["val","test"],
                     help="which fold's metrics to pull (val = default; test = "
                          "held-out leaderboard with test-fold OR/HR refits).")
    args = ap.parse_args()
    global OUT_BASE
    if args.beta_source != "raw":
        OUT_BASE = Path(str(OUT_BASE) + f"_{args.beta_source}")
    print(f"reading from {OUT_BASE}  (split={args.split})")

    TASKS = _build_tasks(args.split)

    # PRS-CS uses the unpruned 115-SNP pool — no LD config loop. Single dir.
    cfg_iter = [None] if args.beta_source == "prscs" else LD_CONFIGS

    frames = []
    for cfg in cfg_iter:
        for task, short, kind, cols in TASKS:
            p = (OUT_BASE / task / "leaderboard.csv") if cfg is None \
                else (OUT_BASE / cfg / task / "leaderboard.csv")
            if not p.exists():
                print(f"[skip] {p}"); continue
            df = pd.read_csv(p)
            # Tolerate missing val-fold CI cols on older leaderboards
            for c in cols:
                if c not in df.columns: df[c] = float("nan")
            df["key"] = df["source"] + "|" + df["covar_mode"]
            if kind == "pm":
                mcol, scol = cols
                df["cell"] = [_fmt_pm(m, s) for m, s in zip(df[mcol], df[scol])]
            else:  # ci
                mcol, lcol, hcol = cols
                df["cell"] = [_fmt_ci(m, lo, hi)
                                for m, lo, hi in zip(df[mcol], df[lcol], df[hcol])]
            tag = cfg if cfg is not None else "strictQC"  # single column tag for prscs
            sub = df[["key", "n_snps", "cell"]].rename(columns={
                "n_snps": f"n_snps__{tag}",
                "cell":   f"{short}__{tag}",
            })
            sub_idx = sub.set_index("key")
            frames.append((tag, short, sub_idx))

    # Outer-join everything
    out = None
    for tag, short, sub in frames:
        cols_to_use = [f"{short}__{tag}"]
        ns_col = f"n_snps__{tag}"
        if ns_col in sub.columns and (out is None or ns_col not in out.columns):
            cols_to_use = [ns_col] + cols_to_use
        joined = sub[cols_to_use]
        out = joined if out is None else out.join(joined, how="outer")

    # Reset key → (source, covar_mode)
    out = out.reset_index()
    out[["source", "covar_mode"]] = out["key"].str.split("|", expand=True)
    out = out.drop(columns=["key"])

    # Column order: source, covar_mode, then per-tag block:
    # n_snps | bacc | OR | cindex | HR | R2
    tag_iter = ["strictQC"] if args.beta_source == "prscs" else LD_CONFIGS
    cols_ordered = ["source", "covar_mode"]
    for tag in tag_iter:
        cols_ordered.append(f"n_snps__{tag}")
        cols_ordered.append(f"bacc__{tag}")
        cols_ordered.append(f"OR__{tag}")
        cols_ordered.append(f"cindex__{tag}")
        cols_ordered.append(f"HR__{tag}")
        cols_ordered.append(f"R2__{tag}")
    cols_ordered = [c for c in cols_ordered if c in out.columns]
    out = out[cols_ordered]

    # Pretty column rename for display (TSV keeps the explicit names).
    # Three-line headers so each column can stay narrow: <window>\n<r²>\n<metric>
    # Shorter labels fit inside narrow metric columns without bleeding.
    SHORT_CFG = {
        "ld_1000kb_r2_0.8": ("1Mb",   "r²<.8"),
        "ld_500kb_r2_0.2":  ("0.5Mb", "r²<.2"),
        "ld_250kb_r2_0.1":  ("0.25Mb","r²<.1"),
        "ld_50_5_r2_0.5":   ("50/5",  "r²<.5"),
        "strictQC":         ("strict-QC", "no LD prune"),
    }
    def _short(c: str) -> str:
        for cfg, (window, r2) in SHORT_CFG.items():
            if c.endswith(f"__{cfg}"):
                metric = c.split("__")[0]
                return f"{window}\n{r2}\n{metric}"
        return c

    # Drop Kosteridis_novel_AD from raw-β display (strict subset of MTAG_AD)
    if args.beta_source == "raw":
        out = out[out["source"] != "Kosteridis_novel_AD"].copy()

    # Sort: covar_mode block; then pin covariates_only at the top of each block
    # (clinical baseline reference), then sort the rest by BalAcc descending.
    sort_col = "bacc__strictQC" if args.beta_source == "prscs" else "bacc__ld_1000kb_r2_0.8"
    if sort_col in out.columns:
        def _bacc_val(v):
            if not isinstance(v, str) or v == "": return -1.0
            try:    return float(v.split("+/-")[0])
            except: return -1.0
        out["_sort"] = out[sort_col].apply(_bacc_val)
        # Lower _baseline value = higher in block. covariates_only=0, else=1.
        out["_baseline"] = (out["source"] != "covariates_only").astype(int)
        out = (out.sort_values(["covar_mode","_baseline","_sort"],
                                 ascending=[True, True, False])
                  .drop(columns=["_sort","_baseline"]))

    tsv_p = OUT_BASE / f"master_leaderboard_{args.split}.tsv"
    out.to_csv(tsv_p, sep="\t", index=False)
    print(f"wrote master TSV: {tsv_p}")

    # PNG — wide landscape, multi-line headers, white header background,
    # top-per-column bolding only.
    n_rows = len(out)
    fig_w = 22 if args.beta_source == "prscs" else 50
    # Add extra vertical space for the (taller) header + 2-line title.
    fig_h = max(2.5, 0.32 * n_rows + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h)); ax.axis("off")
    title = ("Strict-QC PRS — master leaderboard across LD configs and tasks\n"
              f"{args.split.upper()} set; "
              f"bacc / cindex / R² = TRAIN-fitted model evaluated on {args.split} "
              f"(mean+/-std over 3 seeds, held-out generalisation);   "
              f"OR / HR = {args.split.upper()}-fold REFIT "
              f"(coefficient + 95% Wald CI on {args.split} itself)\n"
              "bold = top per column (rank metrics) or CI excludes 1 (OR/HR)")
    if args.beta_source == "prscs":
        title += ("\nn_snps = | sumstats ∩ HM3 LD ref ∩ strict-QC bim | "
                   "(PRS-CS posterior coverage)")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=16, loc="left")
    show = out.copy()
    for c in show.columns:
        if c.startswith("n_snps__"):
            show[c] = show[c].apply(lambda v: "" if pd.isna(v) else str(int(v)))

    display_headers = [_short(c) if "__" in c else c for c in show.columns]

    # Narrower metric columns (multi-line headers let columns stay tight).
    # PRS-CS layout has fewer columns and longer header labels ("strict-QC /
    # no LD prune / <metric>"), so per-column widths are looser.
    is_prscs = args.beta_source == "prscs"
    widths = []
    for c in show.columns:
        if c == "source":              widths.append(0.075)
        elif c == "covar_mode":        widths.append(0.060)
        elif c.startswith("n_snps"):
            widths.append(0.045 if is_prscs else 0.022)
        elif c.startswith("OR__") or c.startswith("HR__"):
            widths.append(0.075 if is_prscs else 0.055)
        else:
            widths.append(0.060 if is_prscs else 0.034)
    s = sum(widths); widths = [w / s for w in widths]

    tbl = ax.table(cellText=show.values.tolist(), colLabels=display_headers,
                    loc="center", cellLoc="center", colWidths=widths)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.45)
    # Header band: white bg, black bold text. PRS-CS has fewer columns so the
    # header text wraps over more visual rows — give it more vertical room.
    HDR_HEIGHT = 0.110 if args.beta_source == "prscs" else 0.075
    for j in range(len(show.columns)):
        c = tbl[(0, j)]
        c.set_facecolor("white")
        c.set_text_props(color="black", weight="bold", fontsize=11)
        c.set_height(HDR_HEIGHT)
        c.set_edgecolor("black")
    # Bold the top value per metric column WITHIN EACH covar_mode block.
    # So each (column × covar_mode) cell-group has one bold winner.
    # Higher-is-better: bacc, cindex, R2, OR, HR. Skip n_snps + ID cols.
    col_idx = {c: i for i, c in enumerate(show.columns)}
    def _mean(s: str) -> float:
        if not isinstance(s, str) or s == "": return float("-inf")
        # Format options: "0.730+/-0.090" (pm) or "1.23 [0.94-1.61]" (ci)
        try:
            if "+/-" in s:           return float(s.split("+/-")[0])
            return float(s.split(" ")[0])
        except (ValueError, IndexError):
            return float("-inf")
    def _ci_excludes_one(s: str) -> bool:
        # Only meaningful for OR/HR ci-format cells.
        if not isinstance(s, str) or "[" not in s: return False
        try:
            inside = s.split("[")[1].rstrip("]")
            lo, hi = (float(x) for x in inside.split("-"))
            return lo > 1.0 or hi < 1.0
        except (ValueError, IndexError):
            return False
    out_idx = out.reset_index(drop=True)
    # Top-per-column bolding applies ONLY to the rank metrics (bacc/cindex/R2).
    # OR/HR cells use a separate, stricter rule below: bold iff 95% CI excludes
    # 1 (a non-significant OR/HR shouldn't be bolded just for being top-of-block).
    rank_prefixes = ("bacc__", "cindex__", "R2__")
    # Display-precision-aware bolding: every row whose mean rounds to the same
    # displayed value as the column max gets bolded (so e.g. 0.7296 / 0.7301 /
    # 0.7297 — all rendered as 0.730 — are bolded together, not just argmax).
    _display_dp = {"bacc__": 3, "cindex__": 3, "R2__": 3}
    for col in show.columns:
        prefix = next((p for p in rank_prefixes if col.startswith(p)), None)
        if prefix is None: continue
        dp = _display_dp[prefix]
        for mode in out_idx["covar_mode"].unique():
            mask = out_idx["covar_mode"] == mode
            sub = out_idx.loc[mask, col]
            if sub.empty: continue
            means = [_mean(v) for v in sub.tolist()]
            if all(m == float("-inf") for m in means): continue
            top_rounded = round(max(means), dp)
            for li, m in enumerate(means):
                if m == float("-inf"): continue
                if round(m, dp) == top_rounded:
                    global_row = int(sub.index[li])
                    tbl[(global_row + 1, col_idx[col])].set_text_props(
                        weight="bold")
    # OR/HR cells: bold iff the 95% CI excludes 1 (statistically significant).
    # Never bold just for being top-of-column — a top OR with CI crossing 1
    # is non-significant and should NOT be highlighted.
    for col in show.columns:
        if not (col.startswith("OR__") or col.startswith("HR__")):
            continue
        for ri, v in enumerate(out_idx[col].tolist()):
            if _ci_excludes_one(v):
                tbl[(ri + 1, col_idx[col])].set_text_props(weight="bold")

    png_p = OUT_BASE / f"master_leaderboard_{args.split}.png"
    fig.savefig(png_p, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote master PNG: {png_p}")
    print("\nMaster leaderboard preview (first 12 rows):")
    print(out.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
