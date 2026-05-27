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

# (task subdir, primary metric mean+std column pair, short label for display)
TASKS = [
    ("classification", ("val_bacc_mean", "val_bacc_std"), "bacc"),
    ("cox",            ("val_cindex_mean", "val_cindex_std"), "cindex"),
    ("aao_regression", ("val_r2_mean", "val_r2_std"),         "R2"),
]


def _fmt_pm(m, s):
    if pd.isna(m): return ""
    if pd.isna(s): return f"{m:.3f}"
    return f"{m:.3f}+/-{s:.3f}"


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-source", default="raw", choices=["raw","prscs"])
    args = ap.parse_args()
    global OUT_BASE
    if args.beta_source != "raw":
        OUT_BASE = Path(str(OUT_BASE) + f"_{args.beta_source}")
    print(f"reading from {OUT_BASE}")

    # PRS-CS uses the unpruned 115-SNP pool — no LD config loop. Single dir.
    cfg_iter = [None] if args.beta_source == "prscs" else LD_CONFIGS

    frames = []
    for cfg in cfg_iter:
        for task, (mcol, scol), short in TASKS:
            p = (OUT_BASE / task / "leaderboard.csv") if cfg is None \
                else (OUT_BASE / cfg / task / "leaderboard.csv")
            if not p.exists():
                print(f"[skip] {p}"); continue
            df = pd.read_csv(p)
            df["key"] = df["source"] + "|" + df["covar_mode"]
            df["cell"] = [_fmt_pm(m, s) for m, s in zip(df[mcol], df[scol])]
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

    # Column order: source, covar_mode, then per-tag block
    tag_iter = ["strictQC"] if args.beta_source == "prscs" else LD_CONFIGS
    cols_ordered = ["source", "covar_mode"]
    for tag in tag_iter:
        cols_ordered.append(f"n_snps__{tag}")
        cols_ordered.append(f"bacc__{tag}")
        cols_ordered.append(f"cindex__{tag}")
        cols_ordered.append(f"R2__{tag}")
    cols_ordered = [c for c in cols_ordered if c in out.columns]
    out = out[cols_ordered]

    # Pretty column rename for display (TSV keeps the explicit names).
    # Three-line headers so each column can stay narrow: <window>\n<r²>\n<metric>
    SHORT_CFG = {
        "ld_1000kb_r2_0.8": ("1000 kb", "r²<0.8"),
        "ld_500kb_r2_0.2":  ("500 kb",  "r²<0.2"),
        "ld_250kb_r2_0.1":  ("250 kb",  "r²<0.1"),
        "ld_50_5_r2_0.5":   ("50 SNP",  "r²<0.5"),
        "strictQC":         ("strict-QC", "no LD prune"),
    }
    def _short(c: str) -> str:
        for cfg, (window, r2) in SHORT_CFG.items():
            if c.endswith(f"__{cfg}"):
                metric = c.split("__")[0]
                return f"{window}\n{r2}\n{metric}"
        return c

    # Sort by covar_mode (alpha) then by best classification BalAcc at default config
    sort_col = "bacc__strictQC" if args.beta_source == "prscs" else "bacc__ld_1000kb_r2_0.8"
    if sort_col in out.columns:
        def _bacc_val(v):
            if not isinstance(v, str) or v == "": return -1.0
            try:    return float(v.split("+/-")[0])
            except: return -1.0
        out["_sort"] = out[sort_col].apply(_bacc_val)
        out = out.sort_values(["covar_mode","_sort"], ascending=[True, False]).drop(columns=["_sort"])

    tsv_p = OUT_BASE / "master_leaderboard.tsv"
    out.to_csv(tsv_p, sep="\t", index=False)
    print(f"wrote master TSV: {tsv_p}")

    # PNG — wide landscape, multi-line headers, white header background,
    # top-per-column bolding only.
    n_rows = len(out)
    fig_w = 30
    fig_h = max(2.5, 0.32 * n_rows + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h)); ax.axis("off")
    ax.set_title("Strict-QC PRS — master leaderboard across LD configs and tasks   "
                 "(VAL set; mean +/- std over 3 seeds; bold = top value in column)",
                 fontsize=12, fontweight="bold", pad=10, loc="left")
    show = out.copy()
    for c in show.columns:
        if c.startswith("n_snps__"):
            show[c] = show[c].apply(lambda v: "" if pd.isna(v) else str(int(v)))

    display_headers = [_short(c) if "__" in c else c for c in show.columns]

    # Narrower metric columns (multi-line headers let columns stay tight).
    widths = []
    for c in show.columns:
        if c == "source":              widths.append(0.085)
        elif c == "covar_mode":        widths.append(0.070)
        elif c.startswith("n_snps"):   widths.append(0.020)
        else:                          widths.append(0.040)
    s = sum(widths); widths = [w / s for w in widths]

    tbl = ax.table(cellText=show.values.tolist(), colLabels=display_headers,
                    loc="center", cellLoc="center", colWidths=widths)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.45)
    # Header band: white bg, black bold text, larger font, COMPACT height
    HDR_HEIGHT = 0.060   # smaller = less vertical space taken by 3-line header
    for j in range(len(show.columns)):
        c = tbl[(0, j)]
        c.set_facecolor("white")
        c.set_text_props(color="black", weight="bold", fontsize=11)
        c.set_height(HDR_HEIGHT)
        c.set_edgecolor("black")
    # Bold the top value per metric column WITHIN EACH covar_mode block.
    # So each (column × covar_mode) cell-group has one bold winner.
    # For higher-is-better metrics (bacc, cindex, R2). Skip n_snps + ID cols.
    col_idx = {c: i for i, c in enumerate(show.columns)}
    def _mean(s: str) -> float:
        if not isinstance(s, str) or s == "": return float("-inf")
        try:    return float(s.split("+/-")[0])
        except: return float("-inf")
    out_idx = out.reset_index(drop=True)
    for col in show.columns:
        if not any(col.startswith(p) for p in ("bacc__", "cindex__", "R2__")):
            continue
        for mode in out_idx["covar_mode"].unique():
            mask = out_idx["covar_mode"] == mode
            sub = out_idx.loc[mask, col]
            if sub.empty: continue
            means = [_mean(v) for v in sub.tolist()]
            if all(m == float("-inf") for m in means): continue
            local_top = int(np.argmax(means))
            global_row = int(sub.index[local_top])
            tbl[(global_row + 1, col_idx[col])].set_text_props(weight="bold")

    png_p = OUT_BASE / "master_leaderboard.png"
    fig.savefig(png_p, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote master PNG: {png_p}")
    print("\nMaster leaderboard preview (first 12 rows):")
    print(out.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
