"""
43b_render_prscs_comparison.py
==============================
Read per-source raw β (from <SRC>_full_snps_resolution.tsv) and PRS-CS
posterior β (from <SRC>_prscs_posterior_beta.tsv) for the 115 strict-QC SNPs.
Render a comparison TSV + PNG showing:
  - raw β vs posterior β scatter, per source
  - shrinkage magnitude (|posterior| / |raw|) per source
  - n_snps with posterior available
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path("D:/ADNI_SNP_Omni2.5M_20140220")
SRC_PRS_DIR  = ROOT / "source_prs"
POSTERIOR    = ROOT / "source_prs" / "prscs_posterior"
QC_STRICT    = ROOT / "GWAS_comprehensive_v2" / "QC_strict"
OUT = POSTERIOR / "raw_vs_prscs_comparison"

# Sources that ran PRS-CS successfully (matches 43_run_prscs_shrinkage.SOURCES)
PRSCS_SOURCES = ["Bellenguez", "Kunkle", "Kosteridis", "Schwanzentruber"]


def _load_raw(src: str) -> pd.DataFrame:
    """Lead-SNP resolution table with columns rsID, beta_A1 (or similar)."""
    p = SRC_PRS_DIR / f"{src}_full_snps_resolution.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    # Look for beta column
    beta_col = next((c for c in df.columns
                      if c in ("beta_A1","beta","Beta","BETA","GWAS_BETA")), None)
    if beta_col is None:
        return pd.DataFrame()
    df["raw_beta"] = pd.to_numeric(df[beta_col], errors="coerce")
    return df[["rsID", "raw_beta"]].dropna()


def _load_posterior(src: str) -> pd.DataFrame:
    p = POSTERIOR / f"{src}_prscs_posterior_beta.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", low_memory=False)
    # columns: chrom rsID bp A1 A2 posterior_beta
    return df[["rsID", "posterior_beta"]]


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # Target SNPs (strict-QC subset)
    bim = pd.read_csv(QC_STRICT / "recover_all_pool_strictQC.bim",
                       sep=r"\s+", header=None,
                       names=["chrom","rsID","cM","bp","A1","A2"])
    target = set(bim["rsID"])
    print(f"Target SNPs (strict-QC pool): {len(target)}")

    fig, axes = plt.subplots(1, len(PRSCS_SOURCES),
                              figsize=(4*len(PRSCS_SOURCES), 4.2), squeeze=False)
    summary_rows = []
    for j, src in enumerate(PRSCS_SOURCES):
        raw  = _load_raw(src)
        post = _load_posterior(src)
        if raw.empty or post.empty:
            print(f"  {src}: raw or posterior empty — skipping plot")
            axes[0, j].set_title(f"{src} (no data)"); axes[0, j].axis("off")
            summary_rows.append({"source": src, "n_overlap": 0})
            continue
        m = raw.merge(post, on="rsID")
        m = m[m["rsID"].isin(target)]
        ax = axes[0, j]
        if len(m):
            ax.scatter(m["raw_beta"], m["posterior_beta"], s=24, alpha=0.7)
            lim = max(m["raw_beta"].abs().max(), m["posterior_beta"].abs().max()) * 1.1
            ax.plot([-lim, lim], [-lim, lim], "--", color="grey", lw=0.6)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            shrink = (m["posterior_beta"].abs() / m["raw_beta"].abs().replace(0, np.nan)).median()
            ax.set_title(f"{src}  (n={len(m)})\n|post|/|raw| median = {shrink:.2f}")
        else:
            ax.set_title(f"{src} (0 overlap)")
        ax.axhline(0, color="grey", lw=0.4); ax.axvline(0, color="grey", lw=0.4)
        ax.set_xlabel("raw β (published lead)")
        if j == 0: ax.set_ylabel("PRS-CS posterior β")
        summary_rows.append({"source": src, "n_overlap": len(m),
                              "median_shrink_ratio": float(
                                  (m["posterior_beta"].abs() / m["raw_beta"].abs().replace(0, np.nan)).median()
                                  if len(m) else np.nan)})
        m.to_csv(OUT / f"{src}_raw_vs_prscs.tsv", sep="\t", index=False)

    fig.suptitle("Raw published β vs PRS-CS posterior β  —  strict-QC pool",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png = OUT / "raw_vs_prscs_comparison.png"
    fig.savefig(png, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {png}")

    sum_df = pd.DataFrame(summary_rows)
    sum_p = OUT / "raw_vs_prscs_summary.tsv"
    sum_df.to_csv(sum_p, sep="\t", index=False)
    print(f"wrote {sum_p}")
    print(sum_df.to_string(index=False))


if __name__ == "__main__":
    main()
