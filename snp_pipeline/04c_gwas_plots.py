"""
gwas_plots.py – Manhattan and QQ plots for ADNI GWAS results.

Inputs
------
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/gwas_CN_vs_AD_no_pruning.assoc.logistic

Outputs
-------
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/plots/manhattan.png
  D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/plots/qq.png

Usage
-----
  pip install pandas numpy matplotlib scipy
  python gwas_plots.py
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import chi2

# ── Paths ─────────────────────────────────────────────────────────────────────
GWAS_FILE = pathlib.Path(
    "D:/ADNI_SNP_Omni2.5M_20140220/results/gwas_no_ld/"
    "gwas_CN_vs_AD_no_pruning.assoc.logistic"
)
OUT_DIR = GWAS_FILE.parent / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Style ─────────────────────────────────────────────────────────────────────
FIGSIZE    = (14, 5)
DPI        = 180
FONT       = "DejaVu Sans"
BG_COLOR   = "white"
TICK_COL   = "#333333"
SPINE_COL  = "#cccccc"
TITLE_COL  = "black"

# Significance thresholds
P_GW   = 5e-8    # genome-wide significance
P_SUG  = 1e-5    # suggestive significance (user's threshold)

# Alternating chromosome colours (muted blue / muted grey-blue)
CHR_COLORS = ["#4a7cb5", "#8eadd4"]
HIT_COLOR  = "#e63946"   # red for hits above suggestive line

plt.rcParams.update({
    "font.family":    FONT,
    "figure.facecolor": BG_COLOR,
    "axes.facecolor":   BG_COLOR,
    "axes.edgecolor":   SPINE_COL,
    "axes.labelcolor":  TITLE_COL,
    "xtick.color":      TICK_COL,
    "ytick.color":      TICK_COL,
    "text.color":       TITLE_COL,
    "grid.color":       SPINE_COL,
    "grid.linewidth":   0.5,
})


# ─────────────────────────────────────────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────────────────────────────────────────
print("Loading GWAS results (large file, this may take 30–60 s) …")
df = pd.read_csv(GWAS_FILE, sep=r"\s+")
df.columns = df.columns.str.strip()

# Keep ADD test only, drop NA p-values
df = df[df["TEST"] == "ADD"].copy()
df["P"] = pd.to_numeric(df["P"], errors="coerce")
df = df.dropna(subset=["P"])
df = df[df["P"] > 0]

# Chromosome as integer (drop X/Y/MT for now)
df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
df = df.dropna(subset=["CHR"])
df["CHR"] = df["CHR"].astype(int)
df["BP"]  = pd.to_numeric(df["BP"], errors="coerce")
df = df.dropna(subset=["BP"])
df["BP"]  = df["BP"].astype(int)

df["-log10P"] = -np.log10(df["P"])

print(f"  Retained {len(df):,} SNPs after filtering.")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Compute cumulative x-axis positions for Manhattan plot
# ─────────────────────────────────────────────────────────────────────────────
df = df.sort_values(["CHR", "BP"])

# Vectorised: compute per-chromosome BP minimum once, then map
chrom_bp_min = df.groupby("CHR")["BP"].transform("min")
chrom_bp_max = df.groupby("CHR")["BP"].max()

chrom_offsets = {}
chrom_mids    = {}
offset = 0
for chrom in sorted(df["CHR"].unique()):
    span = chrom_bp_max[chrom] - df.loc[df["CHR"] == chrom, "BP"].min()
    chrom_offsets[chrom] = offset
    chrom_mids[chrom]    = offset + span / 2
    offset += span + 5_000_000   # 5 Mb gap between chromosomes

# Map each row's chromosome to its offset, then add (BP - chr_min) — fully vectorised
offset_series = df["CHR"].map(chrom_offsets)
df["x_pos"]   = offset_series + (df["BP"] - chrom_bp_min)

# ─────────────────────────────────────────────────────────────────────────────
# 3. Manhattan plot
# ─────────────────────────────────────────────────────────────────────────────
print("Plotting Manhattan …")
fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
fig.patch.set_facecolor(BG_COLOR)
ax.set_facecolor(BG_COLOR)

chroms = sorted(df["CHR"].unique())
for i, chrom in enumerate(chroms):
    sub = df[df["CHR"] == chrom]
    # Colour hits (above suggestive line) bright red, rest alternating
    colors = np.where(sub["-log10P"] >= -np.log10(P_SUG), HIT_COLOR, CHR_COLORS[i % 2])
    ax.scatter(sub["x_pos"], sub["-log10P"],
               c=colors, s=3, alpha=0.75, linewidths=0, rasterized=True)

# Significance lines
gw_y  = -np.log10(P_GW)
sug_y = -np.log10(P_SUG)
ax.axhline(gw_y,  color="#ff006e", linewidth=1.0, linestyle="--",
           label=f"Genome-wide (p=5×10⁻⁸)")
ax.axhline(sug_y, color="#fb8500", linewidth=0.8, linestyle=":",
           label=f"Suggestive (p=1×10⁻⁵)")

# Annotate the top hits (above suggestive threshold)
hits = df[df["-log10P"] >= sug_y].sort_values("-log10P", ascending=False)
print(f"  Hits above p<1e-5: {len(hits)}")
for _, row in hits.iterrows():
    ax.annotate(
        row["SNP"],
        xy=(row["x_pos"], row["-log10P"]),
        xytext=(0, 6), textcoords="offset points",
        fontsize=6, color=HIT_COLOR, ha="center", va="bottom",
        fontweight="bold",
        arrowprops=dict(arrowstyle="-", color=HIT_COLOR, lw=0.5),
    )

# x-axis: chromosome labels
ax.set_xticks([chrom_mids[c] for c in chroms])
ax.set_xticklabels([str(c) for c in chroms], fontsize=7)
ax.set_xlabel("Chromosome", fontsize=10)
ax.set_ylabel("–log₁₀(p-value)", fontsize=10)
ax.set_title("Manhattan Plot: CN vs AD (GWAS, no LD pruning)", fontsize=12,
             color=TITLE_COL, fontweight="bold")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color(SPINE_COL)
ax.spines["bottom"].set_color(SPINE_COL)
ax.set_xlim(left=0)

legend = ax.legend(framealpha=0.85, facecolor=BG_COLOR, edgecolor=SPINE_COL, fontsize=8)

out_path = OUT_DIR / "manhattan.png"
fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
plt.close(fig)
print(f"  Saved → {out_path}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. QQ plot
# ─────────────────────────────────────────────────────────────────────────────
print("Plotting QQ …")

p_vals = df["P"].values
n      = len(p_vals)
observed  = np.sort(-np.log10(p_vals))[::-1]
expected  = -np.log10(np.arange(1, n + 1) / n)

# Genomic inflation factor λ
chi2_obs    = chi2.ppf(1 - p_vals, df=1)
lambda_gc   = np.median(chi2_obs) / chi2.ppf(0.5, df=1)

fig, ax = plt.subplots(figsize=(6, 6), dpi=DPI)
fig.patch.set_facecolor(BG_COLOR)
ax.set_facecolor(BG_COLOR)

ax.scatter(expected, observed, s=4, alpha=0.6, color="#4a7cb5",
           linewidths=0, rasterized=True, label="Observed")
max_val = max(observed.max(), expected.max()) * 1.05
ax.plot([0, max_val], [0, max_val], color="#e63946", linewidth=1.2,
        linestyle="--", label="Expected (null)")

ax.set_xlabel("Expected –log₁₀(p)", fontsize=10)
ax.set_ylabel("Observed –log₁₀(p)", fontsize=10)
ax.set_title("QQ Plot: CN vs AD", fontsize=12, color=TITLE_COL, fontweight="bold")

ax.text(0.05, 0.92, f"λ (genomic inflation) = {lambda_gc:.3f}",
        transform=ax.transAxes, fontsize=9, color=TICK_COL,
        bbox=dict(facecolor=BG_COLOR, edgecolor=SPINE_COL, boxstyle="round,pad=0.3"))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color(SPINE_COL)
ax.spines["bottom"].set_color(SPINE_COL)
ax.legend(framealpha=0.85, facecolor=BG_COLOR, edgecolor=SPINE_COL, fontsize=8)

out_path = OUT_DIR / "qq.png"
fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
plt.close(fig)
print(f"  Saved → {out_path}")

print(f"\nDone! Plots saved to {OUT_DIR}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Print top hits summary table
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n{'─'*70}")
print(f"Top SNPs (p < 1×10⁻⁵):")
print(f"{'─'*70}")
top = (df[df["P"] < P_SUG]
       [["CHR", "SNP", "BP", "A1", "OR", "SE", "P"]]
       .sort_values("P"))
if top.empty:
    print("  None found.")
else:
    print(top.to_string(index=False))
print(f"{'─'*70}")
print(f"Genomic inflation factor λ = {lambda_gc:.4f}")
print("  λ ≈ 1.00 → no inflation (good)")
print("  λ >> 1   → possible population stratification or confounding")
