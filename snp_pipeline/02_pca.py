"""
pca.py – Visualise ADNI SNP PCA results as UMAP plots.

Inputs
------
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/SNP_filtered_with_mri_LD_pruned_pca.eigenvec
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/SNP_filtered_with_mri_LD_pruned_pca.eigenval
  D:/ADNI_SNP_Omni2.5M_20140220/phenotype.pheno      (FID IID PHENO;  1=CN 2=AD -9=missing)
  D:/ADNI_SNP_Omni2.5M_20140220/covariates.cov       (FID IID AGE SEX APOE4; SEX 1=M 2=F)
  D:/ADNI_SNP_Omni2.5M_20140220/results/missingness/missing.imiss   (subject missingness)
  D:/ADNI_SNP_Omni2.5M_20140220/results/maf/maf.frq                 (SNP-level MAF)
  D:/ADNI_SNP_Omni2.5M_20140220/results/hwe/hwe.hwe                 (SNP-level HWE p-values)
  D:/ADNI_BIDS_project/bids/genotype/subjects_with_snp_and_mri.tsv  (site per subject)

Outputs
-------
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/scree.png          ← QC
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/pc1_pc2.png        ← QC
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/maf_dist.png       ← QC
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/hwe_dist.png       ← QC
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/missingness.png    ← QC / UMAP
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/site.png           ← UMAP
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/diagnosis.png
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/age.png
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/sex.png
  D:/ADNI_SNP_Omni2.5M_20140220/results/pca/umaps/apoe4.png

Usage
-----
  pip install umap-learn matplotlib pandas numpy
  python pca.py
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")           # headless – no display needed on HPC
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import umap

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220")
PCA_STEM   = DATA_DIR / "results/pca_after_ld_r0.2/SNP_filtered_with_mri_LD_pruned_pca"
PHENO_FILE = DATA_DIR / "phenotype.pheno"
COVAR_FILE = DATA_DIR / "covariates.cov"
OUT_DIR    = DATA_DIR / "results/pca_after_ld_r0.2/umaps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load eigenvalues → variance explained per PC ──────────────────────────────
eigenval    = pd.read_csv(str(PCA_STEM) + ".eigenval", header=None, names=["eigenvalue"]).dropna()
total_var   = eigenval["eigenvalue"].sum()
pct_var     = (eigenval["eigenvalue"] / total_var * 100).round(2).tolist()
pc1_pct     = pct_var[0]   # % variance explained by PC1
pc2_pct     = pct_var[1]   # % variance explained by PC2

# ── Load eigenvec (PC scores) ─────────────────────────────────────────────────
eigenvec = pd.read_csv(str(PCA_STEM) + ".eigenvec", sep=r"\s+")
# columns: FID IID PC1 … PC20
pc_cols = [c for c in eigenvec.columns if c.startswith("PC")]

# ── Load phenotype ────────────────────────────────────────────────────────────
pheno = pd.read_csv(PHENO_FILE, sep=r"\s+")
pheno.columns = ["FID", "IID", "PHENO"]

# ── Load covariates ───────────────────────────────────────────────────────────
covar = pd.read_csv(COVAR_FILE, sep=r"\s+")
covar.columns = ["FID", "IID", "AGE", "SEX", "APOE4"]

# ── Merge on IID ──────────────────────────────────────────────────────────────
df = (
    eigenvec
    .merge(pheno[["IID", "PHENO"]], on="IID", how="left")
    .merge(covar[["IID", "AGE", "SEX", "APOE4"]], on="IID", how="left")
)

# ── UMAP on PC1 + PC2 ─────────────────────────────────────────────────────────
print("Running UMAP on PC1 & PC2 …")
reducer = umap.UMAP(
    n_components=2,
    n_neighbors=15,
    min_dist=0.3,
    metric="euclidean",
    random_state=42,
)
embedding = reducer.fit_transform(df[["PC1", "PC2"]].values)
df["UMAP1"] = embedding[:, 0]
df["UMAP2"] = embedding[:, 1]

# ── Shared aesthetics ─────────────────────────────────────────────────────────
FIGSIZE   = (8, 6)
DPI       = 180
PT_SIZE   = 18
PT_ALPHA  = 0.80
FONT      = "DejaVu Sans"
BG_COLOR  = "white"
AX_COLOR  = "white"
TICK_COL  = "#333333"
TITLE_COL = "black"
SPINE_COL = "#cccccc"

plt.rcParams.update({
    "font.family":      FONT,
    "figure.facecolor": BG_COLOR,
    "axes.facecolor":   AX_COLOR,
    "axes.edgecolor":   SPINE_COL,
    "axes.labelcolor":  TICK_COL,
    "xtick.color":      TICK_COL,
    "ytick.color":      TICK_COL,
    "grid.color":       SPINE_COL,
    "grid.linewidth":   0.5,
    "text.color":       TITLE_COL,
})

def base_fig(xlabel="UMAP 1", ylabel="UMAP 2"):
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_xlabel(xlabel, fontsize=11, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=11, labelpad=8)
    return fig, ax

def save(fig, name, title):
    fig.suptitle(title, fontsize=13, fontweight="bold", color=TITLE_COL, y=1.01)
    fig.tight_layout()
    path = OUT_DIR / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close(fig)
    print(f"  Saved → {path}")

# ─────────────────────────────────────────────────────────────────────────────
# QC-A. SCREE PLOT  – variance explained per PC
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig(xlabel="Principal Component", ylabel="Variance Explained (%)")

bar_colors = ["#7c83f5" if i < 2 else "#3a3d5c" for i in range(len(pct_var))]
ax.bar(range(1, len(pct_var) + 1), pct_var, color=bar_colors, edgecolor="none", width=0.7)
ax.set_xticks(range(1, len(pct_var) + 1))
ax.set_xticklabels([f"PC{i+1}" for i in range(len(pct_var))], rotation=45, ha="right", fontsize=8)

# Annotate PC1 and PC2 bars
for i in (0, 1):
    ax.text(i + 1, pct_var[i] + 0.15, f"{pct_var[i]:.1f}%",
            ha="center", va="bottom", fontsize=8, color="#c8c8e8")

ax.set_xlim(0.3, len(pct_var) + 0.7)
ax.set_ylim(0, max(pct_var) * 1.18)
save(fig, "scree.png", "Scree Plot: SNP PCA Variance Explained")

# ─────────────────────────────────────────────────────────────────────────────
# QC-B. RAW PC1 vs PC2 SCATTER  (coloured by diagnosis)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig(
    xlabel=f"PC1 ({pc1_pct:.1f}% variance explained)",
    ylabel=f"PC2 ({pc2_pct:.1f}% variance explained)",
)

for code, label in {**{1: "CN", 2: "AD"}, -9: "Missing"}.items():
    col  = {"CN": "#4fc3f7", "AD": "#ef5350", "Missing": "#555577"}.get(label, "#999999")
    mask = df["PHENO"] == code
    ax.scatter(
        df.loc[mask, "PC1"], df.loc[mask, "PC2"],
        c=col, label=f"{label} (n={mask.sum()})",
        s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
    )

legend = ax.legend(
    framealpha=0.25, facecolor=AX_COLOR, edgecolor=SPINE_COL,
    fontsize=9, markerscale=1.5, title="Diagnosis", title_fontsize=9,
)
legend.get_title().set_color(TITLE_COL)
save(fig, "pc1_pc2.png", "Raw PC1 vs PC2: QC Scatter")

# ─────────────────────────────────────────────────────────────────────────────
# 1. DIAGNOSIS
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig(
    xlabel=f"UMAP 1  (input: PC1 {pc1_pct:.1f}% + PC2 {pc2_pct:.1f}%)",
    ylabel="UMAP 2",
)

diag_map   = {1: "CN", 2: "AD"}
diag_color = {"CN": "#4fc3f7", "AD": "#ef5350", "Missing": "#555577"}

for code, label in {**diag_map, -9: "Missing"}.items():
    col   = diag_color.get(label, "#999999")
    mask  = df["PHENO"] == code
    ax.scatter(
        df.loc[mask, "UMAP1"], df.loc[mask, "UMAP2"],
        c=col, label=f"{label} (n={mask.sum()})",
        s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
    )

legend = ax.legend(
    framealpha=0.85, facecolor=AX_COLOR, edgecolor=SPINE_COL,
    fontsize=9, markerscale=1.5, title="Diagnosis", title_fontsize=9,
)
legend.get_title().set_color(TITLE_COL)
save(fig, "diagnosis.png", "UMAP of SNP PCs: Coloured by Diagnosis")

# ─────────────────────────────────────────────────────────────────────────────
# 2. AGE  (continuous colourmap)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig()

valid_age = df["AGE"].notna()
sc = ax.scatter(
    df.loc[valid_age, "UMAP1"], df.loc[valid_age, "UMAP2"],
    c=df.loc[valid_age, "AGE"],
    cmap="plasma", s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
)
if (~valid_age).any():
    ax.scatter(
        df.loc[~valid_age, "UMAP1"], df.loc[~valid_age, "UMAP2"],
        c="#555577", s=PT_SIZE, alpha=0.5, linewidths=0, label="Missing",
    )
    ax.legend(framealpha=0.25, facecolor=AX_COLOR, edgecolor=SPINE_COL,
              fontsize=9, markerscale=1.5)

cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label("Age (years)", color=TICK_COL, fontsize=10)
cbar.ax.yaxis.set_tick_params(color=TICK_COL)
plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=TICK_COL)
save(fig, "age.png", "UMAP of SNP PCs: Coloured by Age")

# ─────────────────────────────────────────────────────────────────────────────
# 3. SEX  (categorical: 1=Male, 2=Female)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig()

sex_map   = {1: "Male", 2: "Female"}
sex_color = {"Male": "#42a5f5", "Female": "#ec407a", "Missing": "#555577"}

for code, label in {**sex_map, -9: "Missing"}.items():
    mask = (df["SEX"] == code) if code != -9 else df["SEX"].isna()
    if not mask.any():
        continue
    ax.scatter(
        df.loc[mask, "UMAP1"], df.loc[mask, "UMAP2"],
        c=sex_color.get(label, "#999999"),
        label=f"{label} (n={mask.sum()})",
        s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
    )

legend = ax.legend(
    framealpha=0.25, facecolor=AX_COLOR, edgecolor=SPINE_COL,
    fontsize=9, markerscale=1.5, title="Sex", title_fontsize=9,
)
legend.get_title().set_color(TITLE_COL)
save(fig, "sex.png", "UMAP of SNP PCs: Coloured by Sex")

# ─────────────────────────────────────────────────────────────────────────────
# 4. APOE4  (categorical: 0 / 1 / 2 copies)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig()

apoe_color = {
    0.0: "#66bb6a",   # 0 copies – green
    1.0: "#ffa726",   # 1 copy   – orange
    2.0: "#ef5350",   # 2 copies – red
}

df["APOE4"] = pd.to_numeric(df["APOE4"], errors="coerce")
apoe_values = sorted(df["APOE4"].dropna().unique())

for val in apoe_values:
    mask  = df["APOE4"] == val
    label = f"ε4×{int(val)}"
    ax.scatter(
        df.loc[mask, "UMAP1"], df.loc[mask, "UMAP2"],
        c=apoe_color.get(val, "#aaaaaa"),
        label=f"{label} (n={mask.sum()})",
        s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
    )

missing_apoe = df["APOE4"].isna()
if missing_apoe.any():
    ax.scatter(
        df.loc[missing_apoe, "UMAP1"], df.loc[missing_apoe, "UMAP2"],
        c="#555577", label=f"Missing (n={missing_apoe.sum()})",
        s=PT_SIZE, alpha=0.5, linewidths=0,
    )

legend = ax.legend(
    framealpha=0.25, facecolor=AX_COLOR, edgecolor=SPINE_COL,
    fontsize=9, markerscale=1.5, title="APOE4 copies", title_fontsize=9,
)

legend.get_title().set_color(TITLE_COL)
save(fig, "apoe4.png", "UMAP of SNP PCs: Coloured by APOE4 Dosage")

# ─────────────────────────────────────────────────────────────────────────────
# Load additional QC files
# ─────────────────────────────────────────────────────────────────────────────
MISS_DIR = DATA_DIR / "results/missingness"
MAF_DIR  = DATA_DIR / "results/maf"
HWE_DIR  = DATA_DIR / "results/hwe"
SITE_TSV = pathlib.Path("D:/ADNI_BIDS_project/bids/genotype/subjects_with_snp_and_mri.tsv")

# Subject-level missingness (.imiss)
imiss = pd.read_csv(MISS_DIR / "missing.imiss", sep=r"\s+")
imiss.columns = imiss.columns.str.strip()
# IID in imiss is in 002_S_0413 format – matches eigenvec IID directly
df = df.merge(imiss[["IID", "F_MISS"]], on="IID", how="left")

# Site (from subjects_with_snp_and_mri.tsv)
site_df = pd.read_csv(SITE_TSV, sep="\t")
# adni_subject_id is '002_S_0413', convert to match IID format
site_df["IID"] = site_df["adni_subject_id"]
df = df.merge(site_df[["IID", "site"]], on="IID", how="left")

# ─────────────────────────────────────────────────────────────────────────────
# QC-C. MAF DISTRIBUTION  (SNP-level histogram)
# ─────────────────────────────────────────────────────────────────────────────
print("Loading MAF data (large file, may take a moment) …")
maf_df = pd.read_csv(MAF_DIR / "maf.frq", sep=r"\s+")
maf_df.columns = maf_df.columns.str.strip()
maf_vals = pd.to_numeric(maf_df["MAF"], errors="coerce").dropna()
maf_vals = maf_vals[maf_vals > 0]   # drop monomorphic SNPs

fig, ax = base_fig(xlabel="Minor Allele Frequency (MAF)", ylabel="Number of SNPs")
ax.hist(maf_vals, bins=50, color="#5c85d6", edgecolor="none", alpha=0.85)
ax.axvline(0.01, color="#ef5350", linewidth=1.2, linestyle="--", label="MAF = 0.01")
ax.axvline(0.05, color="#ffa726", linewidth=1.2, linestyle="--", label="MAF = 0.05")
ax.legend(framealpha=0.85, facecolor=AX_COLOR, edgecolor=SPINE_COL, fontsize=9)
ax.text(0.98, 0.96, f"n = {len(maf_vals):,} SNPs",
        transform=ax.transAxes, ha="right", va="top", fontsize=9, color=TICK_COL)
save(fig, "maf_dist.png", "QC: MAF Distribution")

# ─────────────────────────────────────────────────────────────────────────────
# QC-D. HWE DISTRIBUTION  (–log10 p-value histogram)
# ─────────────────────────────────────────────────────────────────────────────
print("Loading HWE data (large file, may take a moment) …")
hwe_df = pd.read_csv(HWE_DIR / "hwe.hwe", sep=r"\s+")
hwe_df.columns = hwe_df.columns.str.strip()
hwe_p = pd.to_numeric(hwe_df["P"], errors="coerce").dropna()
hwe_p = hwe_p[(hwe_p > 0) & (hwe_p <= 1)]
log_p  = -np.log10(hwe_p)

fig, ax = base_fig(xlabel="–log₁₀(HWE p-value)", ylabel="Number of SNPs")
ax.hist(log_p, bins=60, color="#7c83f5", edgecolor="none", alpha=0.85)
ax.axvline(-np.log10(1e-6), color="#ef5350", linewidth=1.2, linestyle="--",
           label="p = 1×10⁻⁶ (typical threshold)")
ax.legend(framealpha=0.85, facecolor=AX_COLOR, edgecolor=SPINE_COL, fontsize=9)
ax.text(0.98, 0.96, f"n = {len(hwe_p):,} SNPs",
        transform=ax.transAxes, ha="right", va="top", fontsize=9, color=TICK_COL)
save(fig, "hwe_dist.png", "QC: HWE –log₁₀(p) Distribution")

# ─────────────────────────────────────────────────────────────────────────────
# 5. MISSINGNESS  (UMAP coloured by per-subject call rate)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig(
    xlabel=f"UMAP 1  (input: PC1 {pc1_pct:.1f}% + PC2 {pc2_pct:.1f}%)",
    ylabel="UMAP 2",
)

valid_miss = df["F_MISS"].notna()
sc_miss = ax.scatter(
    df.loc[valid_miss, "UMAP1"], df.loc[valid_miss, "UMAP2"],
    c=df.loc[valid_miss, "F_MISS"],
    cmap="YlOrRd", s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
)
if (~valid_miss).any():
    ax.scatter(
        df.loc[~valid_miss, "UMAP1"], df.loc[~valid_miss, "UMAP2"],
        c="#aaaaaa", s=PT_SIZE, alpha=0.5, linewidths=0, label="Missing",
    )
    ax.legend(framealpha=0.85, facecolor=AX_COLOR, edgecolor=SPINE_COL, fontsize=9)
cbar = fig.colorbar(sc_miss, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label("SNP missingness rate (F_MISS)", color=TICK_COL, fontsize=9)
cbar.ax.yaxis.set_tick_params(color=TICK_COL)
plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=TICK_COL)
save(fig, "missingness.png", "UMAP of SNP PCs: Coloured by Subject Missingness")

# ─────────────────────────────────────────────────────────────────────────────
# 6. SITE  (UMAP coloured by acquisition site)
# ─────────────────────────────────────────────────────────────────────────────
fig, ax = base_fig(
    xlabel=f"UMAP 1  (input: PC1 {pc1_pct:.1f}% + PC2 {pc2_pct:.1f}%)",
    ylabel="UMAP 2",
)

unique_sites = sorted(df["site"].dropna().unique().astype(int))
# Use a qualitative colourmap (tab20 handles up to 20 categories well)
cmap_sites = plt.cm.get_cmap("tab20", len(unique_sites))
site_color  = {s: cmap_sites(i) for i, s in enumerate(unique_sites)}

for s in unique_sites:
    mask = df["site"] == s
    ax.scatter(
        df.loc[mask, "UMAP1"], df.loc[mask, "UMAP2"],
        color=site_color[s],
        label=f"Site {s} (n={mask.sum()})",
        s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
    )
missing_site = df["site"].isna()
if missing_site.any():
    ax.scatter(
        df.loc[missing_site, "UMAP1"], df.loc[missing_site, "UMAP2"],
        c="#aaaaaa", label=f"Missing (n={missing_site.sum()})",
        s=PT_SIZE, alpha=0.5, linewidths=0,
    )

legend = ax.legend(
    framealpha=0.85, facecolor=AX_COLOR, edgecolor=SPINE_COL,
    fontsize=7, markerscale=1.4, title="Site", title_fontsize=8,
    ncol=2, loc="upper right",
)
legend.get_title().set_color(TITLE_COL)
save(fig, "site.png", "UMAP of SNP PCs: Coloured by Acquisition Site")

print("\nDone! All plots saved to", OUT_DIR)