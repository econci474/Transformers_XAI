"""
07_analyse_zscore_distributions.py (gwas_signed_regression)
============================================================
Compare z-score distributions across BMFM regression input datasets.

Datasets compared:
  08c  - Base single-SNP forward-only dataset
  08d  - Forward + reverse complement (single SNP, 2x)
  08e  - Chromosome-level 8kb windows, forward only
  08ef - Chromosome-level 8kb windows, forward + reverse
  08f  - Combinatorial EA/OA (if generated)

Usage:
  python snp_pipeline/gwas_signed_regression/07_analyse_zscore_distributions.py
"""
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

BASE = pathlib.Path("D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs")

DATASETS = {
    "08c: single-SNP forward": BASE / "bmfm_gwas_signed_regression_without_ukb",
    "08d: single-SNP fwd+RC":  BASE / "bmfm_gwas_signed_regression_without_ukb_augmented",
    "08e: by-chrom forward":   BASE / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom" / "forward",
    "08e: by-chrom fwd+RC":    BASE / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom" / "forward_and_reverse",
    "08f: combinatorial fwd":  BASE / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos" / "forward",
    "08f: combinatorial fwd+RC": BASE / "bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos" / "forward_and_reverse",
}

OUT_DIR = pathlib.Path("outputs/snp_pipeline/zscore_analysis")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SPLITS = ["train", "dev", "test"]


def load_zscores(dataset_dir: pathlib.Path) -> pd.Series | None:
    """Load and concatenate z_score across all splits."""
    parts = []
    for split in SPLITS:
        p = dataset_dir / f"{split}.csv"
        if p.exists():
            df = pd.read_csv(p, usecols=["z_score"])
            parts.append(df["z_score"])
    if not parts:
        return None
    return pd.concat(parts, ignore_index=True)


# ── Collect data ──────────────────────────────────────────────────────────────
print("\nZ-score distribution summary")
print("=" * 75)
print(f"{'Dataset':<30} {'N':>7} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8} {'|z|>5':>6}")
print("-" * 75)

all_data = {}
for name, path in DATASETS.items():
    zs = load_zscores(path)
    if zs is None:
        print(f"{name:<30}  (not found — skipping)")
        continue
    all_data[name] = zs
    extreme = (zs.abs() > 5).sum()
    print(f"{name:<30} {len(zs):>7,} {zs.mean():>8.2f} {zs.std():>8.2f} "
          f"{zs.min():>8.2f} {zs.max():>8.2f} {extreme:>6}")

print("=" * 75)

# ── Plot ──────────────────────────────────────────────────────────────────────
n_plots = len(all_data)
if n_plots == 0:
    print("No datasets found — check paths.")
    raise SystemExit(1)

fig, axes = plt.subplots(1, n_plots, figsize=(4 * n_plots, 5), sharey=False)
if n_plots == 1:
    axes = [axes]

colors = plt.cm.tab10(np.linspace(0, 0.9, n_plots))

for ax, (name, zs), color in zip(axes, all_data.items(), colors):
    ax.hist(zs, bins=50, color=color, alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_title(name, fontsize=8)
    ax.set_xlabel("z-score")
    ax.set_ylabel("Count" if ax == axes[0] else "")
    ax.text(0.97, 0.97,
            f"N={len(zs):,}\nμ={zs.mean():.2f}\nσ={zs.std():.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

fig.suptitle("BMFM GWAS Regression — Z-score Distributions by Dataset", fontsize=11, y=1.02)
plt.tight_layout()

out_png = OUT_DIR / "zscore_distributions.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nPlot saved → {out_png}")
plt.show()
