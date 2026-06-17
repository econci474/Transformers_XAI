"""
08_plot_training_curves.py (gwas_signed_regression)
====================================================
Read TFEvents files from a Lightning log directory and plot train/val loss.
Saves metrics to JSON and plots to PNG.

Does NOT launch the TensorBoard server — avoids the pkg_resources issue.

Usage
-----
  python 08_plot_training_curves.py
  python 08_plot_training_curves.py --logdir path/to/lightning_logs --out results/
"""
import os
import json
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

LOGDIR_DEFAULT = (
    r"D:\ADNI_SNP_Omni2.5M_20140220\bmfm_outputs_regression"
    r"\bmfm_gwas_signed_regression_without_ukb\lightning_logs\lightning_logs"
)
OUT_DEFAULT = r"C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\outputs\snp_pipeline"

parser = argparse.ArgumentParser()
parser.add_argument("--logdir", default=LOGDIR_DEFAULT)
parser.add_argument("--out",    default=OUT_DEFAULT, help="Output directory for JSON + PNG")
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

# ── Import only the backend (avoids pkg_resources import in tensorboard.main) ─
from tensorboard.backend.event_processing import event_accumulator  # noqa

# ── Read all versions ──────────────────────────────────────────────────────────
all_data = {}   # {version_name: {tag: [(step, value), ...]}}

for version in sorted(os.listdir(args.logdir)):
    vpath = os.path.join(args.logdir, version)
    if not os.path.isdir(vpath):
        continue
    ea = event_accumulator.EventAccumulator(vpath)
    ea.Reload()
    tags = ea.Tags().get("scalars", [])
    if not tags:
        print(f"  {version}: no scalar tags found — skipping")
        continue
    print(f"\n=== {version} ===  tags: {sorted(tags)}")
    all_data[version] = {}
    for tag in sorted(tags):
        events = ea.Scalars(tag)
        vals   = [(e.step, float(e.value)) for e in events]
        all_data[version][tag] = vals
        print(f"  {tag}: {[(s, round(v, 3)) for s, v in vals]}")

if not all_data:
    print("No data found — check your --logdir path.")
    raise SystemExit(1)

# ── Save to JSON ───────────────────────────────────────────────────────────────
json_path = os.path.join(args.out, "metrics.json")
with open(json_path, "w") as f:
    json.dump(all_data, f, indent=2)
print(f"\nMetrics saved → {json_path}")

# ── Plot: train + val only (2 panels); test MSE goes in title ───────────────
LOSS_KEYWORDS = {
    "train": ["train/loss", "train_loss", "train/z_score"],
    "val"  : ["val/loss",   "val_loss",   "validation/loss", "val/z_score",
              "validation/z_score"],
}

# Extract test MSE for title (from version_2 if present)
test_mse = None
for version, tags in all_data.items():
    for tag, vals in tags.items():
        if "test" in tag.lower() and "mse" in tag.lower() and vals:
            test_mse = vals[-1][1]   # last value

title = "BMFM GWAS Regression — Training Curves (MSE loss)"
if test_mse is not None:
    title += f"  [Test MSE: {test_mse:.2f}  |  Test RMSE: {test_mse**0.5:.2f}]"

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
fig.suptitle(title, fontsize=11, y=1.02)

colors = cm.tab10(np.linspace(0, 0.5, max(len(all_data), 1)))

for ax, (split, keywords) in zip(axes, LOSS_KEYWORDS.items()):
    ax.set_title(f"{split.capitalize()} Loss")
    ax.set_xlabel("Step")
    ax.set_ylabel("MSE Loss")
    ax.grid(True, alpha=0.3)
    plotted = False
    for i, (version, tags) in enumerate(all_data.items()):
        for tag, vals in tags.items():
            if any(kw in tag.lower() for kw in keywords):
                steps, values = zip(*vals)
                ax.plot(steps, values, marker="o", markersize=3,
                        color=colors[i], label=f"{version} — {tag}")
                plotted = True
    if not plotted:
        for i, (version, tags) in enumerate(all_data.items()):
            for tag, vals in tags.items():
                if split in tag.lower():
                    steps, values = zip(*vals)
                    ax.plot(steps, values, marker="o", markersize=3,
                            color=colors[i], label=f"{version} — {tag}")
    ax.legend(fontsize=7, loc="upper right")

plt.tight_layout()
png_path = os.path.join(args.out, "training_curves.png")
plt.savefig(png_path, dpi=150, bbox_inches="tight")
print(f"Plot saved → {png_path}")
plt.show()
