"""
01_fuse_T3.py — T3 per-horizon binary late fusion (clinical ⊕ MRI)
==================================================================
Thin orchestrator over the K-class fusion engine
`integration/T2_classification/01_fuse_T2.py`. Runs one binary fusion per
conversion horizon (3yr / 5yr / 7yr), each with its leaderboard-best clinical
encoder (baseline) and the per-task best cached-head MRI model (m12), framing
M12 (CL_bl + MRI_m12). Outputs land under `integration/T3/<horizon>/<arch>/`.

Per-horizon config (clinical pairs by HORIZON/label — note clinical names 7y as
`T3d_conv7y`, MRI as `T3c_conv7y`):

  horizon  clinical encoder (full_ft)        clinical task  MRI task      best MRI
  3yr      BioClinical-ModernBERT-large       T3a_conv3y     T3a_conv3y    ViT-MAE
  5yr      BioClinical-ModernBERT-base        T3b_conv5y     T3b_conv5y    BrainMVP
  7yr      ModernBERT-large                   T3d_conv7y     T3c_conv7y    ViT-MAE

Both MRI archs (ViT-MAE + BrainMVP cached) are run every horizon so the summary
can compare; the table notes the per-horizon best. Variants produced by the
engine: clinical_only, mri_only, weighted_avg (val-bACC), weighted_avg_J
(Youden), elastic_net, lr (+ 0.5 equal-weight is the weighted_avg w=0.5 point).

Run (needs an env with pandas+sklearn+matplotlib+pyarrow, e.g. snp):
  conda run -n snp python integration/T3/01_fuse_T3.py
  conda run -n snp python integration/T3/01_fuse_T3.py --horizons 3yr 7yr --archs vit_mae
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "T2_classification" / "01_fuse_T2.py"
DERIV = Path("D:/ADNI_BIDS_project/derivatives")
CLIN_ROOT = DERIV / "encoder_outputs_no_cdr_post_exclusion"
MRI_ROOT = DERIV / "mri_embeddings"
BMVP_FT_ROOT = DERIV / "brainmvp_embeddings"   # BrainMVP full_ft prob CSVs

# horizon -> (clinical model, clinical task, MRI task, label)
HORIZONS = {
    "3yr": dict(clin_model="BioClinical-ModernBERT-large", clin_task="T3a_conv3y",
                mri_task="T3a_conv3y", label="T3 conv<=3y"),
    "5yr": dict(clin_model="BioClinical-ModernBERT-base",  clin_task="T3b_conv5y",
                mri_task="T3b_conv5y", label="T3 conv<=5y"),
    "7yr": dict(clin_model="ModernBERT-large",             clin_task="T3d_conv7y",
                mri_task="T3c_conv7y", label="T3 conv<=7y"),
}

# MRI source key -> (display label, template builder taking mri_task -> Path with {seed}).
MRI_ARCHS = {
    "vit_mae":     ("ViT-MAE frozen/none",
                    lambda t: MRI_ROOT / t / "vit_mae_frozen_none_cached" / "embeddings_seed_{seed}.csv"),
    "brainmvp":    ("BrainMVP frozen/none",
                    lambda t: MRI_ROOT / t / "brainmvp_frozen_none_cached" / "embeddings_seed_{seed}.csv"),
    "brainmvp_ft": ("BrainMVP full_ft plus_original",
                    lambda t: BMVP_FT_ROOT / t / "aug_plus_original" / "seed_{seed}" / "embeddings_seed_{seed}.csv"),
}


def run_one(horizon: str, arch: str) -> int:
    cfg = HORIZONS[horizon]
    label, tmpl = MRI_ARCHS[arch]
    mri_template = tmpl(cfg["mri_task"])
    # skip (horizon, arch) where the extraction CSV doesn't exist (e.g. no full_ft for T3b)
    if not Path(str(mri_template).format(seed=0)).exists():
        print(f"  [skip] {horizon}/{arch}: no MRI CSV at {str(mri_template).format(seed=0)}")
        return 0
    out_dir = HERE / horizon
    cmd = [
        sys.executable, str(ENGINE),
        "--n-classes", "2", "--class-names", "noConv", "conv",
        "--clin-root", str(CLIN_ROOT),
        "--clin-model", cfg["clin_model"], "--clin-strat", "full_ft",
        "--clin-task", cfg["clin_task"],
        "--mri-name", label,
        "--out-name", arch,
        "--mri-csv-template", str(mri_template).replace("\\", "/"),
        "--framings", "M12",
        "--task-label", cfg["label"],
        "--out-dir", str(out_dir),
    ]
    print("\n" + "=" * 78)
    print(f"  T3 {horizon}  |  clinical={cfg['clin_model']}  |  MRI={label}")
    print("=" * 78)
    return subprocess.run(cmd).returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--horizons", nargs="+", default=list(HORIZONS), choices=list(HORIZONS))
    ap.add_argument("--archs", nargs="+", default=list(MRI_ARCHS), choices=list(MRI_ARCHS))
    args = ap.parse_args()

    if not ENGINE.exists():
        sys.exit(f"engine not found: {ENGINE}")
    fails = []
    for h in args.horizons:
        for a in args.archs:
            rc = run_one(h, a)
            if rc != 0:
                fails.append((h, a, rc))
    print("\n" + "=" * 78)
    if fails:
        print(f"  DONE with {len(fails)} failures: {fails}")
        sys.exit(1)
    print("  DONE — all T3 fusion configs succeeded. Outputs under integration/T3/<horizon>/<arch>/.")
    print("  Run 02_summary.py for the ranked cross-horizon table.")
    print("=" * 78)


if __name__ == "__main__":
    main()
