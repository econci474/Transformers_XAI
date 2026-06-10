#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=xmodal_mamba
# Logs live UNDER hpc-work, NEVER in the repo. Create the dir before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/xmodal_mamba/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/xmodal_mamba/logs/xmodal_mamba_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/xmodal_mamba/logs/xmodal_mamba_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# =============================================================================
# 02c_run_long_mamba_submit_csd3.sh  (Phase 2 — LONG arm, frozen-Mamba pooling)
# Runs the cross_attn_x LONG variant with --clin_pool mamba (54 runs: main lambda=0.01 + sweep),
# to compare frozen-Mamba pooling vs the Phase-1 mean pooling for the 12-month T2 fusion.
#
# WHY HPC: the frozen Mamba pooler reuses the T4 SeqEncoder (state-spaces/mamba-130m-hf via
# transformers). The sequential CPU scan is pathologically slow AND segfaults on long runs (verified
# locally), so this MUST run on a GPU where the CUDA scan/fast-kernels apply. _xmodal_lib auto-uses
# CUDA when torch.cuda.is_available().
#
# PREREQS (do these first):
#   - git push (local) + git pull (CSD3): integration/cross_modal_attention/* and
#     integration/T4/mamba_long_clinical/_mamba_seq_lib.py present under ${REPO}.
#   - The 'clinical' env has torch(+cuda), transformers, sklearn, pandas, pyarrow, numpy.
#   - The required embeddings are on RDS under the two roots below (the preflight checks them).
#   - mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/xmodal_mamba/logs
# SUBMIT:  sbatch integration/cross_modal_attention/02c_run_long_mamba_submit_csd3.sh
#          (extra args pass through, e.g. ... 02c..._submit_csd3.sh --overwrite)
# PULL BACK to local, then regenerate tables:
#   outputs/LONG/*/cross_attn_x/mamba/  and  outputs/sweep/*/LONG/*/cross_attn_x/mamba/
#   -> python integration/cross_modal_attention/03b_summary_lambda_sweep.py ; 03c_summary_joined.py
# =============================================================================
set -u
module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"     # Mamba weights cache (avoid home quota)

REPO="/home/ec474/rds/hpc-work/Transformers_XAI"
# --- Derivative roots on RDS (VERIFY these match your CSD3 layout) -------------------------------
# Clinical embeddings tree (encoder_outputs_no_cdr_post_exclusion{,_longitudinal,_m12}/...):
export XMODAL_CL_ROOT="/home/ec474/rds/hpc-work/ADNI_CL"
# MRI embeddings tree (mri_embeddings/... and brainmvp_embeddings/...):
export XMODAL_MRI_ROOT="/home/ec474/rds/hpc-work/ADNI_MRI"
# -----------------------------------------------------------------------------------------------

# Preflight: env import + CUDA + the exact files the LONG/mamba run needs for seed 0 (both variants).
conda run -n clinical python - <<'PY' || { echo "PREFLIGHT FAILED"; exit 1; }
import os, sys, importlib
for m in ("numpy", "torch", "transformers", "sklearn", "pandas", "pyarrow"):
    importlib.import_module(m)
import torch
assert torch.cuda.is_available(), "no CUDA visible to torch"
print("torch", torch.__version__, "| cuda", torch.version.cuda, "| GPU", torch.cuda.get_device_name(0))
cl, mri = os.environ["XMODAL_CL_ROOT"], os.environ["XMODAL_MRI_ROOT"]
# Required for ANY run: clinical bl + longitudinal (LONG arm) and the variant-B (BrainMVP) npz.
req = [
    os.path.join(cl, "encoder_outputs_no_cdr_post_exclusion", "BioClinical-ModernBERT-large",
                 "T2_multiclass", "seed_0", "full_ft", "embeddings", "embeddings.parquet"),
    os.path.join(cl, "encoder_outputs_no_cdr_post_exclusion_longitudinal", "BioClinical-ModernBERT-large",
                 "T2_long_multiclass", "seed_0", "full_ft", "embeddings", "embeddings.parquet"),
    os.path.join(mri, "brainmvp_embeddings", "T1b_binary", "aug_stochastic", "seed_0",
                 "embeddings_seed_0.npz"),
]
# Optional: variant-A (BrainDINO) npz. If absent, run with --variants B (the script still proceeds).
opt_A = os.path.join(mri, "mri_embeddings", "T2_multiclass", "braindino_frozen_none_cached",
                     "embeddings_seed_0.npz")
miss = [p for p in req if not os.path.exists(p)]
if miss:
    print("MISSING required embeddings (fix XMODAL_CL_ROOT / XMODAL_MRI_ROOT):")
    for p in miss:
        print("   ", p)
    sys.exit(1)
if not os.path.exists(opt_A):
    print("WARNING: variant-A (BrainDINO) npz not found — run with `--variants B` (BrainMVP only):")
    print("   ", opt_A)
print("preflight OK — CL + variant-B present" + (" (+ variant-A)" if os.path.exists(opt_A) else ""))
PY

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
conda run -n clinical python "${REPO}/integration/cross_modal_attention/02c_run_long_mamba.py" "$@" \
    || { echo "ERROR: LONG/mamba run failed (see above)."; exit 1; }
echo "Finished. Outputs under ${REPO}/integration/cross_modal_attention/outputs/{LONG,sweep}/ — pull back and run 03b/03c."
