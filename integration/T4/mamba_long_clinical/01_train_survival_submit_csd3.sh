#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=mamba_surv
# Logs + run outputs live UNDER a hpc-work root — NEVER in the repo/script folder.
# Create the logs dir before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/mamba_survival/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/mamba_survival/logs/mamba_surv_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/mamba_survival/logs/mamba_surv_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# =============================================================================
# 01_train_survival_submit_csd3.sh  (Task 2 — survival head ablation)
# Runs ALL configs × 3 seeds in one job (tiny: frozen 130M over 3 tokens, ~535 pts).
# Logs AND predictions are written under ${SURV_ROOT} (hpc-work), not the repo.
# Then 02_survival_comparison.py is run SEPARATELY in the `survml` env (sksurv/lifelines).
#
# Prereqs:
#   - Task-1 embeddings on RDS:  …/encoder_outputs_no_cdr_post_exclusion_longitudinal/…/seed_*/full_ft/embeddings/
#   - rsync the longitudinal master CSV to RDS (survival labels):
#       …/no_cdr_stratified_post_exclusion/verbose/longitudinal/master_clinical_verbose.csv
#   - git pull (integration/T4/mamba_long_clinical present); mamba-1/2 download to HF_HOME on first run.
#   - mkdir -p ${SURV_ROOT}/logs
#   - (optional, SODEN arm) pip install torchdiffeq into the clinical env; else B_soden is skipped.
# Submit:  sbatch integration/T4/mamba_long_clinical/01_train_survival_submit_csd3.sh
# Smoke one config:  add  --only A_default_mamba1_frozen --seeds 0  to the python call.
# Pull back:  ${SURV_ROOT}/survival/   (then run 02 with --pred_root pointing at it)
# =============================================================================
module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"
# Run EXPLICITLY in the 'clinical' env via `conda run` — independent of whatever env is active on the
# login node at submit time (no activate/inherit/stacking). Fail loudly if it's missing.
conda run -n clinical python -c "import numpy, torch, transformers" || {
    echo "ERROR: 'clinical' env missing/broken. conda env list:"; conda env list; exit 1; }

REPO="/home/ec474/rds/hpc-work/Transformers_XAI"
DERIV="/home/ec474/rds/hpc-work/ADNI_CL"
SURV_ROOT="${DERIV}/mamba_survival"            # logs + predictions root (NOT the repo)
EMB_ROOT="${DERIV}/encoder_outputs_no_cdr_post_exclusion_longitudinal/BioClinical-ModernBERT-large/T2_long_multiclass"
MASTER="${DERIV}/no_cdr_stratified_post_exclusion/verbose/longitudinal/master_clinical_verbose.csv"

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
conda run -n clinical python "${REPO}/integration/T4/mamba_long_clinical/01_train_survival_mamba.py" \
    --seeds 0 1 2 \
    --emb_root "${EMB_ROOT}" \
    --master   "${MASTER}" \
    --out_dir  "${SURV_ROOT}/survival" || { echo "ERROR: training job failed (see above)."; exit 1; }
echo "Finished. Predictions in ${SURV_ROOT}/survival . Next (survml env):"
echo "  conda run -n survml python ${REPO}/integration/T4/mamba_long_clinical/02_survival_comparison.py --pred_root ${SURV_ROOT}/survival"
