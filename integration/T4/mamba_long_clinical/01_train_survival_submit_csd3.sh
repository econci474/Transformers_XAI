#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=mamba_surv
#   mkdir -p /home/ec474/rds/hpc-work/Transformers_XAI/integration/T4/mamba_long_clinical/logs
#SBATCH --output=/home/ec474/rds/hpc-work/Transformers_XAI/integration/T4/mamba_long_clinical/logs/mamba_surv_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/Transformers_XAI/integration/T4/mamba_long_clinical/logs/mamba_surv_%j.err
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
# Then 02_survival_comparison.py is run SEPARATELY in the `survml` env (sksurv/lifelines).
#
# Prereqs:
#   - Task-1 embeddings on RDS:  …/encoder_outputs_no_cdr_post_exclusion_longitudinal/…/seed_*/full_ft/embeddings/
#   - rsync the longitudinal master CSV to RDS (survival labels):
#       …/no_cdr_stratified_post_exclusion/verbose/longitudinal/master_clinical_verbose.csv
#   - git pull (integration/T4/mamba_long_clinical present); mamba-1/2 download to HF_HOME on first run.
#   - (optional, SODEN arm) pip install torchdiffeq into the clinical env; else B_soden is skipped.
# Submit:  sbatch integration/T4/mamba_long_clinical/01_train_survival_submit_csd3.sh
# Smoke one config:  add  --only A_default_mamba1_frozen --seeds 0  to the python call.
# Pull back:  integration/T4/mamba_long_clinical/outputs/survival/
# =============================================================================
module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

REPO="/home/ec474/rds/hpc-work/Transformers_XAI"
DERIV="/home/ec474/rds/hpc-work/ADNI_CL"
EMB_ROOT="${DERIV}/encoder_outputs_no_cdr_post_exclusion_longitudinal/BioClinical-ModernBERT-large/T2_long_multiclass"
MASTER="${DERIV}/no_cdr_stratified_post_exclusion/verbose/longitudinal/master_clinical_verbose.csv"

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
python "${REPO}/integration/T4/mamba_long_clinical/01_train_survival_mamba.py" \
    --seeds 0 1 2 \
    --emb_root "${EMB_ROOT}" \
    --master   "${MASTER}"
echo "Finished. Next: conda run -n survml python ${REPO}/integration/T4/mamba_long_clinical/02_survival_comparison.py"
