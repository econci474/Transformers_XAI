#!/bin/bash
#SBATCH -A COMPUTERLAB-SL2-GPU
#SBATCH --job-name=mamba_clf
# Logs + run outputs live UNDER a hpc-work root — NEVER in the repo/script folder.
# Create the logs dir before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/mamba_classifier/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/mamba_classifier/logs/mamba_clf_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/mamba_classifier/logs/mamba_clf_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# =============================================================================
# 04_train_classifier_submit_csd3.sh  (Task 2b — T4 3-class horizon classifier head)
# Runs the Sweep-A ablation (6 configs x 3 seeds) in one job (tiny: frozen 130M over 3 tokens, 146 pts).
# Predictions + logs land under ${CLF_ROOT} (hpc-work), not the repo.
# Then 05_classifier_ablation_summary.py is run SEPARATELY (env: clinical OR survml — only needs sklearn).
#
# Prereqs:
#   - T4-split embeddings on RDS:  encoder_outputs_..._longitudinal_T4split/.../seed_*/full_ft/embeddings/
#   - longitudinal master CSV on RDS (Label_T4 source).
#   - git pull (integration/T4/mamba_long_clinical present); mamba-1 weights cached in HF_HOME.
#   - mkdir -p ${CLF_ROOT}/logs
# Submit (full grid):  sbatch integration/T4/mamba_long_clinical/04_train_classifier_submit_csd3.sh
# Extra args pass THROUGH to the python (via "$@"), e.g. one config/seed:
#   sbatch integration/T4/mamba_long_clinical/04_train_classifier_submit_csd3.sh --only A_default_mamba1_frozen --seeds 0
# Pull back:  ${CLF_ROOT}/classifier/   (then run 04 with --pred_root pointing at it)
# =============================================================================
module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"
# Run EXPLICITLY in the 'clinical' env via `conda run` — independent of whatever env is active on the
# login node at submit time (no activate/inherit/stacking). Fail loudly if it's missing.
conda run -n clinical python -c "import numpy, torch, transformers, sklearn" || {
    echo "ERROR: 'clinical' env missing/broken. conda env list:"; conda env list; exit 1; }

REPO="/home/ec474/rds/hpc-work/Transformers_XAI"
DERIV="/home/ec474/rds/hpc-work/ADNI_CL"
CLF_ROOT="${DERIV}/mamba_classifier"          # logs + predictions root (NOT the repo)
EMB_ROOT="${DERIV}/encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split/BioClinical-ModernBERT-large/T2_long_multiclass"
MASTER="${DERIV}/no_cdr_stratified_post_exclusion/verbose/longitudinal/master_clinical_verbose.csv"

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
conda run -n clinical python "${REPO}/integration/T4/mamba_long_clinical/04_train_classifier_mamba.py" \
    --seeds 0 1 2 \
    --emb_root "${EMB_ROOT}" \
    --master   "${MASTER}" \
    --out_dir  "${CLF_ROOT}/classifier" \
    "$@" || { echo "ERROR: training job failed (see above)."; exit 1; }
echo "Finished. Predictions in ${CLF_ROOT}/classifier . Next:"
echo "  python ${REPO}/integration/T4/mamba_long_clinical/05_classifier_ablation_summary.py --pred_root ${CLF_ROOT}/classifier"
