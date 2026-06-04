#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_longT4
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split/logs/enc_longT4_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split/logs/enc_longT4_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%4
# =============================================================================
# 03j_finetune_longitudinal_T4split_submit_csd3.sh
# =============================================================================
# Task 1b: re-extract the single all-visits longitudinal T2 encoder, but on the **T4-aware** folds
# (01j) so the downstream T4 horizon classifier is LEAKAGE-SAFE on the balanced baseline_T4 split.
# Identical to 03j_finetune_longitudinal_submit_csd3.sh EXCEPT data_dir + out root.
# SLURM_ARRAY_TASK_ID IS the seed (0/1/2).
#
# Prereqs:
#   - Build verbose/longitudinal_allvisits_T4split locally (run 01j) and rsync to RDS.
#   - git pull on the HPC (03j_finetune_longitudinal_extract.py + this submit present).
#   - mkdir -p ${OUT_DIR}/logs
# Submit:  sbatch clinical_pipeline/03j_finetune_longitudinal_T4split_submit_csd3.sh
# Smoke seed 0:  sbatch --array=0 clinical_pipeline/03j_finetune_longitudinal_T4split_submit_csd3.sh
# Pull back only: ${OUT_DIR}/<model>/T2_long_multiclass/seed_*/full_ft/{metrics.json,embeddings/}
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"
# Run EXPLICITLY in the 'clinical' env via `conda run` — independent of whatever env is active on the
# login node at submit time (no activate/inherit/stacking). Fail loudly if it's missing.
conda run -n clinical python -c "import numpy, torch, transformers" || {
    echo "ERROR: 'clinical' env missing/broken. conda env list:"; conda env list; exit 1; }

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/longitudinal_allvisits_T4split"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal_T4split"

MODEL_ID="thomas-sounack/BioClinical-ModernBERT-large"
MODEL_SLUG="${MODEL_ID##*/}"
SEED=${SLURM_ARRAY_TASK_ID}
METRICS="${OUT_DIR}/${MODEL_SLUG}/T2_long_multiclass/seed_${SEED}/full_ft/metrics.json"

echo "========================================================"
echo "  Longitudinal single-space T2 encoder — T4-AWARE folds (Task 1b)"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}  Seed: ${SEED}"
echo "  Data dir       : ${DATA_DIR}"
echo "  Output dir     : ${OUT_DIR}"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "========================================================"

if [ -f "${METRICS}" ]; then
    echo "  metrics.json exists — skipping."
    exit 0
fi

conda run -n clinical --no-capture-output python "${SCRIPT_DIR}/03j_finetune_longitudinal_extract.py" \
    --model_id  "${MODEL_ID}" \
    --seed      "${SEED}" \
    --data_dir  "${DATA_DIR}" \
    --out_dir   "${OUT_DIR}" \
    --hf_cache  "${HF_HOME}"

echo "  Finished."
