#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_long
# Logs land UNDER the longitudinal output directory. Slurm does NOT expand shell vars here, so
# this absolute path must match OUT_DIR/logs below. Create it before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal/logs/enc_long_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal/logs/enc_long_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%4
# =============================================================================
# 03j_finetune_longitudinal_submit_csd3.sh
# =============================================================================
# TASK 1 of the T4 longitudinal-MAMBA arm. Fine-tune ONE BioClinical-ModernBERT-large model on
# ALL clinical visits (contemporaneous CN/MCI/AD label) and extract per-visit embeddings in ONE
# consistent space for the visits baseline -> 1 year. NO checkpoint is saved (HPC space is tight):
# 03j exports embeddings.parquet + metrics.json directly.
#
# Single model x T2_long_multiclass x full_ft x seeds {0,1,2} = 3 runs (throttled %4).
# SLURM_ARRAY_TASK_ID IS the seed (0/1/2).
#
# Prereqs:
#   - Build verbose/longitudinal_allvisits locally (run 01i) and rsync it to RDS.
#   - git pull on the HPC so 03j_finetune_longitudinal_extract.py is present.
#   - mkdir -p ${OUT_DIR}/logs
# Submit (run from ${OUT_DIR} so any relative paths land there):
#   sbatch clinical_pipeline/03j_finetune_longitudinal_submit_csd3.sh
# Smoke seed 0 only:
#   sbatch --array=0 clinical_pipeline/03j_finetune_longitudinal_submit_csd3.sh
# Pull back ONLY: ${OUT_DIR}/<model>/T2_long_multiclass/seed_*/full_ft/{metrics.json,embeddings/}
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"
# Run EXPLICITLY in the 'clinical' env via `conda run` — independent of whatever env is active on the
# login node at submit time (no activate/inherit/stacking). Fail loudly if it's missing.
conda run -n clinical python -c "import numpy, torch, transformers" || {
    echo "ERROR: 'clinical' env missing/broken. conda env list:"; conda env list; exit 1; }

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# ALL-VISITS longitudinal splits (seed_N/{train,val,test}.csv with Generated_Text per visit,
# concurrent Label_visit_diag, visit_months), built by 01i. 03j appends seed_{N} itself.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/longitudinal_allvisits"

# Separate longitudinal output root (isolated from the baseline / m06 / m12 / T4 roots)
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_longitudinal"

MODEL_ID="thomas-sounack/BioClinical-ModernBERT-large"
MODEL_SLUG="${MODEL_ID##*/}"
SEED=${SLURM_ARRAY_TASK_ID}

METRICS="${OUT_DIR}/${MODEL_SLUG}/T2_long_multiclass/seed_${SEED}/full_ft/metrics.json"

echo "========================================================"
echo "  Longitudinal single-space T2 encoder (Task 1)"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}"
echo "  Model          : ${MODEL_ID}"
echo "  Seed           : ${SEED}"
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
