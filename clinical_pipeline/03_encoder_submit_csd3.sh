#!/bin/bash
#SBATCH --job-name=encoder_llm
#SBATCH --output=logs/encoder_%A_%a.log
#SBATCH --error=logs/encoder_%A_%a.err
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --array=0-95
# =============================================================================
# 03_encoder_submit_csd3.sh
# =============================================================================
# SLURM array job for CSD3 (A100 GPU).
# 4 models × 4 tasks × 3 seeds × 2 strategies = 96 combinations
# Each array task runs one combination.
#
# Estimated time per job on A100:
#   base models  : ~10 min (full FT, 5 epochs)
#   large models : ~25 min (full FT, 5 epochs)
#   → Set --time=01:30:00 as generous buffer for large + frozen (20 epochs)
#
# Setup (run once on login node before submitting):
#   module load cuda/12.1
#   bash setup_env.sh
#   export HF_TOKEN="hf_xxxx"
#   export HF_HOME="/rds/user/USERNAME/hpc-work/hf_cache"
#   bash download_models.sh
#
# Submit:
#   mkdir -p logs
#   sbatch 03_encoder_submit_csd3.sh
#
# Resume failed jobs only:
#   sbatch --array=$(python get_pending_ids.py) 03_encoder_submit_csd3.sh
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

# HF cache on RDS (avoids home directory quota)
# Adjust USERNAME to your CSD3 username
export HF_HOME="/rds/user/${USER}/hpc-work/hf_cache"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="/rds/user/${USER}/hpc-work/ADNI_BIDS_project/derivatives/clinical/verbose/baseline"
OUT_DIR="${SCRIPT_DIR}/outputs"

# ── Combination lookup ────────────────────────────────────────────────────────
MODEL_IDS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)
TASKS=("T1_binary" "T2_multiclass" "T3a_conv3y" "T3b_conv5y")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_MODELS=${#MODEL_IDS[@]}    # 4
N_TASKS=${#TASKS[@]}         # 4
N_SEEDS=${#SEEDS[@]}         # 3
N_STRAT=${#STRATEGIES[@]}    # 2

# Decode SLURM_ARRAY_TASK_ID → (model, task, seed, strategy)
ID=${SLURM_ARRAY_TASK_ID}

STRAT_IDX=$(( ID % N_STRAT ));          ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));           ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ));           ID=$(( ID / N_TASKS ))
MODEL_IDX=$(( ID % N_MODELS ))

MODEL_ID="${MODEL_IDS[$MODEL_IDX]}"
TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

FREEZE_FLAG=""
[[ "${STRATEGY}" == "frozen" ]] && FREEZE_FLAG="--freeze_backbone"

MODEL_SLUG="${MODEL_ID##*/}"
METRICS="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}/metrics.json"

echo "========================================================"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}"
echo "  Model          : ${MODEL_ID}"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "========================================================"

# Skip if already done
if [ -f "${METRICS}" ]; then
    echo "  metrics.json exists — skipping."
    exit 0
fi

python "${SCRIPT_DIR}/03_encoder_finetune.py" \
    --model_id  "${MODEL_ID}" \
    --task      "${TASK}" \
    --seed      "${SEED}" \
    --data_dir  "${DATA_DIR}" \
    --out_dir   "${OUT_DIR}" \
    --hf_cache  "${HF_HOME}" \
    ${FREEZE_FLAG}

echo "  Finished."
