#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_nocdr
#SBATCH --output=logs/enc_nocdr_%A_%a.log
#SBATCH --error=logs/enc_nocdr_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --array=0-95
# =============================================================================
# 03b_encoder_stratified_no_cdr_submit_csd3.sh
# =============================================================================
# Encoder-only LLM fine-tuning on NO-CDR STRATIFIED verbose baseline splits.
#
# Identical model/task/seed/strategy grid as 03_encoder_submit_csd3.sh, but:
#   - DATA_DIR  -> no_cdr_stratified verbose baseline (CDR columns and text
#                  stripped to prevent diagnostic leakage)
#   - OUT_DIR   -> separate output directory for clean comparison
#
# 4 models × 4 tasks × 3 seeds × 2 strategies = 96 combinations.
#
# Submit:
#   mkdir -p logs
#   sbatch clinical_pipeline/03b_encoder_stratified_no_cdr_submit_csd3.sh
#
# Resume failed jobs only:
#   sbatch --array=$(python get_pending_ids.py) clinical_pipeline/03b_encoder_stratified_no_cdr_submit_csd3.sh
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

# HF cache on RDS (avoids home directory quota)
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# KEY CHANGE: no-CDR stratified verbose baseline splits
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline_stratified_no_cdr_verbose"

# KEY CHANGE: separate output directory
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_stratified"

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
echo "  Encoder-only LLM — No-CDR Stratified"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}"
echo "  Model          : ${MODEL_ID}"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Data dir       : ${DATA_DIR}"
echo "  Output dir     : ${OUT_DIR}"
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
