#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_10y
#SBATCH --output=logs/enc_10y_%A_%a.log
#SBATCH --error=logs/enc_10y_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --array=0-23
# =============================================================================
# 03c_encoder_stratified_no_cdr_10y_submit_csd3.sh
# =============================================================================
# Runs the T3c_conv10y task ONLY on the no-CDR stratified verbose baseline.
# Output saves alongside the existing T1/T2/T3a/T3b results in the same
# directory structure (encoder_outputs_no_cdr_stratified/).
#
# 4 models × 1 task × 3 seeds × 2 strategies = 24 combinations.
#
# Submit:
#   sbatch clinical_pipeline/03c_encoder_stratified_no_cdr_10y_submit_csd3.sh
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# SAME data and output dirs as 03b — results slot into existing structure
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline_stratified_no_cdr_verbose"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_stratified"

# ── Combination lookup: T3c_conv10y only ──────────────────────────────────────
MODEL_IDS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)
TASK="T3c_conv10y"
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_MODELS=${#MODEL_IDS[@]}    # 4
N_SEEDS=${#SEEDS[@]}         # 3
N_STRAT=${#STRATEGIES[@]}    # 2

# Decode SLURM_ARRAY_TASK_ID → (model, seed, strategy)
ID=${SLURM_ARRAY_TASK_ID}

STRAT_IDX=$(( ID % N_STRAT ));          ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));           ID=$(( ID / N_SEEDS ))
MODEL_IDX=$(( ID % N_MODELS ))

MODEL_ID="${MODEL_IDS[$MODEL_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

FREEZE_FLAG=""
[[ "${STRATEGY}" == "frozen" ]] && FREEZE_FLAG="--freeze_backbone"

MODEL_SLUG="${MODEL_ID##*/}"
METRICS="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}/metrics.json"

echo "========================================================"
echo "  Encoder-only LLM — T3c_conv10y (No-CDR Stratified)"
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
