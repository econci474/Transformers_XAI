#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_T4
# Logs land UNDER the T4 output directory. Slurm does NOT expand shell vars here, so this
# absolute path must match OUT_DIR/logs below. Create it before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_T4/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_T4/logs/enc_T4_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_T4/logs/enc_T4_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-23%4
# =============================================================================
# 03f_encoder_T4_submit_csd3.sh
# =============================================================================
# Encoder-only LLM fine-tuning for the SINGLE task T4_conv_horizon (3-class
# AD-conversion horizon: <3y / 3-7y / >=7y, converters pMCI + pCN_to_AD),
# matching the MRI T4 task. Kept SEPARATE from the 240-run 03d grid:
#   - its OWN split tree   -> verbose/baseline_T4   (built by 01f from baseline_T4)
#   - its OWN output root   -> encoder_outputs_no_cdr_post_exclusion_T4
#   - its OWN wandb project -> clinical-encoder-T4
#
# 4 models x 1 task x 3 seeds x 2 strategies = 24 combinations (throttled %4).
#
# Prereqs:
#   - Build + upload verbose/baseline_T4 (run 01f locally, rsync to RDS).
#   - mkdir -p ${OUT_DIR}/logs
# Submit (run from ${OUT_DIR} so any relative paths land there):
#   sbatch clinical_pipeline/03f_encoder_T4_submit_csd3.sh
# Smoke one combo first:
#   sbatch --array=0 clinical_pipeline/03f_encoder_T4_submit_csd3.sh
#
# After runs finish, sync the offline W&B runs + build the T4 embeddings manifest:
#   wandb sync --sync-all "${OUT_DIR}/wandb_offline"
#   python clinical_pipeline/03e_build_embeddings_manifest.py \
#       --out_dir "${OUT_DIR}" --data_dir "${DATA_DIR}"
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

# HF cache on RDS (avoids home-dir quota)
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# T4 VERBOSE splits (seed_N/{train,val,test}.csv with Generated_Text + Label_T4),
# built by 01f from tabular/baseline_T4. 03_encoder_finetune.py appends seed_{N} itself,
# so DATA_DIR must be the dir CONTAINING the seed_* folders.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline_T4"

# Separate T4 output root (isolated from the 240-run sweep)
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_T4"

# W&B offline — collect under one dir for a single `wandb sync --sync-all`
export WANDB_MODE="offline"
export WANDB_DIR="${OUT_DIR}/wandb_offline"
mkdir -p "${WANDB_DIR}"
WANDB_PROJECT="clinical-encoder-T4"

# ── Combination lookup (single task) ──────────────────────────────────────────
MODEL_IDS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)
TASK="T4_conv_horizon"
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_MODELS=${#MODEL_IDS[@]}    # 4
N_SEEDS=${#SEEDS[@]}         # 3
N_STRAT=${#STRATEGIES[@]}    # 2

# Decode SLURM_ARRAY_TASK_ID -> (model, seed, strategy)   (no task index — single task)
ID=${SLURM_ARRAY_TASK_ID}
STRAT_IDX=$(( ID % N_STRAT ));   ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
MODEL_IDX=$(( ID % N_MODELS ))

MODEL_ID="${MODEL_IDS[$MODEL_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

FREEZE_FLAG=""
[[ "${STRATEGY}" == "frozen" ]] && FREEZE_FLAG="--freeze_backbone"

MODEL_SLUG="${MODEL_ID##*/}"
METRICS="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}/metrics.json"

echo "========================================================"
echo "  Encoder-only LLM — T4_conv_horizon (3-class, converters)"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}"
echo "  Model          : ${MODEL_ID}"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Data dir       : ${DATA_DIR}"
echo "  Output dir     : ${OUT_DIR}"
echo "  W&B (offline)  : ${WANDB_PROJECT}  (dir: ${WANDB_DIR})"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "========================================================"

# Skip if already done
if [ -f "${METRICS}" ]; then
    echo "  metrics.json exists — skipping."
    exit 0
fi

python "${SCRIPT_DIR}/03_encoder_finetune.py" \
    --model_id       "${MODEL_ID}" \
    --task           "${TASK}" \
    --seed           "${SEED}" \
    --data_dir       "${DATA_DIR}" \
    --out_dir        "${OUT_DIR}" \
    --hf_cache       "${HF_HOME}" \
    --wandb \
    --wandb_project  "${WANDB_PROJECT}" \
    ${FREEZE_FLAG}

echo "  Finished."
