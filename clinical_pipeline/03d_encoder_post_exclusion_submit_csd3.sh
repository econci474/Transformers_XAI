#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_postexcl
# Logs land UNDER the output directory (not the script dir). Slurm does NOT expand
# shell variables here, so this absolute path must match OUT_DIR/logs below.
# Create it before submitting:  mkdir -p <OUT_DIR>/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion/logs/enc_postexcl_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion/logs/enc_postexcl_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-239%4
# =============================================================================
# 03d_encoder_post_exclusion_submit_csd3.sh
# =============================================================================
# Encoder-only LLM fine-tuning on the NO-CDR STRATIFIED *POST-EXCLUSION* verbose
# baseline splits. Fresh runs for the post-exclusion cohort.
#
#   - DATA_DIR -> post-exclusion verbose baseline (UPLOAD this to RDS first; see below)
#   - OUT_DIR  -> separate post-exclusion output directory
#   - 10 tasks (adds T1b CN+MCI/AD, T1c sCN/prog, T1d pMCI/sMCI, T1e sCN/pCN, T3d conv7y)
#   - W&B logging in OFFLINE mode (sync on HPC afterwards)
#
# 4 models x 10 tasks x 3 seeds x 2 strategies = 240 combinations (throttled %4).
#
# Submit (logs go to OUT_DIR/logs):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion/logs
#   sbatch clinical_pipeline/03d_encoder_post_exclusion_submit_csd3.sh
#
# After all runs finish, sync the offline W&B runs:
#   wandb sync --sync-all "${OUT_DIR}/wandb_offline"
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

# HF cache on RDS (avoids home directory quota)
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# Post-exclusion verbose baseline splits (seed_N/{train,val,test}.csv with Generated_Text).
# 03_encoder_finetune.py appends seed_{N} itself, so DATA_DIR must be the dir CONTAINING
# the seed_* folders. Local source: D:\ADNI_BIDS_project\derivatives\clinical\
#   no_cdr_stratified_post_exclusion\verbose\baseline  ->  this RDS path.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline"

# Separate post-exclusion output directory
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion"

# W&B offline — collect all runs under one dir for a single `wandb sync --sync-all`
export WANDB_MODE="offline"
export WANDB_DIR="${OUT_DIR}/wandb_offline"
mkdir -p "${WANDB_DIR}"
WANDB_PROJECT="clinical-encoder-post-exclusion"

# ── Combination lookup ────────────────────────────────────────────────────────
MODEL_IDS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)
TASKS=(
    "T1_binary"
    "T1b_cnmci_ad"
    "T2_multiclass"
    "T1c_scn_prog"
    "T1d_pmci_smci"
    "T1e_scn_pcn"
    "T3a_conv3y"
    "T3b_conv5y"
    "T3d_conv7y"
    "T3c_conv10y"
)
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_MODELS=${#MODEL_IDS[@]}    # 4
N_TASKS=${#TASKS[@]}         # 10
N_SEEDS=${#SEEDS[@]}         # 3
N_STRAT=${#STRATEGIES[@]}    # 2

# Decode SLURM_ARRAY_TASK_ID -> (model, task, seed, strategy)
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
echo "  Encoder-only LLM — No-CDR Stratified POST-EXCLUSION"
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
