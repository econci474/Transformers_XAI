#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_ft
#SBATCH --output=logs/bmvp_ft_%A_%a.log
#SBATCH --error=logs/bmvp_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-17%4
# =============================================================================
# 04g_finetune_BrainMVP_submit_csd3.sh
# =============================================================================
# Fine-tunes BrainMVP UniFormer-Small on ADNI MRI classification tasks.
# Same experimental setup as the ViT pipeline (04g_finetune_ViT_combined).
#
# 3 tasks × 3 seeds × 2 strategies = 18 combinations (array 0-17).
#
# Output layout: brain_mvp_uniformer/aug_${AUGMENT}/BrainMVP_uniformer/...
#
# Usage:
#   sbatch --export=ALL,AUGMENT=none          mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=stochastic     mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=plus_original  mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Run-mode knobs -----------------------------------------------------------
AUGMENT="${AUGMENT:-stochastic}"
AUG_COPIES="${AUG_COPIES:-1}"
USE_WANDB="${USE_WANDB:-1}"
WANDB_PROJECT="${WANDB_PROJECT:-brainmvp_finetune}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"

# VISCODE2-aligned master from 03c
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/viscode_2_aligned/master_mri_clinical_matched_viscode2.csv"

# No-CDR stratified baseline splits
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline_stratified_no_cdr"

# Output root: brain_mvp_uniformer/aug_${AUGMENT}
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brain_mvp_uniformer/aug_${AUGMENT}"

# -- Session: all longitudinal ------------------------------------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup: 3 tasks, 3 seeds, 2 strategies -----------------------
TASKS=("T1_binary" "T1b_binary" "T2_multiclass")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_TASKS=${#TASKS[@]}
N_SEEDS=${#SEEDS[@]}
N_STRAT=${#STRATEGIES[@]}

ID=${SLURM_ARRAY_TASK_ID}
STRAT_IDX=$(( ID % N_STRAT ));   ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

RUN_DIR="${OUT_DIR}/BrainMVP_uniformer/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"
CKPT="${RUN_DIR}/last_checkpoint.pt"

mkdir -p logs
mkdir -p "${OUT_DIR}"

# -- wandb flags --------------------------------------------------------------
WANDB_FLAGS=""
if [ "${USE_WANDB}" = "1" ]; then
    WANDB_FLAGS="--wandb --wandb_project ${WANDB_PROJECT}"
fi

echo "============================================================"
echo "  BrainMVP UniFormer fine-tuning"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Augment        : ${AUGMENT}  (aug_copies=${AUG_COPIES})"
echo "  Output dir     : ${RUN_DIR}"
echo "  Started        : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"
    exit 1
fi

# -- Resume support ------------------------------------------------------------
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."
    exit 0
fi
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- will AUTO-RESUME."
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_BrainMVP.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --augment             "${AUGMENT}" \
    --aug_copies          "${AUG_COPIES}" \
    ${WANDB_FLAGS} \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --brainmvp_inputs_dir "${BRAINMVP_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         2

EXIT_CODE=$?

echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "============================================================"

exit ${EXIT_CODE}
