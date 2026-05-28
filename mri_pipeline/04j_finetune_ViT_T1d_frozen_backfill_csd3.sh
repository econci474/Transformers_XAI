#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_T1d_frzn
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_T1d_frzn_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_T1d_frzn_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%3
# =============================================================================
# 04i_finetune_ViT_T1d_frozen_backfill_csd3.sh
# =============================================================================
# One-off backfill for the T1d_binary gap in the ViT-MAE75 / frozen / random
# sweep. The original 04_finetune_ViT_submit_csd3.sh sweep used
# TASKS=(T1_binary T1b_binary T1c_binary T2_multiclass) -- it never ran T1d.
# This script runs the missing 3 cells (T1d × 3 seeds × frozen × random)
# and writes to the canonical vit_outputs_debug tree so the rows slot into
# the cross-model aggregator alongside the existing 4-task frozen runs.
#
# Note: LLRD is a no-op for frozen training (only the head sees gradient),
# so the choice between --llrd_gamma=0.7 and 1.0 is mathematically irrelevant
# here. Using 0.7 to match the original sweep's recorded llrd_gamma for
# config consistency in metrics.json.
#
# Array decoder: ID = seed (3 cells, all T1d_binary + frozen + random).
#
#   ID 0 -> seed_0
#   ID 1 -> seed_1
#   ID 2 -> seed_2
#
# Output: vit_outputs_debug/ViT_B_mae75/T1d_binary/seed_<n>/frozen/
# W&B project: vit_debugging (matches the original sweep; the existing
# frozen rows are already in this project)
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs
#
# Submit:  sbatch mri_pipeline/04i_finetune_ViT_T1d_frozen_backfill_csd3.sh
# Smoke :  sbatch --array=0 mri_pipeline/04i_finetune_ViT_T1d_frozen_backfill_csd3.sh
# Already-completed runs (metrics.json present) auto-skip.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

export WANDB_MODE="${WANDB_MODE:-offline}"
WANDB_PROJECT="${WANDB_PROJECT:-vit_debugging}"

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug"

PY_SESSION_FLAGS="--long all"

# -- Fixed for this one-off backfill ------------------------------------------
TASK="T1d_binary"
STRATEGY="frozen"
AUGMENT="random"
LLRD_GAMMA="0.7"            # no-op for frozen, kept for config consistency

# -- Combination lookup (just seeds) -----------------------------------------
SEEDS=(0 1 2)
SEED="${SEEDS[$SLURM_ARRAY_TASK_ID]}"

MODEL_SLUG="ViT_B_mae75"
RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"
CKPT="${RUN_DIR}/last_checkpoint.pt"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-MAE75 frozen T1d backfill"
echo "  Job ID      : $SLURM_JOB_ID  Array: $SLURM_ARRAY_TASK_ID"
echo "  Task        : ${TASK}    Seed: ${SEED}    Strategy: ${STRATEGY}"
echo "  Recipe      : aug=${AUGMENT}  llrd_gamma=${LLRD_GAMMA} (no-op for frozen)"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : ${WANDB_PROJECT}"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs dir not found: ${VIT_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"; exit 1
fi

if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- python will AUTO-RESUME."
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --llrd_gamma          "${LLRD_GAMMA}" \
    --augment             "${AUGMENT}" \
    ${PY_SESSION_FLAGS} \
    --wandb \
    --wandb_project       "${WANDB_PROJECT}" \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
