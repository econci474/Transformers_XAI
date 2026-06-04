#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_T3_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_T3_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_T3_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-8%4
# =============================================================================
# 04j_finetune_BrainMVP_T3abc_full_ft_submit_csd3.sh
# =============================================================================
# BrainMVP UniFormer FULL fine-tuning on the 3 binary conversion windows
# T3a/T3b/T3c (conversion to AD within 3/5/7y). T3d/10y is omitted
# (degenerate). Counterpart to the cached-head T3 sweep + the T4 full_ft
# (04i); use this to test whether adapting the encoder beats the frozen
# cached-head linear probe on T3 (the cached-head T3 results sit near chance).
#
# Array covers 3 tasks x 3 seeds = 9 cells, submitted ONCE per augment mode:
#   id = TASK_IDX + 3*SEED_IDX  (TASK fastest-varying)
#
# T3 uses the STANDARD baseline split tree (NOT baseline_T4) and
# session_policy=baseline_anchored with label_anchor_max_months=12 (the
# abundant m12 scans, handled internally by the trainer; --long all is fine).
#
# Submit per augment:
#   sbatch --export=ALL,AUGMENT=plus_original mri_pipeline/brain_mvp/04j_finetune_BrainMVP_T3abc_full_ft_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=stochastic    mri_pipeline/brain_mvp/04j_finetune_BrainMVP_T3abc_full_ft_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=none          mri_pipeline/brain_mvp/04j_finetune_BrainMVP_T3abc_full_ft_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

AUGMENT="${AUGMENT:-plus_original}"
AUG_COPIES="${AUG_COPIES:-1}"
USE_WANDB="${USE_WANDB:-1}"
WANDB_PROJECT="${WANDB_PROJECT:-brainmvp_T3_ft}"
export WANDB_MODE="${WANDB_MODE:-offline}"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"

MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"

OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/aug_${AUGMENT}"

# T3 = baseline_anchored (Label_Ny broadcast on every visit row, capped to
# <=12 months internally); --long all is harmless.
PY_SESSION_FLAGS="--long all"

TASKS=("T3a_conv3y" "T3b_conv5y" "T3c_conv7y")
SEEDS=(0 1 2)
STRATEGY="full_ft"
N_TASKS=${#TASKS[@]}

ID=${SLURM_ARRAY_TASK_ID}
TASK_IDX=$(( ID % N_TASKS ));  ID=$(( ID / N_TASKS ))
SEED_IDX=$(( ID % ${#SEEDS[@]} ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

RUN_DIR="${OUT_DIR}/BrainMVP_uniformer/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${RUN_DIR}"
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs"
export WANDB_DIR="${RUN_DIR}"

WANDB_FLAGS=""
if [ "${USE_WANDB}" = "1" ]; then
    WANDB_FLAGS="--wandb --wandb_project ${WANDB_PROJECT}"
fi

echo "============================================================"
echo "  BrainMVP full_ft -- T3 binary conversion windows"
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Array ID  : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*${#SEEDS[@]}-1)))"
echo "  GPU       : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Task      : ${TASK}   Seed: ${SEED}   Aug: ${AUGMENT}"
echo "  Output    : ${RUN_DIR}"
echo "  wandb     : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started   : $(date)"
echo "============================================================"

if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Splits dir not found: ${DATA_DIR}/seed_${SEED}"; exit 1
fi
if [ ! -d "${BRAINMVP_INPUTS_DIR}" ]; then
    echo "[ERROR] BrainMVP inputs dir not found: ${BRAINMVP_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] BrainMVP pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_BrainMVP.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --augment             "${AUGMENT}" \
    --aug_copies          "${AUG_COPIES}" \
    ${PY_SESSION_FLAGS} \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --brainmvp_inputs_dir "${BRAINMVP_INPUTS_DIR}" \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         "${SLURM_CPUS_PER_TASK}" \
    ${WANDB_FLAGS}

EXIT_CODE=$?
echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
