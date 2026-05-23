#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-23%4
# =============================================================================
# 04g_finetune_BrainMVP_submit_csd3.sh   --   BrainMVP "brainmvp_debug" re-run
# =============================================================================
# Fine-tunes BrainMVP UniFormer-Small on ADNI MRI classification tasks.
#
# 4 tasks x 3 seeds x 2 strategies = 24 combinations (array 0-23).
# Array index decoder: STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task     : T1_binary | T1b_binary | T1c_binary | T2_multiclass
#   Seed     : 0 | 1 | 2
#   Strategy : full_ft | frozen
#
# Post-fix re-run: full_ft LR lowered to 5e-5 (the paper's 3e-4 NaN'd every
# full_ft run ~epoch 27), abort-on-NaN, balanced-accuracy checkpoint selection,
# and the post-exclusion data/splits (so it is comparable to the ViT sweep).
# Everything for this sweep -- checkpoints, metrics.json, offline-wandb, and the
# SLURM .log/.err -- lives under /home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/,
# separate from the old broken brain_mvp_uniformer/ tree.
#
# Output layout: brainmvp_debug/aug_${AUGMENT}/BrainMVP_uniformer/<task>/seed_<n>/<strategy>/
#
# Pre-flight (run once before submitting):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs
#   (the #SBATCH --output dir must exist BEFORE submitting).
#
# Usage (one submission per augmentation mode):
#   sbatch --export=ALL,AUGMENT=none          mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=stochastic    mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=plus_original mri_pipeline/brain_mvp/04g_finetune_BrainMVP_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
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
WANDB_PROJECT="${WANDB_PROJECT:-brainmvp_debug}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"

# Post-exclusion extended master + splits (same as the ViT debugging sweep).
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"

# Output root -- fresh tree, separate from the old broken brain_mvp_uniformer/.
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/aug_${AUGMENT}"

# -- Session: all longitudinal ------------------------------------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup: 4 tasks, 3 seeds, 2 strategies -----------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T2_multiclass")
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

mkdir -p "${RUN_DIR}"
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs"

# Offline wandb data for this run lands beside its checkpoints (output tree).
export WANDB_DIR="${RUN_DIR}"

# -- wandb flags --------------------------------------------------------------
WANDB_FLAGS=""
if [ "${USE_WANDB}" = "1" ]; then
    WANDB_FLAGS="--wandb --wandb_project ${WANDB_PROJECT}"
fi

echo "============================================================"
echo "  BrainMVP UniFormer fine-tuning -- brainmvp_debug re-run"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Augment        : ${AUGMENT}  (aug_copies=${AUG_COPIES})"
echo "  Output dir     : ${RUN_DIR}"
echo "  wandb          : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started        : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"
    exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"
    exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"
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
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "  Output         : ${RUN_DIR}/"
echo "============================================================"

exit ${EXIT_CODE}
