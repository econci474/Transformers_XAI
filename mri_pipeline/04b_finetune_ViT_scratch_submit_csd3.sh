#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_scratch
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/slurm_logs/vit_scratch_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/slurm_logs/vit_scratch_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-11%4
# =============================================================================
# 04b_finetune_ViT_scratch_submit_csd3.sh  --  ViT-from-scratch baseline sweep
# =============================================================================
# SLURM array job: trains the ViT-B/3D FROM SCRATCH (random init, no MAE
# pretraining) -- the non-pretrained baseline for the thesis comparison.
#
# Kept fully SEPARATE from the debugging sweep: everything for this job --
# checkpoints, metrics.json, offline-wandb data, and the SLURM .log/.err --
# lives under  /home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/  so it never
# mixes with the full_ft/frozen runs in vit_outputs_debug/ (those are run by
# 04_finetune_ViT_submit_csd3.sh).
#
# 1 model x 4 tasks x 3 seeds x 1 strategy (scratch) = 12 combinations.
# Array index decoder: SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task : T1_binary | T1b_binary | T1c_binary | T2_multiclass
#   Seed : 0 | 1 | 2
#
# Runs offline wandb -> project 'vit_baseline'. wandb sync after the sweep:
#   wandb sync <OUT_DIR>/ViT_B_scratch/*/seed_*/scratch/wandb/offline-run-*
#
# Pre-flight (run once before submitting):
#   1. mri conda env.
#   2. POST-EXCLUSION matched labels CSV + clinical splits at the paths below.
#   3. ViT inputs at ${VIT_INPUTS_DIR} (from 03, --long all).
#   4. mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
#   No pretrained checkpoint is needed -- scratch = random init.
#
# Submit:  sbatch mri_pipeline/04b_finetune_ViT_scratch_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Offline wandb -- the compute nodes have no internet; sync after the sweep.
export WANDB_MODE=offline

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline"

# -- Session selection (must match what 03 produced) --------------------------
PY_SESSION_FLAGS="--long all"      # every available session (bl..mAll)

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T2_multiclass")
SEEDS=(0 1 2)
N_TASKS=${#TASKS[@]}    # 4
N_SEEDS=${#SEEDS[@]}    # 3

# Decode SLURM_ARRAY_TASK_ID -> (task, seed)
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="scratch"

RUN_DIR="${OUT_DIR}/ViT_B_scratch/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"

# Offline wandb data for this run lands beside its checkpoints (output tree).
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-from-scratch baseline sweep"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY} (random init, no MAE checkpoint)"
echo "  Sessions       : ${PY_SESSION_FLAGS}"
echo "  Output dir     : ${RUN_DIR}"
echo "  wandb          : offline -> project vit_baseline (WANDB_DIR=${WANDB_DIR})"
echo "  Started        : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs directory not found: ${VIT_INPUTS_DIR}"
    echo "        Run 03_prepare_ViT_submit_csd3.sh first."
    exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"
    echo "        Upload the post-exclusion extended master to HPC."
    exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"
    echo "        Upload the no_cdr_stratified_post_exclusion splits to HPC."
    exit 1
fi

# -- Resume support: skip if metrics.json already exists ----------------------
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."
    exit 0
fi

# No --pretrained_ckpt: scratch = random init (04_supervised_finetuning_ViT.py
# requires the checkpoint only for full_ft / frozen).
python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --grad_accum_steps    8 \
    --wandb \
    --wandb_project       vit_baseline \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "  Output         : ${RUN_DIR}/"
echo "============================================================"

exit ${EXIT_CODE}
