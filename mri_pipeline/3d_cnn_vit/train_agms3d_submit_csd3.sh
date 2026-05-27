#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=agms3d
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs/slurm_logs/agms3d_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs/slurm_logs/agms3d_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%4
# =============================================================================
# train_agms3d_submit_csd3.sh   --   AG-MS3D-CNN vanilla variant + rescue config
# =============================================================================
# SLURM array job: trains the AG-MS3D-CNN classifier on the 5 diagnostic tasks.
#
# This is the RESCUE sweep -- the prior vanilla+lr=1e-3+ls=0.1 L4 smoke
# collapsed (val_bacc 0.5, train_loss flat at ln 2). This sweep bundles the
# higher-LR + no-label-smoothing diagnostic the user approved, plus a batch
# size bump from 4 -> 8 since A100-80GB has ample VRAM (vs the shared L4):
#   --backbone vanilla       (nn.Conv3d instead of depthwise-separable)
#   --lr 3e-3                (bigger steps to escape the majority-class minimum)
#   --label_smoothing 0.0    (LS+weighted CE was flattening the loss landscape)
#   --batch_size 8           (A100-80GB headroom; trainer default is 4 for L4)
#   weighted CE kept on      (was already enabled in the prior run)
# Outputs land under AGMS3DCNN_vanilla/ -- sibling to the existing
# AGMS3DCNN_separable/ (formerly AGMS3DCNN/) tree, no overwriting.
#
# 5 tasks x 3 seeds = 15 combinations.
# Array index decoder: SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task  : T1_binary | T1b_binary | T1c_binary | T1d_binary | T2_multiclass
#   Seed  : 0 | 1 | 2
#
# CPU allocation per CSD3 ampere docs: slurm allocates 32 CPUs per GPU
# automatically; we request 8 (the dataloader sweet spot for 193x229x193
# single-channel volumes -- diminishing returns above 8). num_workers
# scales from $SLURM_CPUS_PER_TASK so any future bump is one number.
#
# Runs offline wandb -> sync to project 'agms3d_vanilla_rescue' after the sweep:
#   wandb sync <OUT_DIR>/AGMS3DCNN_vanilla/*/seed_*/wandb/offline-run-*
#
# SLURM .log/.err and per-run offline-wandb data land under the OUTPUT tree,
# never in the scripts/git directory.
#
# Pre-flight (run once before submitting):
#   1. mri conda env available.
#   2. CNN inputs at ${CNN_INPUTS_DIR} (produced by 00_prepare_CNN_submit_csd3.sh).
#   3. Post-exclusion matched labels CSV + clinical splits at the paths below.
#   4. mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
#
# Submit:  sbatch mri_pipeline/3d_cnn_vit/train_agms3d_submit_csd3.sh
# Smoke 1: sbatch --array=0 mri_pipeline/3d_cnn_vit/train_agms3d_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Offline wandb -- compute nodes have no internet; sync after the sweep.
export WANDB_MODE=offline

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/3d_cnn_vit"
CNN_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs"

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
N_TASKS=${#TASKS[@]}      # 5
N_SEEDS=${#SEEDS[@]}      # 3

# Decode SLURM_ARRAY_TASK_ID -> (task, seed)
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

RUN_DIR="${OUT_DIR}/AGMS3DCNN_vanilla/${TASK}/seed_${SEED}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"

# Offline wandb data for this run lands beside its checkpoints (output tree).
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  AG-MS3D-CNN -- array sweep"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS-1)))"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  CPUs/task   : $SLURM_CPUS_PER_TASK   (-> num_workers)"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : offline -> project agms3d_vanilla_rescue"
echo "  Started     : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -d "${CNN_INPUTS_DIR}" ]; then
    echo "[ERROR] CNN inputs dir not found: ${CNN_INPUTS_DIR}"
    echo "        Run 00_prepare_CNN_submit_csd3.sh first."
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

# -- Resume support: skip if metrics.json already exists ----------------------
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."
    exit 0
fi

python "${SCRIPT_DIR}/train_agms3d.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --backbone            vanilla \
    --lr                  3e-3 \
    --label_smoothing     0.0 \
    --batch_size          8 \
    --cnn_inputs_dir      "${CNN_INPUTS_DIR}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --wandb \
    --wandb_project       agms3d_vanilla_rescue \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
