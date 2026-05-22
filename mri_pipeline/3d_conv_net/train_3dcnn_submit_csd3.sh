#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=cnn3d
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/cnn3d_outputs/slurm_logs/cnn3d_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/cnn3d_outputs/slurm_logs/cnn3d_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-5%4
# =============================================================================
# train_3dcnn_submit_csd3.sh   --   Spasov 3D CNN baseline sweep
# =============================================================================
# SLURM array job: trains the 3D CNN baselines on CN-vs-AD (T1c_binary).
#
# 2 models x 3 seeds = 6 combinations.
# Array index decoder: SEED_IDX -> MODEL_IDX (low-to-high stride).
#
#   Model : vanilla (MRI3DCNN) | separable (MRI3DSeparableCNN)
#   Seed  : 0 | 1 | 2
#
# Runs offline wandb -> sync to project 'cnn3d_baseline' after the sweep:
#   wandb sync <OUT_DIR>/Spasov3DCNN_*/*/seed_*/wandb/offline-run-*
#
# SLURM .log/.err and per-run offline-wandb data land under the OUTPUT tree,
# never in the scripts/git directory.
#
# Pre-flight (run once before submitting):
#   1. mri conda env available; pretrained data uploaded.
#   2. CNN inputs at ${CNN_INPUTS_DIR} (produced by 00_prepare_CNN_submit_csd3.sh).
#   3. Post-exclusion matched labels CSV + clinical splits at the paths below.
#   4. mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/cnn3d_outputs/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
#
# Submit:  sbatch mri_pipeline/3d_conv_net/train_3dcnn_submit_csd3.sh
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
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/3d_conv_net"
CNN_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/cnn3d_outputs"

# -- Combination lookup -------------------------------------------------------
MODELS=("vanilla" "separable")
SEEDS=(0 1 2)
N_MODELS=${#MODELS[@]}   # 2
N_SEEDS=${#SEEDS[@]}     # 3

# Decode SLURM_ARRAY_TASK_ID -> (model, seed)
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
MODEL_IDX=$(( ID % N_MODELS ))

MODEL="${MODELS[$MODEL_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

RUN_DIR="${OUT_DIR}/Spasov3DCNN_${MODEL}/T1c_binary/seed_${SEED}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"

# Offline wandb data for this run lands beside its checkpoints (output tree).
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  Spasov 3D CNN -- array sweep"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_MODELS*N_SEEDS-1)))"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Model       : ${MODEL}"
echo "  Seed        : ${SEED}"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : offline -> project cnn3d_baseline"
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

python "${SCRIPT_DIR}/train_3dcnn.py" \
    --model               "${MODEL}" \
    --seed                "${SEED}" \
    --cnn_inputs_dir      "${CNN_INPUTS_DIR}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --wandb \
    --wandb_project       cnn3d_baseline \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
