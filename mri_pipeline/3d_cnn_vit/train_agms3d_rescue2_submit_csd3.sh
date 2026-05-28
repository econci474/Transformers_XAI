#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=agms3d_r2
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs_rescue2/slurm_logs/agms3d_r2_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs_rescue2/slurm_logs/agms3d_r2_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%4
# =============================================================================
# train_agms3d_rescue2_submit_csd3.sh   --   AG-MS3D-CNN rescue attempt #2
# =============================================================================
# RESCUE #1 (lr=3e-3, ls=0.0, batch=8, vanilla Conv3d, basic flips only)
# collapsed on 14/15 cells -- val_bacc stuck at chance, train_loss frozen at
# ln(K). See wandb project `agms3d_vanilla_rescue` for evidence. Diagnosis:
# from-scratch 3D CNN on ~580 train subjects, head too big (~130k params on
# top of a tiny encoder), too little augmentation, lr too aggressive.
#
# RESCUE #2 bundles three single-variable changes designed to break the
# constant-predictor fixed point without touching the architecture topology:
#
#   --lr 5e-4                  (vs rescue1's 3e-3; matches Spasov-CNN's safe LR)
#   --label_smoothing 0.1      (vs rescue1's 0.0; softens majority-class minimum)
#   --head slim                (256 -> 128 -> n_outputs; ~6x fewer head params,
#                               66k vs 395k -- reduces capacity for the
#                               constant-predictor lock-in)
#   --strong_aug               (adds 4 MONAI augmentations on top of axis
#                               flips: RandAffine, RandGaussianNoise,
#                               RandBiasField, RandAdjustContrast)
#
# Plus the rescue1-confirmed defaults that stay in place:
#   --backbone vanilla         (standard nn.Conv3d, not depthwise-separable)
#   --batch_size 8             (A100-80GB headroom)
#
# Output slug: AGMS3DCNN_vanilla_slim/  (the slim-head triggers the new slug
# in train_agms3d.py:454; vanilla+large stays at AGMS3DCNN_vanilla/ for
# rescue1 compatibility).
#
# Output root: agms3d_outputs_rescue2/  (separate tree -- does NOT collide
# with the rescue1 sweep at agms3d_outputs/).
#
# 5 tasks x 3 seeds = 15 combinations (array 0-14%4).
# Array decoder: SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task  : T1_binary | T1b_binary | T1c_binary | T1d_binary | T2_multiclass
#   Seed  : 0 | 1 | 2
#
# What to expect
# --------------
# - The "slim head + LS + lower LR" combination is the standard rescue
#   recipe for from-scratch 3D on tiny data; if it still collapses on >5
#   of 15 cells, AG-MS3D-vanilla has earned its "negative result" badge
#   for the thesis and the writeup should attribute it to no-pretraining
#   plus small ADNI sample (a thesis-relevant finding in itself --
#   contrast with the pretrained-model rows which all clear 0.5).
# - If T1c_binary alone escapes (matching rescue1 s0), that's a sign the
#   easier task is on the edge of trainable; not a rescue.
# - Any cell crossing val_bacc > 0.55 with finite train_loss decrease is
#   a real escape -- promote that recipe.
#
# W&B project: agms3d_vanilla_rescue2   (NEW project, not shared with
# agms3d_vanilla_rescue so the two attempts visually separate in the UI).
#
# Pre-flight (run once before submitting):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs_rescue2/slurm_logs
#
# Submit:  sbatch mri_pipeline/3d_cnn_vit/train_agms3d_rescue2_submit_csd3.sh
# Smoke :  sbatch --array=0 mri_pipeline/3d_cnn_vit/train_agms3d_rescue2_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Offline wandb -- compute nodes have no internet; sync after the sweep.
export WANDB_MODE=offline

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/3d_cnn_vit"
CNN_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs_rescue2"

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

# Output slug for the slim-head variant -- matches what train_agms3d.py emits.
MODEL_SLUG="AGMS3DCNN_vanilla_slim"
RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"

# Offline wandb data for this run lands beside its checkpoints (output tree).
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  AG-MS3D-CNN -- rescue #2"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS-1)))"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  CPUs/task   : $SLURM_CPUS_PER_TASK   (-> num_workers)"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : offline -> project agms3d_vanilla_rescue2"
echo "  Recipe      : vanilla Conv3d | head=slim | strong_aug | lr=5e-4 | ls=0.1 | bs=8"
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
    --head                slim \
    --strong_aug \
    --lr                  5e-4 \
    --label_smoothing     0.1 \
    --batch_size          8 \
    --cnn_inputs_dir      "${CNN_INPUTS_DIR}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --wandb \
    --wandb_project       agms3d_vanilla_rescue2 \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
