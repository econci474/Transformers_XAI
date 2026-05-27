#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_lr_sweep
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_lr_sweep/slurm_logs/bmvp_lr_sweep_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_lr_sweep/slurm_logs/bmvp_lr_sweep_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-44%4
# =============================================================================
# 04h_finetune_BrainMVP_lr_sweep_submit_csd3.sh  --  3-point LR sweep, full_ft
# =============================================================================
# QUEUE THIS *AFTER* the current sweeps finish on HPC (AG-MS3D rescue,
# BrainDINO LoRA, BrainDINO full_ft, BrainDINO frozen+aug, ViT-MAE hi_lr).
# It is a follow-up to the existing brainmvp_debug sweep (commit 3ca882a),
# not a replacement.
#
# What this sweep tests
# ---------------------
# The original BrainMVP UniFormer fine-tuning at lr=3e-4 NaN'd every full_ft
# run ~epoch 27. The post-NaN fix dropped LR to 5e-5. That was the largest
# stable value we hand-picked; we never re-swept around it to confirm 5e-5
# is the optimum vs a conservative ceiling.
#
# This sweep brackets 5e-5 with one smaller (3e-5) and one larger (7e-5)
# value, on the strategy that's already strongest (full_ft) at the best-
# performing augmentation PER TASK from the prior brainmvp_debug run:
#
#   T1_binary      -> plus_original   (current best 0.580)
#   T1b_binary     -> stochastic      (current best 0.720)
#   T1c_binary     -> stochastic      (current best 0.792)
#   T1d_binary     -> plus_original   (current best 0.712)
#   T2_multiclass  -> plus_original   (current best 0.489)
#
# Grid: 3 LRs x 5 tasks x 3 seeds = 45 cells. Strategy fixed to full_ft
# (the frozen+head variant consistently underperformed full_ft by 0.05-0.15
# bacc; no value in re-sweeping that).
#
# Array decoder: LR_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   LR     : 3e-5 | 5e-5 | 7e-5
#   Seed   : 0 | 1 | 2
#   Task   : T1 | T1b | T1c | T1d | T2  (augment derived from task above)
#
# Output layout (LR-scoped subroot so the 3 LR variants don't collide):
#   brainmvp_lr_sweep/lr_${LR}/aug_${AUGMENT}/BrainMVP_uniformer/<task>/seed_<n>/full_ft/
#
# What to expect
# --------------
# - 7e-5 might NaN on some cells -- the abort-on-non-finite guard from the
#   post-NaN fix is already in the trainer, so failures cost minutes, not
#   the full 12 h wallclock. metrics.json simply won't be written; check
#   the .err log for the abort line.
# - 3e-5 should be stable but may underfit the head, especially on T1d/T2.
# - If 7e-5 works AND matches/beats 5e-5, that's evidence the encoder has
#   headroom to train faster. Worth a small follow-up at 8-9e-5.
# - If neither 3e-5 nor 7e-5 beats 5e-5, BrainMVP is at its full_ft ceiling
#   for this recipe and further LR sweeping is over.
#
# W&B project: brainmvp_lr_sweep  (NEW project, not shared with brainmvp_debug)
#
# Pre-flight (run once before submitting):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_lr_sweep/slurm_logs
#
# Submit:  sbatch mri_pipeline/brain_mvp/04h_finetune_BrainMVP_lr_sweep_submit_csd3.sh
# Smoke :  sbatch --array=0 mri_pipeline/brain_mvp/04h_finetune_BrainMVP_lr_sweep_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Run-mode knobs -----------------------------------------------------------
AUG_COPIES="${AUG_COPIES:-1}"
USE_WANDB="${USE_WANDB:-1}"
WANDB_PROJECT="${WANDB_PROJECT:-brainmvp_lr_sweep}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"

OUT_ROOT="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_lr_sweep"

# -- Session: all longitudinal ------------------------------------------------
PY_SESSION_FLAGS="--long all"

# -- Sweep grid: 3 LRs x 5 tasks x 3 seeds = 45 -------------------------------
LRS=("3e-5" "5e-5" "7e-5")
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
STRATEGY="full_ft"   # fixed; frozen variant not under test

N_LRS=${#LRS[@]}
N_TASKS=${#TASKS[@]}
N_SEEDS=${#SEEDS[@]}

ID=${SLURM_ARRAY_TASK_ID}
LR_IDX=$(( ID % N_LRS ));      ID=$(( ID / N_LRS ))
SEED_IDX=$(( ID % N_SEEDS ));  ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

LR="${LRS[$LR_IDX]}"
TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

# -- Per-task best augment (from prior brainmvp_debug results) ----------------
case "${TASK}" in
    T1_binary)      AUGMENT="plus_original" ;;
    T1b_binary)     AUGMENT="stochastic"    ;;
    T1c_binary)     AUGMENT="stochastic"    ;;
    T1d_binary)     AUGMENT="plus_original" ;;
    T2_multiclass)  AUGMENT="plus_original" ;;
    *) echo "[ERROR] unknown TASK=${TASK}"; exit 1 ;;
esac

# -- Output paths (LR-scoped subroot so the 3 LRs don't collide) --------------
OUT_DIR="${OUT_ROOT}/lr_${LR}/aug_${AUGMENT}"
RUN_DIR="${OUT_DIR}/BrainMVP_uniformer/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"
CKPT="${RUN_DIR}/last_checkpoint.pt"

mkdir -p "${RUN_DIR}"
mkdir -p "${OUT_ROOT}/slurm_logs"
export WANDB_DIR="${RUN_DIR}"

# -- wandb flags --------------------------------------------------------------
WANDB_FLAGS=""
if [ "${USE_WANDB}" = "1" ]; then
    WANDB_FLAGS="--wandb --wandb_project ${WANDB_PROJECT}"
fi

echo "============================================================"
echo "  BrainMVP UniFormer -- 3-point LR sweep (full_ft, best-aug-per-task)"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_LRS*N_SEEDS*N_TASKS-1)))"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  Strategy    : ${STRATEGY}"
echo "  LR          : ${LR}"
echo "  Augment     : ${AUGMENT}  (aug_copies=${AUG_COPIES})"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started     : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${BRAINMVP_INPUTS_DIR}" ]; then
    echo "[ERROR] BrainMVP inputs dir not found: ${BRAINMVP_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"; exit 1
fi

# -- Resume support ------------------------------------------------------------
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- will AUTO-RESUME."
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_BrainMVP.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --lr                  "${LR}" \
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
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
if [ ! -f "${METRICS}" ] && [ -f "${CKPT}" ]; then
    echo "  NOTE: no metrics.json yet (likely hit 12 h cap). Resubmit to AUTO-RESUME."
fi
if [ ! -f "${METRICS}" ] && [ ! -f "${CKPT}" ]; then
    echo "  NOTE: no metrics.json AND no checkpoint -- likely NaN-aborted early."
    echo "        Check the .err log for the abort line. (LR=${LR} may be too high.)"
fi
echo "============================================================"

exit ${EXIT_CODE}
