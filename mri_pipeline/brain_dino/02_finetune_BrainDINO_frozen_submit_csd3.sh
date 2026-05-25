#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bdino_frozen
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_frozen_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_frozen_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%4
# =============================================================================
# 02_finetune_BrainDINO_frozen_submit_csd3.sh  --  encoder frozen, head only
# =============================================================================
# Fine-tunes BrainDINO ViT-B/16 with the encoder frozen -- only the
# classification head trains. Cheapest baseline for direct comparison with
# the BrainMVP frozen rows. Augmentation is fixed to 'stochastic' (skipping
# none/plus_original this round per the plan -- BrainMVP's frozen rows
# showed minimal sensitivity to the augment axis).
#
# 5 tasks x 3 seeds x 1 augment = 15 combinations (array 0-14).
# Array decoder: SEED_IDX -> TASK_IDX (low-to-high stride).
#
# Output layout: braindino_outputs/aug_stochastic/BrainDINO_vitb16_frozen/<task>/seed_<n>/
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_dino"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/brain_dino_model.pth"
BRAINDINO_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs"

PY_SESSION_FLAGS="--long all"
STRATEGY="frozen"
AUGMENT="stochastic"
WANDB_PROJECT="${WANDB_PROJECT:-braindino_frozen}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Combination lookup: 4 tasks x 3 seeds = 12 -------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)

N_TASKS=${#TASKS[@]}
N_SEEDS=${#SEEDS[@]}

ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));   ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

RUN_OUT_DIR="${OUT_DIR}/aug_${AUGMENT}/BrainDINO_vitb16_${STRATEGY}/${TASK}/seed_${SEED}"
METRICS="${RUN_OUT_DIR}/metrics.json"
CKPT="${RUN_OUT_DIR}/last_checkpoint.pt"

mkdir -p "${RUN_OUT_DIR}"
mkdir -p "${OUT_DIR}/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  BrainDINO ViT-B/16 frozen (head only)"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS-1)))"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  Strategy    : ${STRATEGY}"
echo "  Augment     : ${AUGMENT}"
echo "  Output dir  : ${RUN_OUT_DIR}"
echo "  wandb       : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${BRAINDINO_INPUTS_DIR}" ]; then
    echo "[ERROR] BrainDINO inputs dir not found: ${BRAINDINO_INPUTS_DIR}"
    echo "        Run 01_prepare_braindino_inputs.py first."; exit 1
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
    echo "  last_checkpoint.pt found -- will AUTO-RESUME."
fi

python "${SCRIPT_DIR}/02_supervised_finetuning_BrainDINO.py" \
    --task                  "${TASK}" \
    --seed                  "${SEED}" \
    --strategy              "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --augment               "${AUGMENT}" \
    --wandb \
    --wandb_project         "${WANDB_PROJECT}" \
    --pretrained_ckpt       "${PRETRAINED_CKPT}" \
    --matched_labels_csv    "${MATCHED_LABELS_CSV}" \
    --data_dir              "${DATA_DIR}" \
    --braindino_inputs_dir  "${BRAINDINO_INPUTS_DIR}" \
    --out_dir               "${OUT_DIR}" \
    --num_workers           "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_OUT_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
