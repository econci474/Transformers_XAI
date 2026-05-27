#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bdino_full_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/ft/slurm_logs/bdino_full_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/ft/slurm_logs/bdino_full_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-44%4
# =============================================================================
# 02_finetune_BrainDINO_full_ft_submit_csd3.sh  --  end-to-end full fine-tune
# =============================================================================
# Fine-tunes BrainDINO ViT-B/16 end-to-end. 86 M trainable params with no
# LoRA adapters; higher overfit risk than LoRA but the strongest comparison
# vs. the BrainMVP full_ft sweep.
#
# 5 tasks x 3 seeds x 3 augments = 45 combinations (array 0-44).
# Same array decoder as the LoRA sweep.
#
# Submit AFTER the LoRA sweep has validated the harness end-to-end -- this is
# the most GPU-expensive of the three BrainDINO sweeps (50 epochs vs. 8 for
# LoRA, 86 M vs. 0.5 M trainable params).
#
# Output layout: braindino_outputs/ft/aug_${AUGMENT}/BrainDINO_vitb16_full_ft/<task>/seed_<n>/
# (nested under ft/ so this sweep is kept separate from the existing
#  frozen + head sweep at braindino_outputs/aug_none/... and the
#  LoRA sweep at braindino_outputs/lora/...)
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/ft/slurm_logs
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
# Strategy-scoped subroot: keeps this full-finetune sweep's outputs + slurm
# logs separate from the existing frozen + head sweep (OUT_DIR/aug_*) and
# the LoRA sweep (OUT_DIR/lora/...).
STRATEGY_OUT_DIR="${OUT_DIR}/ft"

PY_SESSION_FLAGS="--long all"
STRATEGY="full_ft"
WANDB_PROJECT="${WANDB_PROJECT:-braindino_full_ft}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Combination lookup: 4 tasks x 3 seeds x 3 augments = 36 ------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
AUGMENTS=("none" "stochastic" "plus_original")

N_TASKS=${#TASKS[@]}
N_SEEDS=${#SEEDS[@]}
N_AUG=${#AUGMENTS[@]}

ID=${SLURM_ARRAY_TASK_ID}
AUG_IDX=$(( ID % N_AUG ));      ID=$(( ID / N_AUG ))
SEED_IDX=$(( ID % N_SEEDS ));   ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
AUGMENT="${AUGMENTS[$AUG_IDX]}"

RUN_OUT_DIR="${STRATEGY_OUT_DIR}/aug_${AUGMENT}/BrainDINO_vitb16_${STRATEGY}/${TASK}/seed_${SEED}"
METRICS="${RUN_OUT_DIR}/metrics.json"
CKPT="${RUN_OUT_DIR}/last_checkpoint.pt"

mkdir -p "${RUN_OUT_DIR}"
mkdir -p "${STRATEGY_OUT_DIR}/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  BrainDINO ViT-B/16 full_ft (end-to-end)"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_AUG-1)))"
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
    --out_dir               "${STRATEGY_OUT_DIR}" \
    --num_workers           "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_OUT_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
