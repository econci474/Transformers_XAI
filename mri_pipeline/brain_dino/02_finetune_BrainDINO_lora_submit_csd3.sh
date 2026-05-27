#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bdino_lora
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/lora/slurm_logs/bdino_lora_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/lora/slurm_logs/bdino_lora_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-44%4
# =============================================================================
# 02_finetune_BrainDINO_lora_submit_csd3.sh  --  LoRA rank-8 finetune (paper)
# =============================================================================
# Fine-tunes BrainDINO ViT-B/16 with LoRA-rank-8 adapters injected into every
# attention layer's qkv + proj -- the paper's "runners/classification.py"
# canonical recipe (lora_rank=8, frozen backbone, ReduceLROnPlateau).
#
# 5 tasks x 3 seeds x 3 augments = 45 combinations (array 0-44).
# Array index decoder: AUG_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task    : T1_binary | T1b_binary | T1c_binary | T1d_binary | T2_multiclass
#   Seed    : 0 | 1 | 2
#   Augment : none | stochastic | plus_original
#
# This sweep is the CHEAPEST of the three strategies (LoRA = ~0.5 M trainable
# params, 8 epochs vs. 50 for full_ft). Submit this FIRST to validate the
# harness end-to-end before committing GPU-hours to the full_ft sweep.
#
# Output layout: braindino_outputs/lora/aug_${AUGMENT}/BrainDINO_vitb16_lora/<task>/seed_<n>/
# (nested under lora/ so this sweep is kept separate from the existing
#  frozen + head sweep at braindino_outputs/aug_none/...)
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/lora/slurm_logs
#   # Ensure peft is installed in the mri env:
#   conda activate mri && pip install peft>=0.7
#   # Ensure the checkpoint is at the canonical HPC path:
#   ls /home/ec474/rds/hpc-work/brain_dino_model.pth
#
# Usage:  sbatch mri_pipeline/brain_dino/02_finetune_BrainDINO_lora_submit_csd3.sh
# Smoke:  sbatch --array=0 mri_pipeline/brain_dino/02_finetune_BrainDINO_lora_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
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
# Strategy-scoped subroot: keeps this LoRA sweep's outputs + slurm logs
# separate from the existing frozen + head sweep that lives under OUT_DIR/aug_*.
STRATEGY_OUT_DIR="${OUT_DIR}/lora"

PY_SESSION_FLAGS="--long all"
STRATEGY="lora"
LORA_RANK="${LORA_RANK:-8}"
WANDB_PROJECT="${WANDB_PROJECT:-braindino_lora}"
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

# Per-aug output root so the trainer's path schema and the aggregator glob
# (braindino_outputs/aug_*/BrainDINO_vitb16_*/...) line up.
RUN_OUT_DIR="${STRATEGY_OUT_DIR}/aug_${AUGMENT}/BrainDINO_vitb16_${STRATEGY}/${TASK}/seed_${SEED}"
METRICS="${RUN_OUT_DIR}/metrics.json"
CKPT="${RUN_OUT_DIR}/last_checkpoint.pt"

mkdir -p "${RUN_OUT_DIR}"
mkdir -p "${STRATEGY_OUT_DIR}/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  BrainDINO ViT-B/16 LoRA-${LORA_RANK} finetune"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_AUG-1)))"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  Strategy    : ${STRATEGY}  (rank=${LORA_RANK})"
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
    --lora_rank             "${LORA_RANK}" \
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
