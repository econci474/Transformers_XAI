#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_T4_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_T4_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_T4_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%4
# =============================================================================
# 04i_finetune_BrainMVP_T4_full_ft_submit_csd3.sh
# =============================================================================
# BrainMVP UniFormer full fine-tuning on T4_conv_horizon (3-class ordinal:
# AD-conversion-horizon bucket for pMCI + pCN_to_AD, 146 subjects).
# Submitted ONCE per augment mode; --array=0-2 covers the 3 seeds.
#
# Total per submission: 1 task x 3 seeds x 1 strategy = 3 array cells.
# Wall-clock per cell: ~30 min (small cohort) on A100.
#
# Submit per augment:
#   sbatch --export=ALL,AUGMENT=none          mri_pipeline/brain_mvp/04i_finetune_BrainMVP_T4_full_ft_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=stochastic    mri_pipeline/brain_mvp/04i_finetune_BrainMVP_T4_full_ft_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=plus_original mri_pipeline/brain_mvp/04i_finetune_BrainMVP_T4_full_ft_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

AUGMENT="${AUGMENT:-stochastic}"
AUG_COPIES="${AUG_COPIES:-1}"
USE_WANDB="${USE_WANDB:-1}"
WANDB_PROJECT="${WANDB_PROJECT:-brainmvp_T4}"
export WANDB_MODE="${WANDB_MODE:-offline}"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"

# T4-specific master + splits (Label_T4 added in-place by 01e_build_T4_labels_and_splits.py)
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline_T4"

OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/aug_${AUGMENT}"

# T4 uses session_policy=baseline_anchored with label_anchor_max_months=0
# (bl scans only), but --long all is harmless because the trainer's
# session-policy logic filters to <=0 months internally.
PY_SESSION_FLAGS="--long all"

TASK="T4_conv_horizon"
SEEDS=(0 1 2)
STRATEGY="full_ft"

SEED="${SEEDS[$SLURM_ARRAY_TASK_ID]}"

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
echo "  BrainMVP full_ft -- T4 (3-class ordinal)"
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Array ID  : $SLURM_ARRAY_TASK_ID  (seed ${SEED} of ${#SEEDS[@]})"
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
    echo "[ERROR] T4 splits dir not found: ${DATA_DIR}/seed_${SEED}"
    echo "        Run clinical_pipeline/01e_build_T4_labels_and_splits.py first."
    exit 1
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
