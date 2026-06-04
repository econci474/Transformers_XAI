#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=vit_head_T4
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_head_T4_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_head_T4_%A_%a.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-53%4
# =============================================================================
# 04_vit_mae_head_sweep_T4_submit_csd3.sh
# =============================================================================
# Cached-head HP sweep on T4_conv_horizon (3-class ordinal: <3y / 3-7y / >=7y
# conversion bucket for pMCI + pCN_to_AD) using ViT-MAE 768-d frozen embeddings
# (from cached_embeddings/vit_mae/vit_mae_pooled.pt). Mirrors the BrainMVP /
# BrainDINO T4 cached sweeps -- only model_name + embed_dim + paths differ.
#
# 1 task x 3 seeds x 3 LRs x 3 drops x 2 LSs = 54 cells.
#
# T4-specific inputs:
#   --data_dir  = baseline_T4 (NOT the default baseline), built by
#                 clinical_pipeline/01e_build_T4_labels_and_splits.py.
#   --matched_labels_csv has the Label_T4 column (added by 01e). T4 is
#                 baseline_anchored with label_anchor_max_months=12, so the
#                 abundant m12 scans (not just sparse bl) are used.
#
# Submit:  sbatch mri_pipeline/cached_head_sweep/04_vit_mae_head_sweep_T4_submit_csd3.sh
# =============================================================================

module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/cached_head_sweep"
EMBEDDINGS_PT="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/vit_mae/vit_mae_pooled.pt"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
# T4-specific split tree (built by 01e_build_T4_labels_and_splits.py):
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline_T4"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/aug_none/ViT_B_mae75_frozen_cached"

MODEL_NAME="ViT-MAE75"
EMBED_DIM=768
WANDB_PROJECT="${WANDB_PROJECT:-vit_mae_frozen_cached_T4}"
export WANDB_MODE="${WANDB_MODE:-offline}"

TASK="T4_conv_horizon"
SEEDS=(0 1 2)
DROPOUTS=(0.1 0.2 0.3)
LABEL_SMOOTHINGS=(0.0 0.1)
LRS=(1e-3 1e-4 1e-5)

N_SEEDS=${#SEEDS[@]};   N_DROP=${#DROPOUTS[@]}
N_LS=${#LABEL_SMOOTHINGS[@]}; N_LR=${#LRS[@]}

ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));   ID=$(( ID / N_SEEDS ))
DROP_IDX=$(( ID % N_DROP ));    ID=$(( ID / N_DROP ))
LS_IDX=$(( ID % N_LS ));        ID=$(( ID / N_LS ))
LR_IDX=$(( ID % N_LR ))

SEED="${SEEDS[$SEED_IDX]}"
DROPOUT="${DROPOUTS[$DROP_IDX]}"
LS="${LABEL_SMOOTHINGS[$LS_IDX]}"
LR="${LRS[$LR_IDX]}"

RUN_OUT_DIR="${OUT_DIR}/${TASK}/seed_${SEED}/lr$(printf '%.0e' $LR)_d${DROPOUT}_ls${LS}"
METRICS="${RUN_OUT_DIR}/metrics.json"

mkdir -p "${RUN_OUT_DIR}"
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  ${MODEL_NAME} cached-head sweep -- T4 (3-class ordinal)"
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Array ID  : $SLURM_ARRAY_TASK_ID (of 0-$((N_SEEDS*N_DROP*N_LS*N_LR-1)))"
echo "  Task      : ${TASK}    Seed: ${SEED}"
echo "  HP        : lr=${LR}  dropout=${DROPOUT}  ls=${LS}  wd=1e-5"
echo "  Output    : ${RUN_OUT_DIR}"
echo "  wandb     : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started   : $(date)"
echo "============================================================"

if [ ! -f "${EMBEDDINGS_PT}" ]; then
    echo "[ERROR] Cached embeddings not found: ${EMBEDDINGS_PT}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] T4 splits dir not found: ${DATA_DIR}/seed_${SEED}"
    echo "        Run clinical_pipeline/01e_build_T4_labels_and_splits.py first."
    exit 1
fi
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi

python "${SCRIPT_DIR}/04_head_finetune_from_embeddings.py" \
    --model_name           "${MODEL_NAME}" \
    --embeddings_pt        "${EMBEDDINGS_PT}" \
    --embed_dim            "${EMBED_DIM}" \
    --task                 "${TASK}" \
    --seed                 "${SEED}" \
    --lr                   "${LR}" \
    --weight_decay         1e-5 \
    --drop_rate            "${DROPOUT}" \
    --label_smoothing      "${LS}" \
    --long                 all \
    --matched_labels_csv   "${MATCHED_LABELS_CSV}" \
    --data_dir             "${DATA_DIR}" \
    --out_dir              "${OUT_DIR}" \
    --wandb \
    --wandb_project        "${WANDB_PROJECT}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
