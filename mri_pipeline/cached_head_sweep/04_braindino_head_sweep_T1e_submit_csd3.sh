#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=bdino_head_T1e
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_head_T1e_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_head_T1e_%A_%a.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-53%4
# =============================================================================
# 04_braindino_head_sweep_T1e_submit_csd3.sh
# =============================================================================
# Cached-embedding head HP sweep on the NEW T1e_pcn_vs_scn task. Identical
# scope as 04_brainmvp_head_sweep_T1e_submit_csd3.sh but with BrainDINO
# 768-d frozen-encoder embeddings.
#
# 1 task x 3 seeds x 3 LRs x 3 drops x 2 LSs = 54 cells.
# T1e is registered in mri_pipeline/04_supervised_finetuning.py
# TASK_CONFIG with pos_cols=["pCN_to_AD","pCN_to_MCI"] / neg_col="sCN".
#
# Submit:  sbatch mri_pipeline/cached_head_sweep/04_braindino_head_sweep_T1e_submit_csd3.sh
# =============================================================================

module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/cached_head_sweep"
EMBEDDINGS_PT="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/braindino_pooled.pt"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/aug_none_hp_tuned/BrainDINO_vitb16_frozen_cached"

MODEL_NAME="BrainDINO"
EMBED_DIM=768
WANDB_PROJECT="${WANDB_PROJECT:-braindino_frozen_cached_T1e}"
export WANDB_MODE="${WANDB_MODE:-offline}"

TASK="T1e_pcn_vs_scn"
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
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  ${MODEL_NAME} cached-head sweep -- T1e"
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
