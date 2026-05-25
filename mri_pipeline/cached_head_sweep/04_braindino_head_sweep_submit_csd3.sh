#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=bdino_head
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_head_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_head_%A_%a.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-269%16
# =============================================================================
# 04_braindino_head_sweep_submit_csd3.sh  --  cached-embedding head HP sweep
# =============================================================================
# Trains a 2-layer MLP head on the BrainDINO cached pooled-CLS embeddings
# (extracted by 03_extract_braindino_embeddings.py). Sweeps lr x dropout
# x label_smoothing across all 5 tasks x 3 seeds. WD fixed at 1e-5
# (paper default; not swept).
#
# 5 tasks x 3 seeds x 3 LRs x 3 dropouts x 2 LSs = 270 cells (array 0-269).
# Decoder ordering (task varies fastest -- results trickle across tasks
# in the first 16 array elements):
#   TASK_IDX -> SEED_IDX -> DROPOUT_IDX -> LS_IDX -> LR_IDX
#
# Each cell trains in ~2-5 min CPU. Full sweep ~50 min at %16 throttle.
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs
#   ls /home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/braindino_pooled.pt
#
# Submit: sbatch mri_pipeline/cached_head_sweep/04_braindino_head_sweep_submit_csd3.sh
# Smoke : sbatch --array=0 ...                # one cell first
# Resume: re-running auto-skips cells with metrics.json present.
# =============================================================================

module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/cached_head_sweep"
EMBEDDINGS_PT="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/braindino_pooled.pt"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/aug_none/BrainDINO_vitb16_frozen_cached"

MODEL_NAME="BrainDINO"
EMBED_DIM=768
WANDB_PROJECT="${WANDB_PROJECT:-braindino_frozen_cached}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- HP grid (WD fixed at 1e-5, not swept) ------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
DROPOUTS=(0.1 0.2 0.3)
LABEL_SMOOTHINGS=(0.0 0.1)
LRS=(1e-3 1e-4 1e-5)

N_TASKS=${#TASKS[@]};       N_SEEDS=${#SEEDS[@]};       N_DROP=${#DROPOUTS[@]}
N_LS=${#LABEL_SMOOTHINGS[@]}; N_LR=${#LRS[@]}
# Total: N_TASKS*N_SEEDS*N_DROP*N_LS*N_LR = 5*3*3*2*3 = 270

# Decoder: task fastest -> seed -> dropout -> ls -> lr
ID=${SLURM_ARRAY_TASK_ID}
TASK_IDX=$(( ID % N_TASKS ));   ID=$(( ID / N_TASKS ))
SEED_IDX=$(( ID % N_SEEDS ));   ID=$(( ID / N_SEEDS ))
DROP_IDX=$(( ID % N_DROP ));    ID=$(( ID / N_DROP ))
LS_IDX=$(( ID % N_LS ));        ID=$(( ID / N_LS ))
LR_IDX=$(( ID % N_LR ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
DROPOUT="${DROPOUTS[$DROP_IDX]}"
LS="${LABEL_SMOOTHINGS[$LS_IDX]}"
LR="${LRS[$LR_IDX]}"

# Output slug matches the trainer's internal layout:
#   <out_dir>/<task>/seed_<n>/lr<LR>_d<D>_ls<LS>/metrics.json
RUN_OUT_DIR="${OUT_DIR}/${TASK}/seed_${SEED}/lr$(printf '%.0e' $LR)_d${DROPOUT}_ls${LS}"
METRICS="${RUN_OUT_DIR}/metrics.json"

mkdir -p "${RUN_OUT_DIR}"
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs"
export WANDB_DIR="${RUN_OUT_DIR}"

echo "============================================================"
echo "  ${MODEL_NAME} cached-head sweep"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_DROP*N_LS*N_LR-1)))"
echo "  Task        : ${TASK}    Seed: ${SEED}"
echo "  HP          : lr=${LR}  dropout=${DROPOUT}  ls=${LS}  wd=1e-5 (fixed)"
echo "  Output dir  : ${RUN_OUT_DIR}"
echo "  wandb       : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${EMBEDDINGS_PT}" ]; then
    echo "[ERROR] Cached embeddings not found: ${EMBEDDINGS_PT}"
    echo "        Run 03_extract_braindino_embeddings.py first."; exit 1
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
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
