#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=bdino_T2wide
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_T2wide_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs/bdino_T2wide_%A_%a.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-575%16
# =============================================================================
# 04_braindino_head_sweep_T2_wide_submit_csd3.sh
#   WIDER frozen linear-probe HP sweep for BrainDINO on T2_multiclass.
# =============================================================================
# BrainDINO is used FROZEN (the published ~90% used a frozen encoder + a
# lightweight head). Our 3-class T2 cached-head plateaued at val-bACC ~0.60,
# so this widens the head/optimisation search on the SAME cached pooled-CLS
# embeddings (cheap, CPU): head {mlp, linear} x lr {1e-2,1e-3,1e-4} x
# weight_decay {1e-5,1e-2} x drop {0.1,0.3} x label_smoothing {0.0,0.1} x
# loss {ce_weighted, focal} x standardize {off,on}  =  192 cells x 3 seeds.
# Run ALL 3 seeds directly; the winner is selected by MEAN val-bACC across the
# 3 seeds (05b_select_best_hp_per_task.py) -> avoids HP overfitting to one seed.
#
# Cells write metrics.json ONLY (--metrics_only): no checkpoints/manifests for
# losing cells. The winning HP's checkpoint + embeddings are produced later by
# the finalize step.
#
# Decoder ordering (seed fastest):
#   SEED -> HEAD -> LR -> WD -> DROP -> LS -> LOSS -> STD
# Total = 3*2*3*2*2*2*2*2 = 576 (array 0-575).
#
# Pre-flight (once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs
#   ls /home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/braindino_pooled.pt
# Submit : sbatch mri_pipeline/cached_head_sweep/04_braindino_head_sweep_T2_wide_submit_csd3.sh
# Smoke  : sbatch --array=0 ...           # one cell first
# Resume : re-running auto-skips cells whose metrics.json already exists.
# =============================================================================

module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/cached_head_sweep"
EMBEDDINGS_PT="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/braindino_pooled.pt"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/aug_none_T2_wide/BrainDINO_vitb16_frozen_cached"

MODEL_NAME="BrainDINO"
EMBED_DIM=768
TASK="T2_multiclass"
WANDB_PROJECT="${WANDB_PROJECT:-braindino_T2_wide}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- WIDER HP grid ------------------------------------------------------------
SEEDS=(0 1 2)
HEADS=(mlp linear)
LRS=(1e-2 1e-3 1e-4)
WDS=(1e-5 1e-2)
DROPS=(0.1 0.3)
LSS=(0.0 0.1)
LOSSES=(ce_weighted focal)
STDS=(0 1)

N_SEEDS=${#SEEDS[@]}; N_HEAD=${#HEADS[@]}; N_LR=${#LRS[@]}; N_WD=${#WDS[@]}
N_DROP=${#DROPS[@]}; N_LS=${#LSS[@]}; N_LOSS=${#LOSSES[@]}; N_STD=${#STDS[@]}

# Decoder: seed fastest -> head -> lr -> wd -> drop -> ls -> loss -> std
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS )); ID=$(( ID / N_SEEDS ))
HEAD_IDX=$(( ID % N_HEAD ));  ID=$(( ID / N_HEAD ))
LR_IDX=$((   ID % N_LR ));    ID=$(( ID / N_LR ))
WD_IDX=$((   ID % N_WD ));    ID=$(( ID / N_WD ))
DROP_IDX=$(( ID % N_DROP ));  ID=$(( ID / N_DROP ))
LS_IDX=$((   ID % N_LS ));    ID=$(( ID / N_LS ))
LOSS_IDX=$(( ID % N_LOSS ));  ID=$(( ID / N_LOSS ))
STD_IDX=$((  ID % N_STD ))

SEED="${SEEDS[$SEED_IDX]}";  HEAD="${HEADS[$HEAD_IDX]}";  LR="${LRS[$LR_IDX]}"
WD="${WDS[$WD_IDX]}";        DROP="${DROPS[$DROP_IDX]}";  LS="${LSS[$LS_IDX]}"
LOSS="${LOSSES[$LOSS_IDX]}"; STD="${STDS[$STD_IDX]}"
STD_FLAG=""; [ "$STD" = "1" ] && STD_FLAG="--standardize"

mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/braindino_outputs/slurm_logs"
export WANDB_DIR="${OUT_DIR}/${TASK}/seed_${SEED}"
mkdir -p "${WANDB_DIR}"

echo "============================================================"
echo "  BrainDINO T2 WIDE frozen-probe sweep"
echo "  Array ID : ${SLURM_ARRAY_TASK_ID} (of 0-$((N_SEEDS*N_HEAD*N_LR*N_WD*N_DROP*N_LS*N_LOSS*N_STD-1)))"
echo "  seed=${SEED} head=${HEAD} lr=${LR} wd=${WD} drop=${DROP} ls=${LS} loss=${LOSS} std=${STD}"
echo "  Started  : $(date)"
echo "============================================================"

if [ ! -f "${EMBEDDINGS_PT}" ]; then
    echo "[ERROR] Cached embeddings not found: ${EMBEDDINGS_PT}"; exit 1
fi

# The trainer self-skips when its computed metrics.json already exists, so no
# bash-side slug duplication is needed here.
python "${SCRIPT_DIR}/04_head_finetune_from_embeddings.py" \
    --model_name           "${MODEL_NAME}" \
    --embeddings_pt        "${EMBEDDINGS_PT}" \
    --embed_dim            "${EMBED_DIM}" \
    --task                 "${TASK}" \
    --seed                 "${SEED}" \
    --head                 "${HEAD}" \
    --lr                   "${LR}" \
    --weight_decay         "${WD}" \
    --drop_rate            "${DROP}" \
    --label_smoothing      "${LS}" \
    --loss                 "${LOSS}" \
    ${STD_FLAG} \
    --metrics_only \
    --long                 all \
    --matched_labels_csv   "${MATCHED_LABELS_CSV}" \
    --data_dir             "${DATA_DIR}" \
    --out_dir              "${OUT_DIR}" \
    --wandb \
    --wandb_project        "${WANDB_PROJECT}"

EXIT_CODE=$?
echo "  Finished : $(date)   exit=${EXIT_CODE}"
exit ${EXIT_CODE}
