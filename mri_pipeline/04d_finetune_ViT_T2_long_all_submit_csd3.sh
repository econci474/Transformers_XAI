#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft_t2all
#SBATCH --output=logs/vit_ft_t2all_%A_%a.log
#SBATCH --error=logs/vit_ft_t2all_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-5
# =============================================================================
# 04d_finetune_ViT_T2_long_all_submit_csd3.sh
# =============================================================================
# Trains ViT_B_mae75 on the FULL longitudinal cohort (--long all, every
# available session per subject) for the 3-class CN/MCI/AD task only.
# Companion to:
#   04_finetune_ViT_submit_csd3.sh    (master sweep, --long 1, 24 jobs)
#   04b_finetune_ViT_anchored_submit_csd3.sh (T3a/T3b bl+m12, 12 jobs)
#   04c_finetune_ViT_bl_submit_csd3.sh       (T1/T2 bl-only, 12 jobs)
#
# Why this exists
#   T2 bl-only is impossible (N=28 across 3 classes -> classes missing in
#   train/val/test). T2 bl+m12 has ~526 scans which is workable but still
#   underpowered for a 3-class task at this batch size. Using --long all
#   scales the (image, label) pair count by ~3x without leaking subject
#   identity across the 80/10/10 splits, since the splits are subject-level
#   (Patient_ID-keyed at seeds 0/1/2 in clinical_pipeline/03b_make_splits).
#   T2 uses session_policy='current': each scan gets its own visit's
#   diagnosis label, so non-bl scans are valid training signal (no
#   baseline-anchor cap is needed for diagnosis classification).
#
# 1 task x 3 seeds x 2 strategies = 6 combinations.
# Array index decoder: STRAT_IDX -> SEED_IDX (low-to-high stride).
#
# Resume / partial sweeps:
#   sbatch mri_pipeline/04d_finetune_ViT_T2_long_all_submit_csd3.sh             # full 6-job sweep
#   sbatch --array=0  mri_pipeline/04d_finetune_ViT_T2_long_all_submit_csd3.sh  # smoke (seed0/full_ft)
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths ---------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/master_mri_clinical_matched_labels.csv"
# Sister output dir for the all-sessions sweep (kept separate from bl/bl+m12)
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_long_all"

# -- Session selection: every available session per subject ------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup: T2_multiclass only ----------------------------------
TASK="T2_multiclass"
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_SEEDS=${#SEEDS[@]}        # 3
N_STRAT=${#STRATEGIES[@]}   # 2

ID=${SLURM_ARRAY_TASK_ID}
STRAT_IDX=$(( ID % N_STRAT ));   ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ))

SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

MODEL_SLUG="ViT_B_mae75"
RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p logs
mkdir -p "${OUT_DIR}"

echo "============================================================"
echo "  ViT supervised fine-tuning -- T2 all-sessions sweep (--long all)"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_SEEDS*N_STRAT-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Sessions       : ${PY_SESSION_FLAGS}  (every available session, ~1500+ scans)"
echo "  Pretrained     : ${PRETRAINED_CKPT}"
echo "  Output dir     : ${RUN_DIR}"
echo "  Started        : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"
    exit 1
fi
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs directory not found: ${VIT_INPUTS_DIR}"
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

python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "  Output         : ${RUN_DIR}/"
echo "============================================================"

exit ${EXIT_CODE}
