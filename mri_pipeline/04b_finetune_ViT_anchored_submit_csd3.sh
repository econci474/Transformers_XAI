#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft_anch
#SBATCH --output=logs/vit_ft_anch_%A_%a.log
#SBATCH --error=logs/vit_ft_anch_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-11
# =============================================================================
# 04b_finetune_ViT_anchored_submit_csd3.sh
# =============================================================================
# Re-runs T3a/T3b only with --long 1 (bl + m12).
#
# Why a sister script instead of the master 04 sweep:
#   The master sweep (04_finetune_ViT_submit_csd3.sh) used --long 1, but T3a/T3b
#   were silently restricted to ses-bl by the legacy session_policy='baseline_only'.
#   Now that 04_supervised_finetuning_ViT.py uses session_policy='baseline_anchored'
#   with a m12 cap, T3a/T3b respect --long 1 (bl + m12 paired with the same
#   subject-level Label_3y / Label_5y). To preserve the old T3a/T3b bl-only
#   results for comparison, this script writes to a sister OUT_DIR.
#
# Why bl + m12 is safe:
#   Label_3y / Label_5y are baseline-anchored ("did the subject convert within
#   3y / 5y of bl"). At m12 there are still 24 / 48 months of headroom, so the
#   m12 scan does not reveal the answer. --long 3 (bl..m36) and --long 5
#   (bl..m60) WOULD cause leakage and are explicitly NOT supported here.
#
# 2 tasks x 3 seeds x 2 strategies = 12 combinations.
# Array index decoder: STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
# Resume / partial sweeps:
#   sbatch mri_pipeline/04b_finetune_ViT_anchored_submit_csd3.sh                 # full 12-job sweep
#   sbatch --array=0  mri_pipeline/04b_finetune_ViT_anchored_submit_csd3.sh      # smoke (T3a/seed0/full_ft)
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
# Sister output dir so old bl-only T3a/T3b runs are preserved
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_bl_m12"

# -- Session selection: bl + m12 only (do NOT change without recomputing leakage) --
PY_SESSION_FLAGS="--long 1"

# -- Combination lookup: T3a / T3b only ---------------------------------------
TASKS=("T3a_conv3y" "T3b_conv5y")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_TASKS=${#TASKS[@]}        # 2
N_SEEDS=${#SEEDS[@]}        # 3
N_STRAT=${#STRATEGIES[@]}   # 2

ID=${SLURM_ARRAY_TASK_ID}
STRAT_IDX=$(( ID % N_STRAT ));   ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

MODEL_SLUG="ViT_B_mae75"
RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p logs
mkdir -p "${OUT_DIR}"

echo "============================================================"
echo "  ViT supervised fine-tuning -- baseline-anchored sweep (bl+m12)"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Sessions       : ${PY_SESSION_FLAGS}  (bl + m12, capped at m12 for safety)"
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
