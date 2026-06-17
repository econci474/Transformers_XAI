#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft_v2s
#SBATCH --output=logs/vit_ft_v2s_%A_%a.log
#SBATCH --error=logs/vit_ft_v2s_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-11
# =============================================================================
# 04e_finetune_ViT_T1_T2_long_all_stratified_viscode2_submit_csd3.sh
# =============================================================================
# Fine-tunes ViT_B_mae75 on T1 (binary CN vs MCI+AD) and T2 (3-class CN/MCI/AD)
# using ALL longitudinal sessions (--long all) with:
#
#   1. VISCODE2-aligned MRI<->clinical matching  (03c master, viscode2_exact)
#   2. No-CDR stratified baseline splits (Label_visit_diag per visit)
#   3. Subject-level 80/10/10 train/val/test from clinical_pipeline stratified
#      splits at baseline (seeds 0, 1, 2)
#
# Label resolution (session_policy='current' in 04_supervised_finetuning.py):
#   - ses-bl   -> Label_bl_multi  (baseline diagnosis)
#   - ses-m12+ -> Label_visit_diag (per-visit diagnosis)
#
# 2 tasks x 3 seeds x 2 strategies = 12 combinations.
# Array index decoder:
#   STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride)
#
# Companion scripts:
#   04d  = T2-only, --long all, VISCODE_long alignment (old)
#   04e  = T1 + T2, --long all, VISCODE2 alignment, no-CDR stratified (this)
#
# Resume / partial sweeps:
#   sbatch mri_pipeline/04e_finetune_ViT_T1_T2_long_all_stratified_viscode2_submit_csd3.sh
#   sbatch --array=0 mri_pipeline/04e_finetune_ViT_T1_T2_long_all_stratified_viscode2_submit_csd3.sh
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"

# KEY CHANGE: VISCODE2-aligned master from 03c
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/viscode_2_aligned/master_mri_clinical_matched_viscode2.csv"

# KEY CHANGE: no-CDR stratified baseline splits
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline_stratified_no_cdr"

# Separate output dir for VISCODE2-aligned stratified runs
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_viscode2_stratified"

# -- Session selection: every available session per subject --------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup: T1 + T2, 3 seeds, 2 strategies -----------------------
TASKS=("T1_binary" "T2_multiclass")
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
echo "  ViT fine-tuning — T1+T2, VISCODE2 stratified, --long all"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Sessions       : ${PY_SESSION_FLAGS}"
echo "  Matched CSV    : ${MATCHED_LABELS_CSV}"
echo "  Split dir      : ${DATA_DIR}"
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

python "${SCRIPT_DIR}/04_supervised_finetuning.py" \
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
