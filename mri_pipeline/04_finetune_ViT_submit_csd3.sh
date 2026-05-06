#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft
#SBATCH --output=logs/vit_ft_%A_%a.log
#SBATCH --error=logs/vit_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-23
# =============================================================================
# 04_finetune_ViT_submit_csd3.sh
# =============================================================================
# SLURM array job for supervised fine-tuning of the MAE-pretrained ViT-B/3D
# on ADNI T1w MRIs. Mirrors clinical_pipeline/03_encoder_submit_csd3.sh.
#
# 1 model x 4 tasks x 3 seeds x 2 strategies = 24 combinations.
# Array index decoder: STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task     : T1_binary | T2_multiclass | T3a_conv3y | T3b_conv5y
#   Seed     : 0 | 1 | 2
#   Strategy : full_ft | frozen
#
# Resume / partial sweeps:
#   sbatch mri_pipeline/04_finetune_ViT_submit_csd3.sh                       # full sweep
#   sbatch --array=0     mri_pipeline/04_finetune_ViT_submit_csd3.sh         # smoke (1 config)
#   sbatch --array=0,3,7 mri_pipeline/04_finetune_ViT_submit_csd3.sh         # specific configs
# Already-completed runs (metrics.json present) are auto-skipped.
#
# Pre-flight on a CSD3 login node (run once before submitting):
#   1. mri conda env (torch 2.4.1+cu121, monai, nibabel, pandas, scikit-learn,
#      einops, tqdm, pyyaml, matplotlib).
#   2. Pretrained ckpt at ${PRETRAINED_CKPT}.
#   3. Matched labels CSV at ${MATCHED_LABELS_CSV} (built by 03b on local or HPC).
#   4. ViT inputs at ${VIT_INPUTS_DIR} populated by 03 with the appropriate
#      session coverage (--long N) for the tasks being run.
#   5. Clinical splits at ${DATA_DIR}/seed_{0,1,2}/{train,val,test}.csv.
#   6. mkdir -p logs
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/master_mri_clinical_matched_labels.csv"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs"

# -- Session selection (must match what 03 produced) --------------------------
# Uncomment ONE. Tasks with session_policy='baseline_only' (T3a, T3b)
# internally restrict to ses-bl rows regardless of this flag — safe to leave
# longitudinal here.
#PY_SESSION_FLAGS="--session bl"   # ses-bl only (smoke / baseline-only sweep)
PY_SESSION_FLAGS="--long 1"        # bl + m12  (matches 03 --long 1 output)
#PY_SESSION_FLAGS="--long 3"       # bl through m36
#PY_SESSION_FLAGS="--long 5"       # bl through m60
#PY_SESSION_FLAGS="--long all"     # every available session

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T2_multiclass" "T3a_conv3y" "T3b_conv5y")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_TASKS=${#TASKS[@]}        # 4
N_SEEDS=${#SEEDS[@]}        # 3
N_STRAT=${#STRATEGIES[@]}   # 2

# Decode SLURM_ARRAY_TASK_ID -> (task, seed, strategy)
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
echo "  ViT supervised fine-tuning -- array sweep"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Sessions       : ${PY_SESSION_FLAGS}"
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
    echo "        Run 03_prepare_ViT_submit_csd3.sh first."
    exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"
    echo "        Run 03b_match_mri_to_clinical.py first."
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
