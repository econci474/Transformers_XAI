#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs/vit_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-35%4
# =============================================================================
# 04_finetune_ViT_submit_csd3.sh   --   ViT debugging re-run
# =============================================================================
# SLURM array job for supervised fine-tuning of the MAE-pretrained ViT-B/3D on
# ADNI T1w MRIs, with the debugging fixes (balanced-acc model selection, LLRD,
# gradient accumulation, pre-flight input gate, per-subject eval).
#
# 1 model x 4 tasks x 3 seeds x 3 strategies = 36 combinations.
# Array index decoder: STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task     : T1_binary | T1b_binary | T1c_binary | T2_multiclass
#   Seed     : 0 | 1 | 2
#   Strategy : full_ft | frozen | scratch
# scratch = random init, no MAE checkpoint -> output slug ViT_B_scratch.
#
# Runs offline wandb -> sync to project 'vit_debugging' after the sweep:
#   wandb sync <OUT_DIR>/ViT_B_*/*/seed_*/*/wandb/offline-run-*
#
# SLURM .log/.err and the per-run offline-wandb data both land under the
# OUTPUT tree, never in the scripts/git directory.
#
# Resume / partial sweeps:
#   sbatch mri_pipeline/04_finetune_ViT_submit_csd3.sh                # full sweep
#   sbatch --array=0     mri_pipeline/04_finetune_ViT_submit_csd3.sh  # smoke (1 config)
#   sbatch --array=0,3,7 mri_pipeline/04_finetune_ViT_submit_csd3.sh  # specific configs
# Already-completed runs (metrics.json present) are auto-skipped.
#
# Pre-flight on a CSD3 login node (run once before submitting):
#   1. mri conda env (torch 2.4.1+cu121, monai, nibabel, pandas, scikit-learn,
#      einops, tqdm, pyyaml, matplotlib, wandb).
#   2. Pretrained ckpt at ${PRETRAINED_CKPT}.
#   3. POST-EXCLUSION matched labels CSV at ${MATCHED_LABELS_CSV} and
#      POST-EXCLUSION clinical splits at ${DATA_DIR}/seed_{0,1,2}/ -- upload the
#      no_cdr_stratified_post_exclusion artifacts to HPC and adjust the paths
#      below if they differ.
#   4. ViT inputs at ${VIT_INPUTS_DIR} populated by 03 (--long all coverage).
#   5. mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Offline wandb -- the compute nodes have no internet; sync after the sweep.
export WANDB_MODE=offline

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
# Post-exclusion clinical splits + matched labels (both uploaded to HPC).
# DATA_DIR must contain seed_{0,1,2}/{train,val,test}.csv.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_debug"

# -- Session selection (must match what 03 produced) --------------------------
PY_SESSION_FLAGS="--long all"      # every available session (bl..mAll)

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T2_multiclass")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen" "scratch")

N_TASKS=${#TASKS[@]}        # 4
N_SEEDS=${#SEEDS[@]}        # 3
N_STRAT=${#STRATEGIES[@]}   # 3

# Decode SLURM_ARRAY_TASK_ID -> (task, seed, strategy)
ID=${SLURM_ARRAY_TASK_ID}
STRAT_IDX=$(( ID % N_STRAT ));   ID=$(( ID / N_STRAT ))
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="${STRATEGIES[$STRAT_IDX]}"

# scratch runs use a distinct model slug (matches 04_supervised_finetuning_ViT.py)
if [ "${STRATEGY}" = "scratch" ]; then
    MODEL_SLUG="ViT_B_scratch"
else
    MODEL_SLUG="ViT_B_mae75"
fi
RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"

# Offline wandb data for this run lands beside its checkpoints (in the output
# tree), not in the scripts directory.
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT supervised fine-tuning -- debugging re-run"
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
echo "  wandb          : offline -> project vit_debugging (WANDB_DIR=${WANDB_DIR})"
echo "  Started        : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ "${STRATEGY}" != "scratch" ] && [ ! -f "${PRETRAINED_CKPT}" ]; then
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
    echo "        Upload the post-exclusion extended master to HPC."
    exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"
    echo "        Upload the no_cdr_stratified_post_exclusion splits to HPC."
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
    --grad_accum_steps    8 \
    --wandb \
    --wandb_project       vit_debugging \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "  Output         : ${RUN_DIR}/"
echo "============================================================"

exit ${EXIT_CODE}
