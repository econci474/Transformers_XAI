#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_v2g
#SBATCH --output=logs/vit_v2g_%A_%a.log
#SBATCH --error=logs/vit_v2g_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-17
# =============================================================================
# 04g_finetune_ViT_combined_submit_csd3.sh
# =============================================================================
# Combined replacement for 04e (T1_binary+T2_multiclass) + 04f (T1b_binary).
# Supersedes both; 04e/04f are kept only for provenance.
#
# Fine-tunes ViT_B_mae75 on T1_binary (CN vs MCI+AD), T1b_binary (CN+MCI vs AD)
# and T2_multiclass (CN/MCI/AD) using ALL longitudinal sessions (--long all):
#
#   1. VISCODE2-aligned MRI<->clinical matching  (03c master, viscode2_exact)
#   2. No-CDR stratified baseline splits (Label_visit_diag per visit)
#   3. Subject-level 80/10/10 train/val/test from clinical_pipeline stratified
#      splits at baseline (seeds 0, 1, 2)
#
# 3 tasks x 3 seeds x 2 strategies = 18 combinations (array 0-17).
# Array index decoder:  STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride)
#
# Publication fine-tuning recipe is kept UNCHANGED (epochs=50 full_ft / 70
# frozen, batch_size=4 — the python script's auto defaults). No --epochs /
# --batch_size override is passed for any augment mode.
#
# Augment modes (script-level env var AUGMENT):
#   none           no augmentation (eval transform on train)
#   random         on-the-fly stochastic aug (default; original behaviour)
#   plus_original  originals (no aug) + AUG_COPIES always-on augmented copies
#
# 12 h wall-time cap (hard CSD3 limit): plus_original full_ft (~2x cost at
# AUG_COPIES=1) may exceed one 12 h window. The python script writes an atomic
# per-epoch last_checkpoint.pt and AUTO-RESUMES on rerun, so simply resubmit
# the same array until every metrics.json exists (skip-if-exists makes
# completed combos no-ops). frozen runs fit one window.
#
# Usage:
#   # log to wandb project vit_orig_vs_aug (offline; wandb sync afterwards):
#   sbatch --export=ALL,AUGMENT=none,USE_WANDB=1            mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=plus_original,USE_WANDB=1   mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh
#   sbatch --export=ALL,AUGMENT=random,USE_WANDB=1          mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh   # optional baseline
#   # resume/finish a partial sweep (auto-resumes from last_checkpoint.pt):
#   sbatch --export=ALL,AUGMENT=plus_original,USE_WANDB=1   mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh
#   # chain windows automatically:
#   sbatch --dependency=afterany:<jobid> --export=ALL,AUGMENT=plus_original,USE_WANDB=1 mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh
#   # single-combo smoke:
#   sbatch --array=0 --export=ALL,AUGMENT=none,USE_WANDB=1  mri_pipeline/04g_finetune_ViT_combined_submit_csd3.sh
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Run-mode knobs (overridable via sbatch --export=ALL,VAR=...) -------------
AUGMENT="${AUGMENT:-random}"              # none | random | plus_original
AUG_COPIES="${AUG_COPIES:-1}"             # plus_original: K aug copies per orig
USE_WANDB="${USE_WANDB:-1}"               # 1 -> enable wandb logging
WANDB_PROJECT="${WANDB_PROJECT:-vit_orig_vs_aug}"
# CSD3 compute nodes have no outbound internet: log offline, `wandb sync` later.
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"

# VISCODE2-aligned master from 03c
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/viscode_2_aligned/master_mri_clinical_matched_viscode2.csv"

# No-CDR stratified baseline splits
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline_stratified_no_cdr"

# Separate output root PER AUGMENT MODE so none/random/plus_original do not
# collide and skip-if-exists stays per-mode. 05/06 aggregators run per root.
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_aug_${AUGMENT}"

# -- Session selection: every available session per subject --------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup: 3 tasks, 3 seeds, 2 strategies -----------------------
TASKS=("T1_binary" "T1b_binary" "T2_multiclass")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_TASKS=${#TASKS[@]}        # 3
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
CKPT="${RUN_DIR}/last_checkpoint.pt"

mkdir -p logs
mkdir -p "${OUT_DIR}"

# -- wandb flags --------------------------------------------------------------
WANDB_FLAGS=""
if [ "${USE_WANDB}" = "1" ]; then
    WANDB_FLAGS="--wandb --wandb_project ${WANDB_PROJECT}"
fi

echo "============================================================"
echo "  ViT fine-tuning — T1/T1b/T2, VISCODE2 stratified, --long all"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Node           : $SLURMD_NODENAME"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV       : $CUDA_VISIBLE_DEVICES"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Augment        : ${AUGMENT}  (aug_copies=${AUG_COPIES})"
echo "  wandb          : USE_WANDB=${USE_WANDB}  project=${WANDB_PROJECT}  WANDB_MODE=${WANDB_MODE}"
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
    echo "  metrics.json already exists -- skipping (run already complete)."
    exit 0
fi
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- python will AUTO-RESUME this combo."
fi

python "${SCRIPT_DIR}/04_supervised_finetuning.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --augment             "${AUGMENT}" \
    --aug_copies          "${AUG_COPIES}" \
    ${WANDB_FLAGS} \
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
if [ ! -f "${METRICS}" ] && [ -f "${CKPT}" ]; then
    echo "  NOTE: no metrics.json yet (likely hit the 12 h cap). Resubmit the"
    echo "        same sbatch line to AUTO-RESUME from last_checkpoint.pt."
fi
echo "============================================================"

exit ${EXIT_CODE}
