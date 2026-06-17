#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_hi_lr
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs/vit_hi_lr_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs/vit_hi_lr_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-29%4
# =============================================================================
# 04h_finetune_ViT_hi_lr_submit_csd3.sh
# =============================================================================
# COMPARISON sweep -- isolated in its own W&B project `vit_mae_hi_lr` so the
# new run set doesn't visually mix with the prior 24 runs in `vit_debugging`.
# Same trainer
# (04_supervised_finetuning.py), same architecture (ViT-B/16 mae75), with
# a single-variable change:
#
#   --lr 1e-3            (vs the trainer's full_ft default of 1e-4)
#   --llrd_gamma 0.7     (UNCHANGED -- keeps the layer-wise LR shield
#                         that was added in commit 2d74f99 to prevent
#                         full_ft collapse; do NOT drop this without a
#                         separate ablation)
#
# Motivation
# ----------
# The 24 prior runs in vit_debugging used lr=1e-4 + llrd_gamma=0.7. Across
# all tasks those runs underperformed the cached-head sweep at lr=1e-3
# (T1c: 0.844 head-sweep vs 0.782 vit_debugging). The most parsimonious
# hypothesis is that the HEAD was training too slowly -- not that LLRD was
# wrong. With lr=1e-3 + llrd=0.7:
#
#   head            -> lr=1e-3                (matches winning head-sweep)
#   last block (12) -> lr=1e-3
#   block 6         -> lr=1e-3 * 0.7^6  = 1.2e-4
#   patch_embed     -> lr=1e-3 * 0.7^12 = 1.4e-5  (still protected)
#
# So the head trains as fast as the winning head-only sweep, but the deep
# encoder layers still get the protective LR decay that prevented collapse
# in commit 2d74f99. If this works we get a single-variable result;
# if it collapses again the issue isn't LR, it's something else.
#
# 5 tasks x 3 seeds x 2 strategies = 30 combinations (array 0-29).
# Array decoder: STRAT_IDX -> SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task    : T1_binary | T1b_binary | T1c_binary | T1d_binary | T2_multiclass
#   Seed    : 0 | 1 | 2
#   Strat   : full_ft | frozen
#
# Strategy notes:
#   - full_ft: lr=1e-3 + llrd_gamma=0.7. This is the actual variable test.
#   - frozen: lr=1e-3 is already the trainer's frozen default; LLRD doesn't
#     apply when the encoder is frozen. These rows will look very similar
#     to the prior frozen runs in vit_debugging (same recipe). Re-running
#     them gets us T1d coverage which the prior 24 runs skipped, and lands
#     them in the new output tree alongside the full_ft sister cells.
#
# Output layout: vit_outputs_hi_lr/ViT_B_mae75/<task>/seed_<n>/<strategy>/
# (new tree to avoid colliding with the prior runs' output dirs)
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs
#
# Submit:  sbatch mri_pipeline/04h_finetune_ViT_hi_lr_submit_csd3.sh
# Smoke :  sbatch --array=0 mri_pipeline/04h_finetune_ViT_hi_lr_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

export WANDB_MODE="${WANDB_MODE:-offline}"
WANDB_PROJECT="${WANDB_PROJECT:-vit_mae_hi_lr}"

# -- Hardcoded paths (post-exclusion, current canonical) ----------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"

PY_SESSION_FLAGS="--long all"
LLRD_GAMMA="${LLRD_GAMMA:-1.0}"   # the variable under test; override allowed

# -- Run-time aug knob (baked into OUT_DIR so multiple augs coexist) ----------
# Default 'random' keeps existing aug=random rsync paths intact:
#   vit_outputs_hi_lr/ViT_B_mae75/...   <- aug=random (legacy default, kept)
# For any other aug, the path gains an aug_<augment>/ layer:
#   vit_outputs_hi_lr/aug_none/ViT_B_mae75/...
#   vit_outputs_hi_lr/aug_plus_original/ViT_B_mae75/...
# Override at submit time:
#   sbatch --export=ALL,AUGMENT=none          ...
#   sbatch --export=ALL,AUGMENT=plus_original ...
AUGMENT="${AUGMENT:-random}"
if [ "${AUGMENT}" = "random" ]; then
    OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr"
else
    OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/aug_${AUGMENT}"
fi
AUG_COPIES="${AUG_COPIES:-1}"     # plus_original only -- ignored for random/none

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
STRATEGIES=("full_ft" "frozen")

N_TASKS=${#TASKS[@]}        # 5
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

# Centralised slurm-log dir under the un-suffixed base so logs from all aug
# variants land together (easier to grep cross-aug failures from one place).
mkdir -p "/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs"
mkdir -p "${RUN_DIR}"
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-MAE75 finetune -- no-LLRD comparison sweep"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS*N_STRAT-1)))"
echo "  Node        : $SLURMD_NODENAME"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  Strategy    : ${STRATEGY}"
echo "  LLRD gamma  : ${LLRD_GAMMA}  (1.0 = no decay; prior runs used 0.7)"
echo "  Augment     : ${AUGMENT}     (aug_copies=${AUG_COPIES} -- only used for plus_original)"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : ${WANDB_PROJECT}  (WANDB_MODE=${WANDB_MODE})"
echo "  Started     : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs dir not found: ${VIT_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"; exit 1
fi

if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- python will AUTO-RESUME this combo."
fi

python "${SCRIPT_DIR}/04_supervised_finetuning.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --llrd_gamma          "${LLRD_GAMMA}" \
    --augment             "${AUGMENT}" \
    --aug_copies          "${AUG_COPIES}" \
    ${PY_SESSION_FLAGS} \
    --wandb \
    --wandb_project       "${WANDB_PROJECT}" \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
if [ ! -f "${METRICS}" ] && [ -f "${CKPT}" ]; then
    echo "  NOTE: no metrics.json yet (likely hit 12h cap). Resubmit to AUTO-RESUME."
fi
echo "============================================================"
exit ${EXIT_CODE}
