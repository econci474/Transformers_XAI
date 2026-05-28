#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_tiny
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_tiny_baseline/slurm_logs/vit_tiny_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_tiny_baseline/slurm_logs/vit_tiny_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%4
# =============================================================================
# 04i_finetune_ViT_tiny_scratch_submit_csd3.sh
# =============================================================================
# ViT-Tiny from-scratch ablation. ViT-Base from-scratch collapsed on multiple
# tasks (see `ViT-scratch / scratch` row in mri_pipeline/outputs/cross_model/
# cross_model_table.png -- 3 degenerate cells on T1, T1b, T1c). The hypothesis
# is that 86M params is far too much capacity for the ~580 train subjects.
# ViT-Tiny is ~5.5M params (~15x smaller), which sits closer to the regime
# where transformers can train from random init on small medical datasets
# (DeiT-Tiny on tiny medical cohorts is the canonical reference).
#
# Architecture (via --vit_size tiny in 04_supervised_finetuning_ViT.py):
#   embed_dim = 192,  depth = 12,  n_heads = 3,  mlp_ratio = 4.0
#   patch_size = 16,  img_size = 128^3,  global_avg_pool = False (CLS token)
#
# This sweep is STRICTLY from-scratch. The MAE-pretrained .pth is ViT-B only,
# so the trainer hard-errors if you combine --vit_size tiny with --strategy
# full_ft / frozen. The Tiny output tree is intentionally separate from the
# ViT-B-scratch tree (vit_baseline/) and from the MAE-pretrained tree
# (vit_outputs_*/) -- nothing collides.
#
# HP picks (informed by what we've seen so far in this thesis)
# -------------------------------------------------------------
# Pattern from prior collapses (AG-MS3D-vanilla at lr=3e-3, ViT-B-scratch):
# from-scratch 3D nets on 580 subjects need LR conservative, augmentation
# strong, and regularisation tuned to make label-smoothing + dropout do
# work the data alone can't. ViT-Tiny is small enough that we can afford
# aggressive augmentation + LS without overfitting; the bigger risk here
# is underfitting (head locks onto majority class in early epochs).
#
#   --lr 5e-4            scratch default is 1e-3; AG-MS3D rescue1 at 3e-3
#                        collapsed; pick a safer middle. DeiT-Tiny on
#                        ImageNet uses 1e-3 but ImageNet has 1000x more
#                        data than ADNI.
#   --weight_decay 5e-2  ~500x higher than the trainer default 1e-4. DeiT-
#                        Tiny canonical wd; without this, a from-scratch
#                        transformer on tiny data has nothing reining in
#                        the attention weight norms.
#   --label_smoothing 0.1   DeiT default; softens the majority-class minimum.
#   --augment plus_original  doubles the effective train set with deterministic
#                        K=1 augmented copies. The trainer's plus_original
#                        outperformed `random` on BrainMVP T1 / T1d / T2,
#                        and ViT-Tiny needs every effective-sample lift.
#   --aug_copies 1       K=1 -> 2x train set (1160 train items at T1).
#   --warmup_epochs 10   scratch default; with a 192-dim ViT a longer warmup
#                        helps the attention temperature settle before the LR
#                        bites.
#   --epochs 100         scratch default; Tiny converges faster than Base in
#                        wall-clock but the epoch count is fine.
#   --patience 15        slight bump over default 10; from-scratch transformers
#                        sometimes plateau then break through.
#   --llrd_gamma 1.0     uniform LR (scratch default). LLRD is a fine-tuning
#                        technique; doesn't apply here.
#   --drop_path_rate 0.1 / --attn_dropout 0.1   trainer defaults; preserve.
#   --grad_accum_steps 8  effective batch = 4 * 8 = 32 (trainer default).
#   --batch_size 4        Tiny fits easily on A100-80GB at bs=4; bumping
#                        wouldn't change effective batch (grad_accum
#                        compensates) so left at default.
#
# What to watch in wandb
# ----------------------
# - First 10 epochs: train_loss must drop below ~0.55 on binary, ~0.95 on
#   T2_multiclass. Anything still at ln(K) (~0.69 binary / ~1.10 multi) at
#   epoch 10 means the head locked onto majority -- same collapse mode as
#   AG-MS3D rescue1 / ViT-B-scratch.
# - val_bacc trajectory: bouncing between 0.50 and ~0.55 is fine in epochs
#   1-20 (transformer warmup); a real escape pushes val_bacc > 0.55 by
#   epoch ~25 and continues climbing.
# - Compare to ViT-B-scratch in the existing `vit_baseline` W&B project:
#   ViT-Tiny *should* outperform ViT-Base from-scratch on every task. If
#   it doesn't, the data really is the bottleneck, not the architecture.
#
# 5 tasks x 3 seeds = 15 combinations (array 0-14%4).
# Array decoder: SEED_IDX -> TASK_IDX (low-to-high stride).
#
#   Task : T1_binary | T1b_binary | T1c_binary | T1d_binary | T2_multiclass
#   Seed : 0 | 1 | 2
#
# W&B project: vit_tiny_scratch
# Output root: vit_tiny_baseline/ViT_T_scratch/<task>/seed_<n>/scratch/
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/vit_tiny_baseline/slurm_logs
#
# Submit:  sbatch mri_pipeline/baseline_vit/04i_finetune_ViT_tiny_scratch_submit_csd3.sh
# Smoke :  sbatch --array=0 mri_pipeline/baseline_vit/04i_finetune_ViT_tiny_scratch_submit_csd3.sh
# Already-completed runs (metrics.json present) are auto-skipped.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Offline wandb -- compute nodes have no internet; sync after the sweep.
export WANDB_MODE=offline

# -- Hardcoded paths ----------------------------------------------------------
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_tiny_baseline"

# -- Session selection (must match what 03 produced) --------------------------
PY_SESSION_FLAGS="--long all"

# -- Combination lookup -------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)
N_TASKS=${#TASKS[@]}    # 5
N_SEEDS=${#SEEDS[@]}    # 3

ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));    ID=$(( ID / N_SEEDS ))
TASK_IDX=$(( ID % N_TASKS ))

TASK="${TASKS[$TASK_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"
STRATEGY="scratch"
VIT_SIZE="tiny"
MODEL_SLUG="ViT_T_scratch"

RUN_DIR="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"
CKPT="${RUN_DIR}/last_checkpoint.pt"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-Tiny from-scratch ablation"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Array ID    : $SLURM_ARRAY_TASK_ID  (of 0-$((N_TASKS*N_SEEDS-1)))"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Task        : ${TASK}"
echo "  Seed        : ${SEED}"
echo "  Size        : ${VIT_SIZE}  (~5.5M params; ViT-B is ~86M)"
echo "  Strategy    : ${STRATEGY} (random init, no MAE checkpoint)"
echo "  Recipe      : lr=5e-4 | wd=5e-2 | ls=0.1 | aug=plus_original K=1"
echo "  Output dir  : ${RUN_DIR}"
echo "  wandb       : offline -> project vit_tiny_scratch"
echo "  Started     : $(date)"
echo "============================================================"

# -- Pre-flight checks --------------------------------------------------------
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs directory not found: ${VIT_INPUTS_DIR}"
    echo "        Run 03_prepare_ViT_submit_csd3.sh first."
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
if [ -f "${CKPT}" ]; then
    echo "  last_checkpoint.pt found -- python will AUTO-RESUME this combo."
fi

# NOTE: --pretrained_ckpt intentionally omitted (scratch = random init).
#       --vit_size tiny + --strategy scratch is the only valid combo for Tiny
#       (the trainer enforces this at parse_args).
python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --vit_size            "${VIT_SIZE}" \
    ${PY_SESSION_FLAGS} \
    --lr                  5e-4 \
    --weight_decay        5e-2 \
    --label_smoothing     0.1 \
    --augment             plus_original \
    --aug_copies          1 \
    --epochs              100 \
    --warmup_epochs       10 \
    --patience            15 \
    --llrd_gamma          1.0 \
    --grad_accum_steps    8 \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --wandb \
    --wandb_project       vit_tiny_scratch \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${RUN_DIR}/"
if [ ! -f "${METRICS}" ] && [ -f "${CKPT}" ]; then
    echo "  NOTE: no metrics.json yet (likely hit 12 h cap). Resubmit to AUTO-RESUME."
fi
echo "============================================================"

exit ${EXIT_CODE}
