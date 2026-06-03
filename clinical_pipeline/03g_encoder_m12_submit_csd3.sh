#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_m12
# Logs land UNDER the m12 output directory. Slurm does NOT expand shell vars here, so this
# absolute path must match OUT_DIR/logs below. Create it before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m12/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m12/logs/enc_m12_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m12/logs/enc_m12_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-17%4
# =============================================================================
# 03g_encoder_m12_submit_csd3.sh
# =============================================================================
# Encoder-only LLM fine-tuning at the 1-YEAR (m12) CROSS-SECTIONAL visit, for the
# late-fusion arm (clinical + the m12-dominant MRI). INPUT = the patient's m12 summary,
# diagnosis LABEL = the concurrent m12 dx (Label_visit_diag). Kept SEPARATE from the
# 240-run 03d grid and the T4 03f job:
#   - its OWN split tree    -> verbose/baseline_m12   (built by 01h, same partition as baseline)
#   - its OWN output root    -> encoder_outputs_no_cdr_post_exclusion_m12
#   - its OWN wandb project  -> clinical-encoder-m12
#
# This is NOT a full grid — it runs only the per-task model picks from
# integration/*/notes.txt, each x 3 seeds (6 combos x 3 seeds = 18 runs, throttled %4):
#   T2_m12_multiclass : BioClinical-ModernBERT-large  full_ft
#   T1_m12_binary     : BioClinical-ModernBERT-base   full_ft
#   T1b_m12_cnmci_ad  : ModernBERT-base               frozen + full_ft
#   T1d_pmci_smci     : ModernBERT-large full_ft + ModernBERT-base full_ft   (m12 data_dir; conversion_group label)
#
# Prereqs:
#   - Build + upload verbose/baseline_m12 (run 01h locally, rsync to RDS).
#   - git pull on the HPC so 03_encoder_finetune.py has the *_m12 TASK_CONFIG entries.
#   - mkdir -p ${OUT_DIR}/logs
# Submit (run from ${OUT_DIR} so any relative paths land there):
#   sbatch clinical_pipeline/03g_encoder_m12_submit_csd3.sh
# Smoke one combo first (combo 0 / seed 0 = BioClinical-large T2_m12 full_ft):
#   sbatch --array=0 clinical_pipeline/03g_encoder_m12_submit_csd3.sh
#
# After runs finish, sync the offline W&B runs + build the m12 embeddings manifest:
#   wandb sync --sync-all "${OUT_DIR}/wandb_offline"
#   python clinical_pipeline/03e_build_embeddings_manifest.py \
#       --out_dir "${OUT_DIR}" --data_dir "${DATA_DIR}"
# =============================================================================

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate clinical

# HF cache on RDS (avoids home-dir quota)
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/clinical_pipeline"

# m12 cross-sectional VERBOSE splits (seed_N/{train,val,test}.csv with m12 Generated_Text,
# concurrent Label_visit_diag, carried conversion_group), built by 01h from verbose/baseline.
# 03_encoder_finetune.py appends seed_{N} itself, so DATA_DIR is the dir CONTAINING seed_*.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline_m12"

# Separate m12 output root (isolated from the 240-run sweep and T4)
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m12"

# W&B offline — collect under one dir for a single `wandb sync --sync-all`
export WANDB_MODE="offline"
export WANDB_DIR="${OUT_DIR}/wandb_offline"
mkdir -p "${WANDB_DIR}"
WANDB_PROJECT="clinical-encoder-m12"

# ── Combination lookup (notes' per-task model picks) ──────────────────────────
# Each entry: "MODEL_ID|TASK|STRATEGY". 6 combos x 3 seeds = 18 runs.
COMBOS=(
    "thomas-sounack/BioClinical-ModernBERT-large|T2_m12_multiclass|full_ft"
    "thomas-sounack/BioClinical-ModernBERT-base|T1_m12_binary|full_ft"
    "answerdotai/ModernBERT-base|T1b_m12_cnmci_ad|frozen"
    "answerdotai/ModernBERT-base|T1b_m12_cnmci_ad|full_ft"
    "answerdotai/ModernBERT-large|T1d_pmci_smci|full_ft"
    "answerdotai/ModernBERT-base|T1d_pmci_smci|full_ft"
)
SEEDS=(0 1 2)

N_COMBOS=${#COMBOS[@]}    # 6
N_SEEDS=${#SEEDS[@]}      # 3

# Decode SLURM_ARRAY_TASK_ID -> (combo, seed)   (seed-minor)
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % N_SEEDS ));     ID=$(( ID / N_SEEDS ))
COMBO_IDX=$(( ID % N_COMBOS ))

IFS='|' read -r MODEL_ID TASK STRATEGY <<< "${COMBOS[$COMBO_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

FREEZE_FLAG=""
[[ "${STRATEGY}" == "frozen" ]] && FREEZE_FLAG="--freeze_backbone"

MODEL_SLUG="${MODEL_ID##*/}"
METRICS="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}/metrics.json"

echo "========================================================"
echo "  Encoder-only LLM — m12 cross-sectional (1-year) fusion arm"
echo "  SLURM array ID : ${SLURM_ARRAY_TASK_ID}"
echo "  Model          : ${MODEL_ID}"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY}"
echo "  Data dir       : ${DATA_DIR}"
echo "  Output dir     : ${OUT_DIR}"
echo "  W&B (offline)  : ${WANDB_PROJECT}  (dir: ${WANDB_DIR})"
echo "  GPU            : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "========================================================"

# Skip if already done
if [ -f "${METRICS}" ]; then
    echo "  metrics.json exists — skipping."
    exit 0
fi

python "${SCRIPT_DIR}/03_encoder_finetune.py" \
    --model_id       "${MODEL_ID}" \
    --task           "${TASK}" \
    --seed           "${SEED}" \
    --data_dir       "${DATA_DIR}" \
    --out_dir        "${OUT_DIR}" \
    --hf_cache       "${HF_HOME}" \
    --wandb \
    --wandb_project  "${WANDB_PROJECT}" \
    ${FREEZE_FLAG}

echo "  Finished."
