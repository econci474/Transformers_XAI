#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=enc_m06
# Logs land UNDER the m06 output directory. Slurm does NOT expand shell vars here, so this
# absolute path must match OUT_DIR/logs below. Create it before submitting:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m06/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m06/logs/enc_m06_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m06/logs/enc_m06_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-17%4
# =============================================================================
# 03h_encoder_m06_submit_csd3.sh
# =============================================================================
# Encoder-only LLM fine-tuning at the 6-MONTH (m06) CROSS-SECTIONAL visit. Together with the
# baseline (bl) runs and the m12 runs (03g), this provides the per-visit clinical probabilities
# bl + m06 + m12 for the TEMPORAL clinical elastic net (EN over a patient's visits up to year 1),
# whose output is then averaged with MRI@m12. Clinical assessments are denser than MRI, so the
# temporal EN exploits more clinical evidence than a single timepoint.
#
# Sibling of 03g (m12); identical structure, m06 paths/tasks:
#   - split tree    -> verbose/baseline_m06   (built by 01h --visit m06)
#   - output root    -> encoder_outputs_no_cdr_post_exclusion_m06
#   - wandb project  -> clinical-encoder-m06
#
# 6 notes'-pick combos x 3 seeds = 18 runs (throttled %4):
#   T2_m06_multiclass : BioClinical-ModernBERT-large  full_ft
#   T1_m06_binary     : BioClinical-ModernBERT-base   full_ft
#   T1b_m06_cnmci_ad  : ModernBERT-base               frozen + full_ft
#   T1d_pmci_smci     : ModernBERT-large full_ft + ModernBERT-base full_ft   (m06 data_dir; conversion_group)
#
# Prereqs:
#   - Build + upload verbose/baseline_m06 (run 01h --visit m06 locally, scp/rsync to RDS).
#   - git pull on the HPC so 03_encoder_finetune.py has the *_m06 TASK_CONFIG entries.
#   - mkdir -p ${OUT_DIR}/logs
# Submit (run from ${OUT_DIR} so any relative paths land there):
#   sbatch clinical_pipeline/03h_encoder_m06_submit_csd3.sh
# Smoke one combo first (combo 0 / seed 0 = BioClinical-large T2_m06 full_ft):
#   sbatch --array=0 clinical_pipeline/03h_encoder_m06_submit_csd3.sh
#
# After runs finish, sync the offline W&B runs + build the m06 embeddings manifest:
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

# m06 cross-sectional VERBOSE splits (seed_N/{train,val,test}.csv with m06 Generated_Text,
# concurrent Label_visit_diag, carried conversion_group), built by 01h --visit m06.
# 03_encoder_finetune.py appends seed_{N} itself, so DATA_DIR is the dir CONTAINING seed_*.
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline_m06"

# Separate m06 output root (isolated from the 240-run sweep, T4, and m12)
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_outputs_no_cdr_post_exclusion_m06"

# W&B offline — collect under one dir for a single `wandb sync --sync-all`
export WANDB_MODE="offline"
export WANDB_DIR="${OUT_DIR}/wandb_offline"
mkdir -p "${WANDB_DIR}"
WANDB_PROJECT="clinical-encoder-m06"

# ── Combination lookup (notes' per-task model picks) ──────────────────────────
# Each entry: "MODEL_ID|TASK|STRATEGY". 6 combos x 3 seeds = 18 runs.
COMBOS=(
    "thomas-sounack/BioClinical-ModernBERT-large|T2_m06_multiclass|full_ft"
    "thomas-sounack/BioClinical-ModernBERT-base|T1_m06_binary|full_ft"
    "answerdotai/ModernBERT-base|T1b_m06_cnmci_ad|frozen"
    "answerdotai/ModernBERT-base|T1b_m06_cnmci_ad|full_ft"
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
echo "  Encoder-only LLM — m06 cross-sectional (6-month) fusion arm"
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
