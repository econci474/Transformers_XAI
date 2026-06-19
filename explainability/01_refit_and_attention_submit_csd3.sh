#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=explain_attn
# Logs land under the checkpoint OUT_DIR. Slurm does NOT expand vars here — must match OUT_DIR/logs.
# Create it before submitting:  mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_explain_checkpoints/logs
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_CL/encoder_explain_checkpoints/logs/explain_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_CL/encoder_explain_checkpoints/logs/explain_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# =============================================================================
# 01_refit_and_attention_submit_csd3.sh   (ONE job, ONE GPU)
# =============================================================================
# Explainability for the best LLM per task (T1d / T1e / T2): the post-exclusion sweep ran WITHOUT
# --save_model, so re-fit the 3 best (model, task) with weights saved, then run attention rollout.
# Both steps share a single GPU in this one job (never more than one GPU — user constraint).
#
# Best LLM per task (val balanced accuracy, full fine-tune):
#   T1d -> answerdotai/ModernBERT-large             (T1d_pmci_smci)
#   T1e -> answerdotai/ModernBERT-base              (T1e_scn_pcn)
#   T2  -> thomas-sounack/BioClinical-ModernBERT-large (T2_multiclass)
#
# Submit:
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_CL/encoder_explain_checkpoints/logs
#   sbatch explainability/01_refit_and_attention_submit_csd3.sh
#
# Account: COMPUTERLAB-SL3-GPU (ampere, 12h cap, free pool). The job is short (~1-2 h) so 12h is ample.
# (LIO-CHARM-SL2-GPU would also work at 36h but SL3-GPU is the intended account here.)
#
# Outputs: checkpoints under OUT_DIR (gitignored); attention_*.png under the repo's
# explainability/<task>/ on RDS -> commit/scp those back to the local repo.
# =============================================================================
set -euo pipefail

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
export HF_HOME="/home/ec474/rds/hpc-work/hf_cache"

REPO="/home/ec474/rds/hpc-work/Transformers_XAI"
FT="${REPO}/clinical_pipeline/03_encoder_finetune.py"
ATTN="${REPO}/explainability/02_attention_rollout.py"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/verbose/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_CL/encoder_explain_checkpoints"
SEED="${SEED:-0}"
mkdir -p "${OUT_DIR}/logs"

# env-independent: use `conda run -n clinical` (NOT activate — avoids inheriting the login env), with
# a fail-fast preflight (matches the repo's submit-script convention).
conda run -n clinical python -c "import torch, transformers, matplotlib; \
print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'tf', transformers.__version__)" || exit 1

echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "SEED=${SEED}  DATA_DIR=${DATA_DIR}  OUT_DIR=${OUT_DIR}"

# (model_id|task_key) — best LLM per task, full fine-tune
CONFIGS=(
  "answerdotai/ModernBERT-large|T1d_pmci_smci"
  "answerdotai/ModernBERT-base|T1e_scn_pcn"
  "thomas-sounack/BioClinical-ModernBERT-large|T2_multiclass"
)

# ── Step 1: re-fit the 3 best LLMs WITH weights (sequential on the one GPU) ──
for cfg in "${CONFIGS[@]}"; do
  IFS='|' read -r model_id task_key <<< "${cfg}"
  slug="${model_id##*/}"
  ckpt="${OUT_DIR}/${slug}/${task_key}/seed_${SEED}/full_ft/best_checkpoint"
  echo ""
  echo "=== re-fit ${task_key} <- ${model_id} (seed ${SEED}) ==="
  if [ -d "${ckpt}" ]; then
    echo "  [skip] checkpoint exists: ${ckpt}"
    continue
  fi
  conda run -n clinical python "${FT}" \
    --model_id "${model_id}" \
    --task     "${task_key}" \
    --seed     "${SEED}" \
    --data_dir "${DATA_DIR}" \
    --out_dir  "${OUT_DIR}" \
    --hf_cache "${HF_HOME}" \
    --save_model || exit 1
done

# ── Step 2: attention rollout overlays (same GPU, after the checkpoints exist) ──
echo ""
echo "=== attention rollout ==="
FIG_DIR="/home/ec474/rds/hpc-work/ADNI_CL/explainability"   # outputs off the repo; scp to D: after
conda run -n clinical python "${ATTN}" \
  --ckpt_root "${OUT_DIR}" \
  --data_dir  "${DATA_DIR}" \
  --out_dir   "${FIG_DIR}" \
  --seed      "${SEED}" || exit 1

echo ""
echo "Done. Attention figures: ${FIG_DIR}/<task>/"
echo "  -> scp back to D:\\ADNI_BIDS_project\\derivatives\\explainability\\<task>\\ (off the repo)."
