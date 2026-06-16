#!/usr/bin/env bash
# 01_refit_best_llms.sh  (department L4, env: clinical) — SINGLE GPU, STRICTLY SEQUENTIAL
# =====================================================================================
# The post-exclusion encoder sweep was run WITHOUT --save_model, so no LLM weights were kept (only
# metrics.json + embeddings). Attention rollout (02_attention_rollout.py) needs the weights, so we re-fit
# the BEST (model, task) per the leaderboard for ONE seed, with --save_model, to recreate the checkpoints.
#
# Best LLM per task (encoder_only_post_exclusions/combined_model_table_val.csv, val balanced accuracy,
# full fine-tune):
#   T1d -> answerdotai/ModernBERT-large             (0.734)
#   T1e -> answerdotai/ModernBERT-base              (0.868)
#   T2  -> thomas-sounack/BioClinical-ModernBERT-large (0.839)
#
# Runs ONE GPU job at a time (no '&'), CUDA_VISIBLE_DEVICES=0 — never more than one GPU (user constraint).
# Env-independent: uses `conda run -n clinical` so it does not depend on the login shell's active env.
# Skips a config whose best_checkpoint/ already exists. ~1 h total for the 3 runs.
#
# Override the data/out/repo locations for the machine you run on (defaults below are the local Windows
# paths; on the dept L4 set DATA_DIR / OUT_DIR / REPO to the Linux paths where the post-exclusion verbose
# baseline splits live, e.g. /home/ec474/.../no_cdr_stratified_post_exclusion/verbose/baseline).
#
# Usage:
#   bash explainability/01_refit_best_llms.sh            # seed 0 (default)
#   SEED=1 bash explainability/01_refit_best_llms.sh     # a different seed
set -euo pipefail

REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_DIR="${DATA_DIR:-D:/ADNI_BIDS_project/derivatives/clinical/no_cdr_stratified_post_exclusion/verbose/baseline}"
OUT_DIR="${OUT_DIR:-D:/ADNI_BIDS_project/derivatives/clinical/encoder_explain_checkpoints}"
SEED="${SEED:-0}"
FT="${REPO}/clinical_pipeline/03_encoder_finetune.py"
export CUDA_VISIBLE_DEVICES=0          # single GPU only

echo "REPO=${REPO}"
echo "DATA_DIR=${DATA_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "SEED=${SEED}"
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
  echo "ERROR: ${DATA_DIR}/seed_${SEED} not found. Set DATA_DIR to the post-exclusion verbose baseline dir on THIS machine." >&2
  exit 1
fi

# preflight: clinical env imports (fail fast, env-independent)
conda run -n clinical python -c "import torch, transformers; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'tf', transformers.__version__)" || exit 1

# (subdir  model_id  task_key)
CONFIGS=(
  "ModernBERT-large|answerdotai/ModernBERT-large|T1d_pmci_smci"
  "ModernBERT-base|answerdotai/ModernBERT-base|T1e_scn_pcn"
  "BioClinical-ModernBERT-large|thomas-sounack/BioClinical-ModernBERT-large|T2_multiclass"
)

for cfg in "${CONFIGS[@]}"; do
  IFS='|' read -r slug model_id task_key <<< "${cfg}"
  ckpt="${OUT_DIR}/${slug}/${task_key}/seed_${SEED}/full_ft/best_checkpoint"
  echo ""
  echo "============================================================"
  echo "  ${task_key}  <-  ${model_id}  (seed ${SEED})"
  echo "============================================================"
  if [ -d "${ckpt}" ]; then
    echo "  [skip] checkpoint exists: ${ckpt}"
    continue
  fi
  conda run -n clinical python "${FT}" \
    --model_id "${model_id}" \
    --task "${task_key}" \
    --seed "${SEED}" \
    --data_dir "${DATA_DIR}" \
    --out_dir "${OUT_DIR}" \
    --save_model || exit 1
done

echo ""
echo "Done. Checkpoints under ${OUT_DIR}/<model>/<task>/seed_${SEED}/full_ft/best_checkpoint/"
echo "Next: conda run -n clinical python explainability/02_attention_rollout.py --ckpt_root \"${OUT_DIR}\" --seed ${SEED}"
