#!/bin/bash
# =============================================================================
# 03_encoder_submit_dept.sh
# =============================================================================
# Runs all combinations sequentially on the department cluster (single L4 GPU).
# 4 models × 4 tasks × 3 seeds × 2 strategies = 96 runs
#
# Estimated time on L4 (full FT base ~20 min, large ~55 min per run):
#   ~8–10 hours total for all 96 combinations
#
# Usage:
#   conda activate clinical
#   bash 03_encoder_submit_dept.sh
#
# Optional: skip already-completed runs (checks for metrics.json)
#   bash 03_encoder_submit_dept.sh --skip-done
# =============================================================================

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKIP_DONE=false
[[ "$1" == "--skip-done" ]] && SKIP_DONE=true

DATA_DIR="D:/ADNI_BIDS_project/derivatives/clinical/verbose/baseline"
OUT_DIR="${SCRIPT_DIR}/outputs"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Optional: point to RDS/network storage if HF home quota is tight
# export HF_HOME="/path/to/hf_cache"

MODEL_IDS=(
    "answerdotai/ModernBERT-base"
    "answerdotai/ModernBERT-large"
    "thomas-sounack/BioClinical-ModernBERT-base"
    "thomas-sounack/BioClinical-ModernBERT-large"
)
TASKS=("T1_binary" "T2_multiclass" "T3a_conv3y" "T3b_conv5y")
SEEDS=(0 1 2)
# strategy: "" = full fine-tune, "--freeze_backbone" = frozen
declare -A STRATEGY_FLAG=( ["full_ft"]="" ["frozen"]="--freeze_backbone" )
STRATEGIES=("full_ft" "frozen")

TOTAL=$(( ${#MODEL_IDS[@]} * ${#TASKS[@]} * ${#SEEDS[@]} * ${#STRATEGIES[@]} ))
COUNT=0
SKIPPED=0

echo "========================================================"
echo "  Encoder LLM sweep — department cluster (L4)"
echo "  Total runs: ${TOTAL}"
echo "========================================================"

for MODEL_ID in "${MODEL_IDS[@]}"; do
    MODEL_SLUG="${MODEL_ID##*/}"
    for TASK in "${TASKS[@]}"; do
        for SEED in "${SEEDS[@]}"; do
            for STRATEGY in "${STRATEGIES[@]}"; do
                COUNT=$((COUNT + 1))
                LOG="${LOG_DIR}/${MODEL_SLUG}_${TASK}_seed${SEED}_${STRATEGY}.log"
                METRICS="${OUT_DIR}/${MODEL_SLUG}/${TASK}/seed_${SEED}/${STRATEGY}/metrics.json"

                echo ""
                echo "  [${COUNT}/${TOTAL}] ${MODEL_SLUG} | ${TASK} | seed=${SEED} | ${STRATEGY}"

                if ${SKIP_DONE} && [ -f "${METRICS}" ]; then
                    echo "  → Skipping (metrics.json exists)"
                    SKIPPED=$((SKIPPED + 1))
                    continue
                fi

                python "${SCRIPT_DIR}/03_encoder_finetune.py" \
                    --model_id  "${MODEL_ID}" \
                    --task      "${TASK}" \
                    --seed      "${SEED}" \
                    --data_dir  "${DATA_DIR}" \
                    --out_dir   "${OUT_DIR}" \
                    ${STRATEGY_FLAG[${STRATEGY}]} \
                    2>&1 | tee "${LOG}"

                EXIT_CODE=${PIPESTATUS[0]}
                if [ ${EXIT_CODE} -ne 0 ]; then
                    echo "  ERROR: run failed (exit ${EXIT_CODE}) — see ${LOG}"
                fi
            done
        done
    done
done

echo ""
echo "========================================================"
echo "  All runs complete. Skipped: ${SKIPPED}/${TOTAL}"
echo "  Logs: ${LOG_DIR}"
echo "  Run 04_encoder_aggregate.py to produce summary table."
echo "========================================================"
