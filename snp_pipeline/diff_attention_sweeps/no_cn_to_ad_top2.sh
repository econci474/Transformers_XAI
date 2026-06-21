#!/usr/bin/env bash
# CN_to_AD ablation on the top-2 MLP2 cells (per 2026-05-27 funcon sweep).
# Drops the CN_to_AD subgroup from POS class. Seed 1 val originally has
# 0 CN_to_AD subjects — this run isolates whether the cross-seed variance
# is driven by that subgroup imbalance.
# 2 configs x 3 seeds = 6 runs. Logged to Transformers_XAI-snp_pipeline.
set -euo pipefail

PROJ="Transformers_XAI-snp_pipeline"
BASE="${BASE:-/content/drive/MyDrive/ADNI_SNP/diff_attn_v2_upload}"
TRAINER="${TRAINER:-snp_pipeline/30v3_train_diff_attention_func.py}"
OUT_ROOT="${OUT_ROOT:-outputs/diff_attn_v3/no_cn_to_ad_top2}"

run () {
    MODEL="$1"; WIDTH="$2"; DROPOUT="$3"; LR="$4"; SEED="$5"
    echo "-- ${MODEL}  w=${WIDTH}  d=${DROPOUT}  lr=${LR}  seed=${SEED}  +no_cn_to_ad --"
    python "$TRAINER" \
        --base "$BASE" \
        --splits-root /nonexistent/force-fallback \
        --wandb-project "$PROJ" \
        --head mlp2 \
        --aggregation chrom_hier \
        --seq-length 101bp \
        --func-integration-mode attn_bias_per_modality \
        --model "$MODEL" \
        --mlp-width "$WIDTH" \
        --mlp-dropout "$DROPOUT" \
        --lr "$LR" \
        --seed "$SEED" \
        --exclude-cn-to-ad \
        --output-root "$OUT_ROOT"
}

for S in 0 1 2; do run bmfm_snp 256 0.3 0.001  "$S"; done
for S in 0 1 2; do run ntv2     512 0.3 0.0003 "$S"; done
