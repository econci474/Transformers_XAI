#!/usr/bin/env bash
# Drop-splice ablation on the top-2 MLP2 cells (per 2026-05-27 funcon sweep).
# 2 configs × 3 seeds = 6 runs. Logged to project Transformers_XAI-snp_pipeline
# (same project as the sweeps) with run names suffixed _no_splice.
#
# Run on Colab (or any GPU host) AFTER `git pull` so the new --drop-splice
# flag is present in 30v3_train_diff_attention_func.py:
#
#   cd /content/drive/MyDrive/ADNI_SNP/diff_attn_v2_upload
#   bash <repo>/snp_pipeline/30a_sweeps/drop_splice_top2.sh
#
# Time: ~6 × ~60s on A100 ≈ 6 minutes.
set -euo pipefail

PROJ="Transformers_XAI-snp_pipeline"
BASE="${BASE:-/content/drive/MyDrive/ADNI_SNP/diff_attn_v2_upload}"
TRAINER="${TRAINER:-snp_pipeline/30v3_train_diff_attention_func.py}"
OUT_ROOT="${OUT_ROOT:-outputs/diff_attn_v3/drop_splice_top2}"

run () {
    MODEL="$1"; WIDTH="$2"; DROPOUT="$3"; LR="$4"; SEED="$5"
    echo "── ${MODEL}  w=${WIDTH}  d=${DROPOUT}  lr=${LR}  seed=${SEED}  +no_splice ──"
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
        --drop-splice \
        --output-root "$OUT_ROOT"
}

# Config 1 — MLP2 rank-1: bmfm_snp / per_modality / w=256 / d=0.3 / lr=1e-3
for S in 0 1 2; do
    run bmfm_snp 256 0.3 0.001 "$S"
done

# Config 2 — MLP2 rank-2 (per-model winner ntv2): ntv2 / per_modality / w=512 / d=0.3 / lr=3e-4
for S in 0 1 2; do
    run ntv2 512 0.3 0.0003 "$S"
done

echo "Done. Check WandB project ${PROJ} for runs named ..._no_splice."
