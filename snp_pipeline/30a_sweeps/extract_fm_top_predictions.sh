#!/usr/bin/env bash
# Run the v3 diff-attention top model 3× (seeds 0/1/2) with the funcon-sweep
# winning config: bmfm_snp / per_modality / mlp2 / 101bp / chrom_hier /
# w=256 / d=0.3 / lr=1e-3. With --save-predictions, each run writes a
# predictions.tsv with per-patient val + test scores.
#
# Run on Colab (or wherever Drive is mounted):
#   bash <repo>/snp_pipeline/30a_sweeps/extract_fm_top_predictions.sh
#
# Outputs at outputs/diff_attn_v3/fm_top_predictions/mlp2/101bp_bmfm_snp_attn_bias_per_modality_chrom_hier_s{0,1,2}/predictions.tsv
set -euo pipefail

PROJ="Transformers_XAI-snp_pipeline"
BASE="${BASE:-/content/drive/MyDrive/ADNI_SNP/diff_attn_v2_upload}"
TRAINER="${TRAINER:-snp_pipeline/30v3_train_diff_attention_func.py}"
OUT_ROOT="${OUT_ROOT:-outputs/diff_attn_v3/fm_top_predictions}"

for S in 0 1 2; do
    echo "-- seed=$S --"
    python "$TRAINER" \
        --base "$BASE" \
        --splits-root /nonexistent/force-fallback \
        --wandb-project "$PROJ" \
        --head mlp2 \
        --aggregation chrom_hier \
        --seq-length 101bp \
        --func-integration-mode attn_bias_per_modality \
        --model bmfm_snp \
        --mlp-width 256 \
        --mlp-dropout 0.3 \
        --lr 0.001 \
        --seed "$S" \
        --save-predictions \
        --output-root "$OUT_ROOT"
done

echo ""
echo "Done. Predictions at $OUT_ROOT/mlp2/*/predictions.tsv"
echo "Concatenate the 3 seed TSVs and upload back to local for the quantile plot."
