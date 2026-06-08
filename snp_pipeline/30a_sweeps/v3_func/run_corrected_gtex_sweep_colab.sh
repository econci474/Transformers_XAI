#!/usr/bin/env bash
# Corrected-GTEx func-integration sweep — Colab A100 version.
#
# 81 runs = 9 modes × 3 families × 3 seeds. Each run logged to wandb
# project=snp_gtex_corrected_func_sweep, group=strict_v2_6mod_gpu (separate
# from any CPU runs).
#
# Required env vars (set before running):
#   SCRIPTS_DIR  — folder containing 30v3_train_diff_attention_func.py
#                  29b_load_functional_features.py, fm_diff_lib.py
#   AG_PARENT    — corrected_gtex/fm_embeddings_short_seq_1kb_with_alphagenome/
#   BASE         — diff_attn_drive_upload/ (has diff_emb + dosage + beta)
#   SPLITS_ROOT  — clinical splits root (per-seed train/val/test CSVs)
#   OUT_ROOT     — where to write per-run output dirs
#   WANDB_PROJECT (optional) — defaults to snp_gtex_corrected_func_sweep
#   WANDB_GROUP   (optional) — defaults to strict_v2_6mod_gpu

set -e

WANDB_PROJECT="${WANDB_PROJECT:-snp_gtex_corrected_func_sweep}"
WANDB_GROUP="${WANDB_GROUP:-strict_v2_6mod_gpu}"

MODES=(
    none
    per_snp_zscored_25
    per_snp_summaries_12
    attn_bias_single
    attn_bias_per_modality
    value_mult
    attn_bias_per_modality_plus_value_mult
    attn_bias_per_modality_plus_per_snp_summaries_12
    attn_bias_single_plus_func_imp_abs
)

total=$(( ${#MODES[@]} * 3 * 3 ))
i=0
t_start=$(date +%s)

for seed in 0 1 2; do
  for mode in "${MODES[@]}"; do
    # BMFM-SNP 2L MLP @ 101bp
    i=$((i+1))
    elapsed=$(( $(date +%s) - t_start ))
    printf '\n[%d/%d] BMFM-SNP 2L MLP  mode=%s seed=%d  (elapsed %ds)\n' \
        "$i" "$total" "$mode" "$seed" "$elapsed"
    python "$SCRIPTS_DIR/30v3_train_diff_attention_func.py" \
        --base "$BASE" --splits-root "$SPLITS_ROOT" \
        --ag-parent "$AG_PARENT" --output-root "$OUT_ROOT" \
        --seq-length 101bp --model bmfm_snp \
        --head mlp2 --aggregation chrom_hier \
        --mlp-width 256 --mlp-dropout 0.1 --lr 0.001 \
        --func-integration-mode "$mode" --seed "$seed" \
        --wandb-project "$WANDB_PROJECT" --wandb-group "$WANDB_GROUP"

    # NTv2 2L MLP @ 101bp
    i=$((i+1))
    elapsed=$(( $(date +%s) - t_start ))
    printf '\n[%d/%d] NTv2 2L MLP  mode=%s seed=%d  (elapsed %ds)\n' \
        "$i" "$total" "$mode" "$seed" "$elapsed"
    python "$SCRIPTS_DIR/30v3_train_diff_attention_func.py" \
        --base "$BASE" --splits-root "$SPLITS_ROOT" \
        --ag-parent "$AG_PARENT" --output-root "$OUT_ROOT" \
        --seq-length 101bp --model ntv2 \
        --head mlp2 --aggregation chrom_hier \
        --mlp-width 512 --mlp-dropout 0.3 --lr 0.0003 \
        --func-integration-mode "$mode" --seed "$seed" \
        --wandb-project "$WANDB_PROJECT" --wandb-group "$WANDB_GROUP"

    # BMFM-SNP XGB @ 1001bp
    i=$((i+1))
    elapsed=$(( $(date +%s) - t_start ))
    printf '\n[%d/%d] BMFM-SNP XGB  mode=%s seed=%d  (elapsed %ds)\n' \
        "$i" "$total" "$mode" "$seed" "$elapsed"
    python "$SCRIPTS_DIR/30v3_train_diff_attention_func.py" \
        --base "$BASE" --splits-root "$SPLITS_ROOT" \
        --ag-parent "$AG_PARENT" --output-root "$OUT_ROOT" \
        --seq-length 1001bp --model bmfm_snp \
        --head xgb --aggregation global_attn \
        --xgb-lr 0.05 --xgb-max-depth 4 --xgb-n-estimators 20 \
        --func-integration-mode "$mode" --seed "$seed" \
        --wandb-project "$WANDB_PROJECT" --wandb-group "$WANDB_GROUP"
  done
done

total_time=$(( $(date +%s) - t_start ))
printf '\n=== Sweep complete: %d runs in %dm %ds ===\n' \
    "$i" "$((total_time / 60))" "$((total_time % 60))"
