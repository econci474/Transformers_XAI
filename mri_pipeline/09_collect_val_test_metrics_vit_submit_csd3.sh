#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=val_test_vit
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/slurm_logs/val_test_vit_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/slurm_logs/val_test_vit_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# =============================================================================
# 09_collect_val_test_metrics_vit_submit_csd3.sh
#   Recompute val + test metrics from saved best_model.pt for the ViT models
#   (ViT-MAE75 + ViT-Base/scratch) that only logged val_bacc, patching a
#   `val_metrics` block into each metrics.json. NO retraining.
#
#   Provenance: the driver discovers runs via 05.MODEL_TREES, so it recomputes
#   from the exact run each table row displays:
#     - ViT-MAE frozen/random  -> vit_outputs_debug/...frozen
#     - ViT-MAE full_ft + frozen/{none,plus_original} -> vit_outputs_hi_lr/...
#     - ViT-Base (scratch)     -> vit_baseline/...
#
# Prereq: the ViT run trees (vit_outputs_debug/, vit_outputs_hi_lr/,
#   vit_baseline/) with best_model.pt must live under ${DERIVS_ROOT} below,
#   and the vit_inputs volumes + MAE pretrained ckpt must be reachable.
#   Verify first, e.g.:
#     ls ${DERIVS_ROOT}/vit_outputs_debug/ViT_B_mae75/T2_multiclass/seed_0/frozen/best_model.pt
#
# Eval-only is fast (well under SL3 12 h). Idempotent via --only-missing.
#
# Run:
#   cd /home/ec474/rds/hpc-work/ADNI_MRI       # so relative SBATCH logs land here
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/09_collect_val_test_metrics_vit_submit_csd3.sh
# Then sync patched metrics.json back to D: and re-run 06 -> 06c -> 07 locally.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
# Env-independent: target the named env explicitly (do NOT `conda activate`).
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps (torch/monai/sklearn/pandas)"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline"
mkdir -p "${ROOT}/ADNI_MRI/slurm_logs"

DERIVS_ROOT="${ROOT}/ADNI_MRI"
VIT_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/vit_inputs"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
VIT_CKPT="${ROOT}/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"

echo "============================================================"
echo "  ViT val/test recompute — $(date)"
echo "  derivs-root: ${DERIVS_ROOT}"
echo "============================================================"

conda run -n mri python "${SCRIPT_DIR}/09_collect_val_test_metrics.py" \
    --derivs-root    "${DERIVS_ROOT}" \
    --vit-inputs     "${VIT_INPUTS}" \
    --vit-ckpt       "${VIT_CKPT}" \
    --data-dir       "${DATA_DIR}" \
    --matched-labels "${MATCHED}" \
    --models vit_mae vit_scratch \
    --num-workers 4 --only-missing
rc=$?

echo "  Finished: $(date)  exit=${rc}"
exit ${rc}
