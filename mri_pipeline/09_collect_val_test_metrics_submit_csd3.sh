#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=val_test_recompute
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/slurm_logs/val_test_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/slurm_logs/val_test_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# =============================================================================
# 09_collect_val_test_metrics_submit_csd3.sh
#   Recompute val + test metrics from saved best_model.pt for BrainMVP /
#   AG-MS3D / 3D-CNN (the supervised trainers that only logged val_bacc),
#   patching a `val_metrics` block into each metrics.json. NO retraining.
#
# Prereqs on HPC (the input volumes already live here; only the checkpoints
# were deleted, so re-push the slim best_model.pt from D: first):
#   rsync the supervised best_model.pt back under
#     /home/ec474/rds/hpc-work/ADNI_MRI/{brainmvp_debug,cnn3d_outputs,agms3d_outputs}/
#   (see the local push commands in the session notes).
#
# Single job (the driver loops all runs sequentially; eval-only is fast — well
# under the SL3 12 h cap). Idempotent via --only-missing.
#
# Run:
#   cd /home/ec474/rds/hpc-work/ADNI_MRI       # so relative SBATCH logs land here
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/09_collect_val_test_metrics_submit_csd3.sh
# Then sync the patched metrics.json back to D: and re-run 06 -> 06c locally.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
# Env-independent: target the named env explicitly (do NOT `conda activate`,
# which stacks on the inherited login-node env). Preflight deps, fail fast.
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps (torch/monai/sklearn/pandas)"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline"
mkdir -p "${ROOT}/ADNI_MRI/slurm_logs"

DERIVS_ROOT="${ROOT}/ADNI_MRI"
BRAINMVP_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
CNN_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/cnn_inputs"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
BRAINMVP_CKPT="${ROOT}/ViT_pretrained/BrainMVP_uniformer.pt"

echo "============================================================"
echo "  val/test recompute — $(date)"
echo "  derivs-root: ${DERIVS_ROOT}"
echo "============================================================"

conda run -n mri python "${SCRIPT_DIR}/09_collect_val_test_metrics.py" \
    --derivs-root   "${DERIVS_ROOT}" \
    --brainmvp-inputs "${BRAINMVP_INPUTS}" \
    --cnn-inputs    "${CNN_INPUTS}" \
    --data-dir      "${DATA_DIR}" \
    --matched-labels "${MATCHED}" \
    --brainmvp-ckpt "${BRAINMVP_CKPT}" \
    --models brainmvp cnn3d agms3d \
    --num-workers 4 --only-missing
rc=$?

echo "  Finished: $(date)  exit=${rc}"
exit ${rc}
