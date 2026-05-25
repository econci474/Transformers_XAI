#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_extract
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/brainmvp/slurm_logs/bmvp_extract_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/brainmvp/slurm_logs/bmvp_extract_%j.err
#SBATCH -p ampere
#SBATCH --exclude=gpu-q-[15,64,68-75,81-82,86-90]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
# =============================================================================
# 03_extract_brainmvp_embeddings_submit_csd3.sh
# =============================================================================
# Single non-array GPU job. Forwards every preprocessed BrainMVP input
# scan (CenterSpatialCrop 96x96x64) through the frozen UniFormer-Small
# encoder, AdaptiveAvgPool3d on stage-4 features (512 channels), saves
# paired (parquet + .pt) keyed by (bids_sub, bids_ses).
#
# Expected wall-clock: ~10 min for ~3000 scans at batch=8 on A100.
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/brainmvp/slurm_logs
#   ls /home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt
#   ls /home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs/
#
# Submit: sbatch mri_pipeline/brain_mvp/03_extract_brainmvp_embeddings_submit_csd3.sh
# Resume: re-running auto-skips scans already in the parquet.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/brainmvp"

mkdir -p "${OUT_DIR}/slurm_logs"

echo "============================================================"
echo "  BrainMVP frozen-encoder embedding extraction"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Out dir     : ${OUT_DIR}"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${BRAINMVP_INPUTS_DIR}" ]; then
    echo "[ERROR] BrainMVP inputs dir not found: ${BRAINMVP_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi

python "${SCRIPT_DIR}/03_extract_brainmvp_embeddings.py" \
    --pretrained_ckpt        "${PRETRAINED_CKPT}" \
    --brainmvp_inputs_dir    "${BRAINMVP_INPUTS_DIR}" \
    --matched_labels_csv     "${MATCHED_LABELS_CSV}" \
    --out_dir                "${OUT_DIR}" \
    --batch_size             8 \
    --num_workers            "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${OUT_DIR}/{brainmvp_pooled.parquet,brainmvp_pooled.pt,manifest.csv}"
echo "============================================================"
exit ${EXIT_CODE}
