#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bdino_extract
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/slurm_logs/bdino_extract_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/slurm_logs/bdino_extract_%j.err
#SBATCH -p ampere
#SBATCH --exclude=gpu-q-[15,64,68-75,81-82,86-90]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
# =============================================================================
# 03_extract_braindino_embeddings_submit_csd3.sh
# =============================================================================
# Single non-array GPU job. Forwards every preprocessed BrainDINO input
# scan through the frozen encoder, mean-pools the 128 axial CLS tokens,
# saves a paired (parquet + .pt) keyed by (bids_sub, bids_ses).
#
# Expected wall-clock: ~30 min for ~3000 scans at batch=4 on A100.
#
# Pre-flight (run once):
#   mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino/slurm_logs
#   ls /home/ec474/rds/hpc-work/brain_dino_model.pth
#   ls /home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs/
#
# Submit: sbatch mri_pipeline/brain_dino/03_extract_braindino_embeddings_submit_csd3.sh
# Resume: re-running auto-skips scans already in the parquet.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_dino"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/brain_dino_model.pth"
BRAINDINO_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/braindino"

mkdir -p "${OUT_DIR}/slurm_logs"

echo "============================================================"
echo "  BrainDINO frozen-encoder embedding extraction"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Inputs      : ${BRAINDINO_INPUTS_DIR}"
echo "  Out dir     : ${OUT_DIR}"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${BRAINDINO_INPUTS_DIR}" ]; then
    echo "[ERROR] BrainDINO inputs dir not found: ${BRAINDINO_INPUTS_DIR}"
    echo "        Run 01_prepare_braindino_inputs.py first."; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi

python "${SCRIPT_DIR}/03_extract_braindino_embeddings.py" \
    --pretrained_ckpt        "${PRETRAINED_CKPT}" \
    --braindino_inputs_dir   "${BRAINDINO_INPUTS_DIR}" \
    --matched_labels_csv     "${MATCHED_LABELS_CSV}" \
    --out_dir                "${OUT_DIR}" \
    --batch_size             4 \
    --num_workers            "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${OUT_DIR}/{braindino_pooled.parquet,braindino_pooled.pt,manifest.csv}"
echo "============================================================"
exit ${EXIT_CODE}
