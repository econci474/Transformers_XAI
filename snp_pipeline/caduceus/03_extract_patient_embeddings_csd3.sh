#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=cad_emb
#SBATCH --output=logs/cad_emb_%A_%a.log
#SBATCH --error=logs/cad_emb_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=0-1
# =============================================================================
# 03_extract_patient_embeddings_csd3.sh
# =============================================================================
# Extract mean-pooled embeddings from frozen Caduceus for all ~180 patients.
#
# Caduceus uses single-nucleotide tokenisation: 8,192 bp → 8,192 tokens
# (well within the 131k context window). Hidden dim = 256.
#
# Array: 0 = ps (parameter sharing), 1 = ph (post-hoc conjoining)
#
# Usage:
#   sbatch snp_pipeline/caduceus/03_extract_patient_embeddings_csd3.sh
#   sbatch --array=0 snp_pipeline/caduceus/03_extract_patient_embeddings_csd3.sh  # ps only
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate caduceus

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline/caduceus"
SEQUENCES="/home/ec474/rds/hpc-work/ADNI_SNP/patient_seqs/all_patients.csv"
OUTDIR="/home/ec474/rds/hpc-work/ADNI_SNP/patient_embeddings"

VARIANTS=("ps" "ph")
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"
OUT_NPZ="${OUTDIR}/caduceus_${VARIANT}/embeddings.npz"

mkdir -p logs
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "  Caduceus — frozen embedding extraction"
echo "  Job ID       : $SLURM_JOB_ID"
echo "  Array ID     : $SLURM_ARRAY_TASK_ID"
echo "  Variant      : ${VARIANT}"
echo "  Node         : $SLURMD_NODENAME"
echo "  GPU          : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Sequences    : ${SEQUENCES}"
echo "  Output       : ${OUTDIR}/caduceus_${VARIANT}/"
echo "  Started      : $(date)"
echo "============================================================"

# Pre-flight
if [ ! -f "${SEQUENCES}" ]; then
    echo "[ERROR] Sequences not found: ${SEQUENCES}"
    exit 1
fi

# Skip if already done
if [ -f "${OUT_NPZ}" ]; then
    echo "  embeddings.npz already exists -- skipping."
    exit 0
fi

python "${SCRIPT_DIR}/03_extract_patient_embeddings.py" \
    --sequences  "${SEQUENCES}" \
    --outdir     "${OUTDIR}" \
    --variant    "${VARIANT}" \
    --device     cuda \
    --batch-size 4

EXIT_CODE=$?

echo "============================================================"
echo "  Finished     : $(date)"
echo "  Exit code    : ${EXIT_CODE}"
echo "============================================================"

exit ${EXIT_CODE}
