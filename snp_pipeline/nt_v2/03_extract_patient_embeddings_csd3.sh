#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=ntv2_emb
#SBATCH --output=logs/ntv2_emb_%j.log
#SBATCH --error=logs/ntv2_emb_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
# =============================================================================
# 03_extract_patient_embeddings_csd3.sh
# =============================================================================
# Extract mean-pooled embeddings from frozen NT v2 500M for all ~180 patients.
#
# NT v2 uses 6-mer tokenisation: 8,192 bp → ~1,366 tokens (within 2,048 limit).
# Hidden dim = 1024. Output: embeddings.npz with legacy format (no prefix).
#
# Usage:
#   sbatch snp_pipeline/nt_v2/03_extract_patient_embeddings_csd3.sh
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ntv2

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline/nt_v2"
SEQUENCES="/home/ec474/rds/hpc-work/ADNI_SNP/patient_seqs/all_patients.csv"
OUTDIR="/home/ec474/rds/hpc-work/ADNI_SNP/patient_embeddings/ntv2"
OUT_NPZ="${OUTDIR}/embeddings.npz"

mkdir -p logs
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "  NT v2 500M — frozen embedding extraction"
echo "  Job ID       : $SLURM_JOB_ID"
echo "  Node         : $SLURMD_NODENAME"
echo "  GPU          : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Sequences    : ${SEQUENCES}"
echo "  Output       : ${OUTDIR}"
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
    --sequences "${SEQUENCES}" \
    --outdir    "${OUTDIR}" \
    --device    cuda \
    --batch-size 8

EXIT_CODE=$?

echo "============================================================"
echo "  Finished     : $(date)"
echo "  Exit code    : ${EXIT_CODE}"
echo "============================================================"

exit ${EXIT_CODE}
