#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=cad131k
#SBATCH --output=logs/cad_131k_emb_%A_%a.log
#SBATCH --error=logs/cad_131k_emb_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
# =============================================================================
# 03b_extract_131k_patient_embeddings_csd3.sh
# =============================================================================
# Extract mean-pooled embeddings from frozen Caduceus using the WIDE 131k
# context windows (built by 02_build_caduceus_patient_window_sequences.py).
#
# These are ~131 kb windows (vs 8 kb standard). Fewer windows per patient
# but much richer genomic context. Sequences will be ~131k tokens each.
#
# Array: 0 = ps (parameter sharing), 1 = ph (post-hoc conjoining)
#
# NOTE: batch-size=1 due to the much longer sequences (~131k tokens).
#       Memory requirement is higher (64G RAM + GPU).
#
# Usage:
#   sbatch snp_pipeline/caduceus/03b_extract_131k_patient_embeddings_csd3.sh
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate caduceus

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline/caduceus"
SEQUENCES="/home/ec474/rds/hpc-work/ADNI_SNP/caduceus_inputs/patient_window_sequences/sequences/all_patients.csv"
OUTDIR="/home/ec474/rds/hpc-work/ADNI_SNP/caduceus_131k_embeddings"

VARIANTS=("ps" "ph")
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"
OUT_NPZ="${OUTDIR}/caduceus_${VARIANT}/embeddings.npz"

mkdir -p logs
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "  Caduceus 131k — frozen embedding extraction (WIDE windows)"
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
    echo "  Build them first with:"
    echo "    python snp_pipeline/caduceus/02_build_caduceus_patient_window_sequences.py"
    exit 1
fi

# Skip if already done
if [ -f "${OUT_NPZ}" ]; then
    echo "  embeddings.npz already exists -- skipping."
    exit 0
fi

# Batch size=1 for 131k token sequences (GPU memory constraint)
python "${SCRIPT_DIR}/03_extract_patient_embeddings.py" \
    --sequences  "${SEQUENCES}" \
    --outdir     "${OUTDIR}" \
    --variant    "${VARIANT}" \
    --device     cuda \
    --batch-size 1

EXIT_CODE=$?

echo "============================================================"
echo "  Finished     : $(date)"
echo "  Exit code    : ${EXIT_CODE}"
echo "============================================================"

exit ${EXIT_CODE}
