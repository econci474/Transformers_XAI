#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=cnn_prep
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs/slurm_logs/cnn_prep_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs/slurm_logs/cnn_prep_%j.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# =============================================================================
# 00_prepare_CNN_submit_csd3.sh   --   Spasov 3D CNN preprocessing
# =============================================================================
# Single CPU job: runs 00_prepare_CNN_inputs.py over the whole sMRIprep tree,
# producing z-scored 3D-CNN inputs (193x229x193, 1 mm, RAS — native MNI grid,
# no resample/crop). Parallelised internally with --n-workers.
#
# No GPU needed — the preprocessing is CPU-only (MONAI transforms).
#
# Pre-flight (run once before submitting):
#   1. mri conda env available.
#   2. sMRIprep tree at ${SMRIPREP_DIR}.
#   3. mkdir -p /home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
#   4. Account: COMPUTERLAB-SL3-CPU. This is a CPU job, so it cannot run on
#      LIO-CHARM (GPU-only); COMPUTERLAB-SL3-CPU is the CPU allocation with
#      hours available.
#
# Submit:  sbatch mri_pipeline/3d_conv_net/00_prepare_CNN_submit_csd3.sh
# Resume:  re-running skips (sub, ses) pairs already in the manifest.
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/3d_conv_net"
SMRIPREP_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/smriprep/smriprep"
OUT_ROOT="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/cnn_inputs"

mkdir -p "${OUT_ROOT}/slurm_logs"

echo "============================================================"
echo "  3D-CNN preprocessing -- 00_prepare_CNN_inputs.py"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Node        : $SLURMD_NODENAME"
echo "  CPUs        : $SLURM_CPUS_PER_TASK"
echo "  sMRIprep    : ${SMRIPREP_DIR}"
echo "  Output      : ${OUT_ROOT}"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -d "${SMRIPREP_DIR}" ]; then
    echo "[ERROR] sMRIprep directory not found: ${SMRIPREP_DIR}"
    exit 1
fi

python "${SCRIPT_DIR}/00_prepare_CNN_inputs.py" \
    --long all \
    --n-workers   "${SLURM_CPUS_PER_TASK}" \
    --smriprep-dir "${SMRIPREP_DIR}" \
    --out-root    "${OUT_ROOT}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${OUT_ROOT}/"
echo "============================================================"
exit ${EXIT_CODE}
