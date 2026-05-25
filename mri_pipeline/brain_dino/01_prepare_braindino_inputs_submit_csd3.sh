#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=bdino_prep
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs/slurm_logs/bdino_prep_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs/slurm_logs/bdino_prep_%j.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
# =============================================================================
# 01_prepare_braindino_inputs_submit_csd3.sh  --  BrainDINO preprocessing
# =============================================================================
# Single CPU job: runs 01_prepare_braindino_inputs.py over the whole sMRIprep
# tree, producing the paper-fidelity finetune-recipe inputs (256, 256, 128) in
# RAS axis order = 128 axial x 256 cor x 256 sag, per arXiv:2604.27277
# Section 4.3.
#
# No GPU needed -- the preprocessing is CPU-only (MONAI transforms +
# nibabel I/O). Parallelised internally with --n-workers.
#
# Pipeline (per volume):
#   1. Load sMRIprep MNI preproc T1w + matching brain mask
#   2. Mask consistency check (info-only; logs mask_diff_frac to manifest;
#      aborts only on >20% catastrophic mismatch e.g. shape/wrong-subject)
#   3. Orient -> RAS, MaskIntensityd (strict skull-strip), CropForeground
#      on the mask, Resize to (256, 256, 128) trilinear
#   4. Percentile clip [0.5, 99.5] on non-zero voxels, then per-volume
#      z-score on the clipped non-zero voxels, restore background zeros
#
# Pre-flight (run once before submitting):
#   1. mri conda env available.
#   2. sMRIprep tree at ${SMRIPREP_DIR} (sessionwise layout).
#   3. mkdir -p /home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs/slurm_logs
#      (the #SBATCH --output dir must exist BEFORE submitting).
#   4. Account: COMPUTERLAB-SL3-CPU. CPU-only icelake job (LIO-CHARM-SL2
#      is GPU-only and cannot run on icelake).
#
# Submit:  sbatch mri_pipeline/brain_dino/01_prepare_braindino_inputs_submit_csd3.sh
# Resume:  re-running auto-skips (sub, ses) pairs already on disk
#          (per-output os.path.isfile check in process_subject).
#
# Post-run verification:
#   conda activate mri
#   python mri_pipeline/brain_dino/01b_check_braindino_inputs.py --sample 20
# =============================================================================

# -- Environment --------------------------------------------------------------
module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- Hardcoded paths (BASH_SOURCE unreliable: SLURM copies .sh to spool dir) --
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_dino"
SMRIPREP_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/smriprep_sessionwise/smriprep"
OUT_ROOT="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/braindino_inputs"

mkdir -p "${OUT_ROOT}/slurm_logs"

echo "============================================================"
echo "  BrainDINO preprocessing -- 01_prepare_braindino_inputs.py"
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

python "${SCRIPT_DIR}/01_prepare_braindino_inputs.py" \
    --long all \
    --n-workers    "${SLURM_CPUS_PER_TASK}" \
    --smriprep-dir "${SMRIPREP_DIR}" \
    --out-root     "${OUT_ROOT}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${OUT_ROOT}/"
echo "============================================================"
exit ${EXIT_CODE}
