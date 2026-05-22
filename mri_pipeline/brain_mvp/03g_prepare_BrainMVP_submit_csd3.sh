#!/bin/bash
# =============================================================================
# 03g_prepare_BrainMVP_submit_csd3.sh
# =============================================================================
# CPU-only SLURM job: applies BrainMVP preprocessing chain (orientation ->
# spacing 1mm -> percentile clip [5th,95th] -> [0,1] -> crop foreground ->
# resize 128x128x64) to sMRIprep MNI T1w volumes, producing:
#   brainmvp_inputs/sub-*/ses-*/sub-*_space-BrainMVP128x64_desc-preproc_T1w.nii.gz
#
# No GPU needed. Parallelised across CPUs with --n-workers.
#
# Submit:
#   mkdir -p logs
#   sbatch mri_pipeline/brain_mvp/03g_prepare_BrainMVP_submit_csd3.sh
# =============================================================================

#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=bmvp_prep
#SBATCH --output=logs/bmvp_prep_%j.log
#SBATCH --error=logs/bmvp_prep_%j.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00

# ── Environment ───────────────────────────────────────────────────────────────
module purge

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Hardcoded paths (BASH_SOURCE doesn't work — SLURM copies .sh to spool dir)
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
SMRIPREP_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/smriprep_sessionwise/smriprep"
OUT_ROOT="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"

# Session selection — uncomment exactly ONE.
# Existing outputs from prior runs are auto-skipped (resume logic), so this
# can be re-run after extending the cutoff.
#PY_SESSION_FLAGS="--session bl"   # ses-bl only
PY_SESSION_FLAGS="--long all"      # every available session

mkdir -p "${OUT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "  BrainMVP preprocessing (128x128x64 @ 1mm) — CSD3 icelake"
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Node      : $SLURMD_NODENAME"
echo "  CPUs      : $SLURM_CPUS_PER_TASK"
echo "  Smriprep  : ${SMRIPREP_DIR}"
echo "  Output    : ${OUT_ROOT}"
echo "  Sessions  : ${PY_SESSION_FLAGS}"
echo "  Started   : $(date)"
echo "============================================================"

python "${SCRIPT_DIR}/03_prepare_BrainMVP.py" \
    --smriprep-dir "${SMRIPREP_DIR}" \
    --out-root     "${OUT_ROOT}" \
    ${PY_SESSION_FLAGS} \
    --n-workers    "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "  Output    : ${OUT_ROOT} (manifest: brainmvp_manifest.csv)"
echo "============================================================"

exit ${EXIT_CODE}
