#!/bin/bash
# =============================================================================
# 03_prepare_ViT_submit_csd3.sh
# =============================================================================
# CPU-only SLURM job: applies MONAI preprocessing chain (orientation -> spacing
# -> z-score -> crop foreground -> resize 128^3) to sMRIprep MNI T1w volumes,
# producing vit_inputs/sub-*/sub-*_space-ViT128_desc-preproc_T1w.nii.gz.
#
# No GPU needed. Parallelised across CPUs with --n-workers.
#
# Submit:
#   mkdir -p logs
#   sbatch mri_pipeline/03_prepare_ViT_submit_csd3.sh
# =============================================================================

#SBATCH -A COMPUTERLAB-SL3-CPU
#SBATCH --job-name=vit_prep
#SBATCH --output=logs/vit_prep_%j.log
#SBATCH --error=logs/vit_prep_%j.err
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# ── Environment ───────────────────────────────────────────────────────────────
module purge

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Hardcoded paths (BASH_SOURCE doesn't work — SLURM copies .sh to spool dir)
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
SMRIPREP_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/smriprep_sessionwise/smriprep"
OUT_ROOT="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"

# Session selection — uncomment exactly ONE.
# Existing outputs from prior runs are auto-skipped (resume logic), so this
# can be re-run after extending the cutoff.
#PY_SESSION_FLAGS="--session bl"   # ses-bl only (~28 scans for SNP+MRI cohort)
PY_SESSION_FLAGS="--long all"      # every available session (~1,595 scans)
#PY_SESSION_FLAGS="--long 3"       # bl through m36 (~830 scans)
#PY_SESSION_FLAGS="--long 5"       # bl through m60 (~1,200 scans)
#PY_SESSION_FLAGS="--long all"     # every available session (~1,595 scans)

mkdir -p "${OUT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "  ViT preprocessing — CSD3 icelake"
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Node      : $SLURMD_NODENAME"
echo "  CPUs      : $SLURM_CPUS_PER_TASK"
echo "  Smriprep  : ${SMRIPREP_DIR}"
echo "  Output    : ${OUT_ROOT}"
echo "  Sessions  : ${PY_SESSION_FLAGS}"
echo "  Started   : $(date)"
echo "============================================================"

python "${SCRIPT_DIR}/03_prepare_ViT.py" \
    --smriprep-dir "${SMRIPREP_DIR}" \
    --out-root     "${OUT_ROOT}" \
    ${PY_SESSION_FLAGS} \
    --n-workers    "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "  Output    : ${OUT_ROOT} (manifest: vit_manifest.csv)"
echo "============================================================"

exit ${EXIT_CODE}
