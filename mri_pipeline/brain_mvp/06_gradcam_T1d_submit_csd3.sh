#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=gradcam_T1d
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d/slurm_logs/gradcam_T1d_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d/slurm_logs/gradcam_T1d_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
# =============================================================================
# 06_gradcam_T1d_submit_csd3.sh
# =============================================================================
# Stage 1 (GPU compute) of the T1d BrainMVP Grad-CAM. Runs 06_gradcam_T1d.py
# over the m12 T1d test cohort for seeds 0/1/2 and writes per-scan 96^3 saliency
# arrays + manifest.csv to gradcam_T1d/. The 128^3 inputs + full_ft checkpoints
# live only on CSD3, so this stage must run here; the figure/atlas rendering is
# done locally afterwards by 07_render_gradcam_T1d.py.
#
# Single short GPU job (~75 forward+backward passes). SL3-GPU (SL2-GPU over
# quota). No array -> no %4 throttle needed.
#
# Submit (from the output folder so any relative logs land beside the outputs):
#   cd /home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp/06_gradcam_T1d_submit_csd3.sh
#
# Then pull results to local D: for rendering (07_render_gradcam_T1d.py):
#   rsync -av ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d \
#       /d/ADNI_BIDS_project/derivatives/
# =============================================================================

# NOTE: no `set -euo pipefail` (CSD3 conda.sh touches unbound PS1 under -u).

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"

# Hardcoded (SLURM copies the .sh to a spool dir, so BASH_SOURCE is unreliable).
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
SCRIPT="${SCRIPT_DIR}/06_gradcam_T1d.py"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/gradcam_T1d"

mkdir -p "${OUT_DIR}/slurm_logs"

# Env-independent: use `conda run -n mri` (do NOT `conda activate`, which stacks
# on whatever env was active at submit time). Preflight the imports and fail loud.
conda run -n mri python -c "import torch, monai, nibabel, numpy, pandas; print('preflight OK; cuda', torch.cuda.is_available())" || exit 1

echo "============================================================"
echo "  T1d BrainMVP Grad-CAM (stage 1: compute)"
echo "  Job ID  : ${SLURM_JOB_ID}"
echo "  Node    : ${SLURMD_NODENAME}"
echo "  GPU     : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Out     : ${OUT_DIR}"
echo "  Started : $(date)"
echo "============================================================"

conda run -n mri python "${SCRIPT}" \
    --seeds 0 1 2 \
    --viscode m12 \
    --target true \
    --out-dir "${OUT_DIR}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
