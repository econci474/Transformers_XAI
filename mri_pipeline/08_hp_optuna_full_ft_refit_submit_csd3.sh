#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=optuna_refit
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs/optuna_refit_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs/optuna_refit_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-8%4
# =============================================================================
# 08_hp_optuna_full_ft_refit_submit_csd3.sh
#   Refit the top-k Optuna configs on 3 seeds (KEEP checkpoints), then pick the
#   winner by mean val-bACC.
# =============================================================================
# Array size = topk * n_seeds (default 3 * 3 = 9 -> array 0-8). The driver
# decodes the array index into (config_k, seed) from <study-dir>/top_configs.json
# (written by `--mode report`). Each cell is a full full_ft run that saves
# best_model.pt under <winner-dir>/config_<k>/.
#
# Run AFTER the search + report:
#   sbatch --export=ALL,ARCH=brainmvp mri_pipeline/08_hp_optuna_full_ft_refit_submit_csd3.sh
#   sbatch --export=ALL,ARCH=vit_mae  mri_pipeline/08_hp_optuna_full_ft_refit_submit_csd3.sh
# Then pick the winner (light, any node):
#   python mri_pipeline/08_hp_optuna_full_ft.py --mode winner --arch $ARCH \
#       --winner-dir /home/ec474/rds/hpc-work/ADNI_MRI/optuna/${ARCH}_T2_winner
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
# Env-independent: do NOT `conda activate` (sbatch inherits the login-node env and
# activation can stack on top of it). Target the named env explicitly with
# `conda run -n mri`, and preflight the deps so a wrong/missing env fails fast.
conda run -n mri python -c "import torch, optuna" \
    || { echo "[ERROR] mri env missing torch/optuna"; exit 1; }

ARCH="${ARCH:?set ARCH=brainmvp or ARCH=vit_mae via --export=ALL,ARCH=...}"
export WANDB_MODE="${WANDB_MODE:-offline}"
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
STUDY_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/optuna/${ARCH}_T2"
WINNER_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/optuna/${ARCH}_T2_winner"
mkdir -p "${WINNER_DIR}" /home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs
export WANDB_DIR="${WINNER_DIR}"

echo "============================================================"
echo "  Optuna full_ft REFIT — arch=${ARCH}  refit-index ${SLURM_ARRAY_TASK_ID}"
echo "  Study  : ${STUDY_DIR}   Winner: ${WINNER_DIR}"
echo "  Started: $(date)"
echo "============================================================"

conda run -n mri python "${SCRIPT_DIR}/08_hp_optuna_full_ft.py" \
    --mode refit --arch "${ARCH}" --task T2_multiclass \
    --study-dir "${STUDY_DIR}" --winner-dir "${WINNER_DIR}" \
    --refit-index "${SLURM_ARRAY_TASK_ID}"
rc=$?

echo "  Finished: $(date)  exit=${rc}"
exit ${rc}
