#!/bin/bash
#SBATCH -A COMPUTERLAB-SL2-GPU
#SBATCH --job-name=optuna_ft
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs/optuna_ft_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs/optuna_ft_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-23%4
# =============================================================================
# 08_hp_optuna_full_ft_search_submit_csd3.sh
#   Optuna Bayesian full_ft HP SEARCH for T2 (BrainMVP or ViT-MAE), seed 0.
# =============================================================================
# Each array task = ONE Optuna trial (one full_ft run at seed 0) against a
# SHARED file-based JournalStorage study, so the %4 concurrent agents jointly
# explore the space (TPE; first ~10 trials random, then Bayesian). 24 trials.
# Each trial drops its *.pt after reading metrics.json (checkpoint-free search).
#
# Pick the model with AUGMENT-less ARCH export:
#   sbatch --export=ALL,ARCH=brainmvp mri_pipeline/08_hp_optuna_full_ft_search_submit_csd3.sh
#   sbatch --export=ALL,ARCH=vit_mae  mri_pipeline/08_hp_optuna_full_ft_search_submit_csd3.sh
# Smoke one trial:  sbatch --array=0 --export=ALL,ARCH=brainmvp ...
# Cap search epochs (faster trials):  --export=ALL,ARCH=brainmvp,SEARCH_EPOCHS=40
#
# Pre-flight (once):  pip install optuna   (into the `mri` env)
#                     mkdir -p /home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs
# After the array finishes -> report (writes top_configs.json):
#   python mri_pipeline/08_hp_optuna_full_ft.py --mode report --arch $ARCH \
#       --study-dir /home/ec474/rds/hpc-work/ADNI_MRI/optuna/${ARCH}_T2 --topk 3
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
mkdir -p "${STUDY_DIR}" /home/ec474/rds/hpc-work/ADNI_MRI/optuna/slurm_logs
export WANDB_DIR="${STUDY_DIR}"

SEARCH_EPOCHS_FLAG=""
[ -n "${SEARCH_EPOCHS}" ] && SEARCH_EPOCHS_FLAG="--search-epochs ${SEARCH_EPOCHS}"

echo "============================================================"
echo "  Optuna full_ft SEARCH — arch=${ARCH}  trial agent ${SLURM_ARRAY_TASK_ID}"
echo "  Study  : ${STUDY_DIR}"
echo "  GPU    : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Started: $(date)"
echo "============================================================"

conda run -n mri python "${SCRIPT_DIR}/08_hp_optuna_full_ft.py" \
    --mode search --arch "${ARCH}" --task T2_multiclass \
    --study-dir "${STUDY_DIR}" --n-trials 1 ${SEARCH_EPOCHS_FLAG}
rc=$?

echo "  Finished: $(date)  exit=${rc}"
exit ${rc}
