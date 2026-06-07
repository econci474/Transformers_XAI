#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=agms3d_van
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs/slurm_logs/agms3d_van_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/agms3d_outputs/slurm_logs/agms3d_van_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-14%4
# =============================================================================
# 04k_finetune_AGMS3D_vanilla_submit_csd3.sh
#   RE-TRAIN AG-MS3D (vanilla backbone) to fill its cross-model val AUC/F1 — the
#   original run saved no checkpoint, so it can't be recomputed. 15 runs:
#     5 diagnostic tasks {T1,T1b,T1c,T1d,T2} × seeds {0,1,2}.
#   Recipe mirrors the original vanilla run (backbone=vanilla, head=large,
#   base_filters=32, lr=3e-3, ls=0.0, epochs=60, bs=8). The trainer now saves a
#   `val_metrics` block at final eval (Part A), so each run stores BOTH val + test.
#
#   NOTE: AG-MS3D collapsed to chance (bACC 0.500) in every original cell — expect
#   near-chance val AUC/F1 (this is an evidence-of-failure row).
#
#   Output: agms3d_outputs/AGMS3DCNN_vanilla/<task>/seed_<s>/  (05.MODEL_TREES glob).
#
# Run (from the output folder):
#   cd /home/ec474/rds/hpc-work/ADNI_MRI
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/3d_cnn_vit/04k_finetune_AGMS3D_vanilla_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline/3d_cnn_vit"
CNN_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/cnn_inputs"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="${ROOT}/ADNI_MRI/agms3d_outputs"
mkdir -p "${OUT_DIR}/slurm_logs"

TASKS=(T1_binary T1b_binary T1c_binary T1d_binary T2_multiclass)
SEEDS=(0 1 2)
ID=${SLURM_ARRAY_TASK_ID}
SEED_IDX=$(( ID % 3 )); TASK_IDX=$(( ID / 3 ))
TASK="${TASKS[$TASK_IDX]}"; SEED="${SEEDS[$SEED_IDX]}"

RUN_DIR="${OUT_DIR}/AGMS3DCNN_vanilla/${TASK}/seed_${SEED}"
export WANDB_MODE=offline
export WANDB_DIR="${RUN_DIR}"
echo "=== AG-MS3D vanilla | array ${ID} | task=${TASK} seed=${SEED} | $(date) ==="
if [ -f "${RUN_DIR}/metrics.json" ]; then echo "  exists -- skipping."; exit 0; fi

conda run -n mri python "${SCRIPT_DIR}/train_agms3d.py" \
    --task "${TASK}" --seed "${SEED}" \
    --backbone vanilla --head large --base_filters 32 \
    --lr 3e-3 --label_smoothing 0.0 --epochs 60 --batch_size 8 \
    --cnn_inputs_dir "${CNN_INPUTS}" \
    --matched_labels_csv "${MATCHED}" --data_dir "${DATA_DIR}" \
    --out_dir "${OUT_DIR}" \
    --num_workers "${SLURM_CPUS_PER_TASK:-4}"
rc=$?
echo "  Finished ${ID}: $(date)  exit=${rc}"
exit ${rc}
