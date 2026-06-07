#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=bmvp_conv
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_conv_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/slurm_logs/bmvp_conv_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-32%4
# =============================================================================
# 04i_finetune_BrainMVP_conv_submit_csd3.sh
#   BrainMVP (UniFormer) fine-tuning to COMPLETE the conversion tables. Matrix
#   (33 runs):
#     full_ft × {none, stochastic, plus_original} × T3b only  (T3a/T3c/T4 done)   = 9
#     frozen  × {stochastic, plus_original} × {T3a, T3b, T3c, T4} × seeds {0,1,2} = 24
#   Each run saves BOTH val_metrics and test_metrics (trainer final eval).
#
#   Output tree: brainmvp_debug/aug_<aug>/BrainMVP_uniformer/<task>/seed_<s>/<strategy>/
#   (already globbed by 05.MODEL_TREES "brainmvp_debug/aug_*/..." — auto-discovered).
#
# Run (from the output folder):
#   cd /home/ec474/rds/hpc-work/ADNI_MRI
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp/04i_finetune_BrainMVP_conv_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps (torch/monai/sklearn/pandas)"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline/brain_mvp"
PRETRAINED_CKPT="${ROOT}/ViT_pretrained/BrainMVP_uniformer.pt"
BRAINMVP_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
mkdir -p "${ROOT}/ADNI_MRI/brainmvp_debug/slurm_logs"

# ── Build the (strategy, augment, task, seed) job list ───────────────────────
SEEDS=(0 1 2)
JOBS=()
# full_ft T3b only (3 augments)
for aug in none stochastic plus_original; do
  for s in "${SEEDS[@]}"; do JOBS+=("full_ft ${aug} T3b_conv5y ${s}"); done
done
# frozen, 2 augments, all four conversion tasks
for aug in stochastic plus_original; do
  for t in T3a_conv3y T3b_conv5y T3c_conv7y T4_conv_horizon; do
    for s in "${SEEDS[@]}"; do JOBS+=("frozen ${aug} ${t} ${s}"); done
  done
done
echo "[info] ${#JOBS[@]} total jobs (array 0-$(( ${#JOBS[@]} - 1 )))"
read -r STRATEGY AUGMENT TASK SEED <<< "${JOBS[$SLURM_ARRAY_TASK_ID]}"

OUT_DIR="${ROOT}/ADNI_MRI/brainmvp_debug/aug_${AUGMENT}"
RUN_DIR="${OUT_DIR}/BrainMVP_uniformer/${TASK}/seed_${SEED}/${STRATEGY}"
export WANDB_MODE=offline
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  BrainMVP conv finetune | array ${SLURM_ARRAY_TASK_ID}"
echo "  strategy=${STRATEGY} augment=${AUGMENT} task=${TASK} seed=${SEED}"
echo "  out: ${RUN_DIR}   $(date)"
echo "============================================================"

if [ -f "${RUN_DIR}/metrics.json" ]; then
  echo "  metrics.json exists -- skipping."; exit 0
fi

conda run -n mri python "${SCRIPT_DIR}/04_supervised_finetuning_BrainMVP.py" \
    --task "${TASK}" --seed "${SEED}" --strategy "${STRATEGY}" \
    --augment "${AUGMENT}" --long all \
    --pretrained_ckpt "${PRETRAINED_CKPT}" \
    --brainmvp_inputs_dir "${BRAINMVP_INPUTS}" \
    --matched_labels_csv "${MATCHED}" --data_dir "${DATA_DIR}" \
    --out_dir "${OUT_DIR}" \
    --num_workers "${SLURM_CPUS_PER_TASK:-4}"
rc=$?
echo "  Finished ${SLURM_ARRAY_TASK_ID}: $(date)  exit=${rc}"
exit ${rc}
