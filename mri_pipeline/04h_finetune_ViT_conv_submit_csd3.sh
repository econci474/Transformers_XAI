#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=vit_conv
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_conv/slurm_logs/vit_conv_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_conv/slurm_logs/vit_conv_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-59%4
# =============================================================================
# 04h_finetune_ViT_conv_submit_csd3.sh
#   ViT-MAE (ViT-B/3D, MAE-pretrained) fine-tuning for the CONVERSION tasks
#   T3a/b/c + T4, to complete the cross-model conversion tables. Matrix (60 runs):
#     full_ft × {none, random, plus_original}   (ViT 'random' renders as stochastic¹)
#     frozen  × {random, plus_original}
#     × tasks {T3a_conv3y, T3b_conv5y, T3c_conv7y, T4_conv_horizon} × seeds {0,1,2}
#   Each run saves BOTH val_metrics and test_metrics (trainer final eval).
#
#   Output tree: vit_outputs_conv/aug_<aug>/ViT_B_mae75/<task>/seed_<s>/<strategy>/
#   (picked up by the 05.MODEL_TREES "vit_outputs_conv/aug_*/..." glob).
#
# Run (from the output folder so relative paths resolve there):
#   cd /home/ec474/rds/hpc-work/ADNI_MRI
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/04h_finetune_ViT_conv_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
# Env-independent: target the named env explicitly (do NOT `conda activate`).
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps (torch/monai/sklearn/pandas)"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="${ROOT}/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
mkdir -p "${ROOT}/ADNI_MRI/vit_outputs_conv/slurm_logs"

# ── Build the (strategy, augment, task, seed) job list ───────────────────────
TASKS=(T3a_conv3y T3b_conv5y T3c_conv7y T4_conv_horizon)
SEEDS=(0 1 2)
CONFIGS=("full_ft none" "full_ft random" "full_ft plus_original" \
         "frozen random" "frozen plus_original")
JOBS=()
for cfg in "${CONFIGS[@]}"; do
  for t in "${TASKS[@]}"; do
    for s in "${SEEDS[@]}"; do
      JOBS+=("${cfg} ${t} ${s}")
    done
  done
done
echo "[info] ${#JOBS[@]} total jobs (array 0-$(( ${#JOBS[@]} - 1 )))"
read -r STRATEGY AUGMENT TASK SEED <<< "${JOBS[$SLURM_ARRAY_TASK_ID]}"

OUT_DIR="${ROOT}/ADNI_MRI/vit_outputs_conv/aug_${AUGMENT}"
RUN_DIR="${OUT_DIR}/ViT_B_mae75/${TASK}/seed_${SEED}/${STRATEGY}"
export WANDB_MODE=offline
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-MAE conv finetune | array ${SLURM_ARRAY_TASK_ID}"
echo "  strategy=${STRATEGY} augment=${AUGMENT} task=${TASK} seed=${SEED}"
echo "  out: ${RUN_DIR}   $(date)"
echo "============================================================"

if [ -f "${RUN_DIR}/metrics.json" ]; then
  echo "  metrics.json exists -- skipping."; exit 0
fi

conda run -n mri python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task "${TASK}" --seed "${SEED}" --strategy "${STRATEGY}" \
    --vit_size base --augment "${AUGMENT}" --long all \
    --pretrained_ckpt "${PRETRAINED_CKPT}" \
    --vit_inputs_dir "${VIT_INPUTS}" \
    --matched_labels_csv "${MATCHED}" --data_dir "${DATA_DIR}" \
    --out_dir "${OUT_DIR}" \
    --num_workers "${SLURM_CPUS_PER_TASK:-4}"
rc=$?
echo "  Finished ${SLURM_ARRAY_TASK_ID}: $(date)  exit=${rc}"
exit ${rc}
