#!/bin/bash
#SBATCH -A COMPUTERLAB-SL3-GPU
#SBATCH --job-name=vit_diag
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs/vit_diag_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/slurm_logs/vit_diag_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-44%4
# =============================================================================
# 04j_finetune_ViT_diag_redo_submit_csd3.sh
#   RE-RUN the ViT-MAE DIAGNOSTIC cells whose checkpoints were deleted (so their
#   val AUC/F1 can't be recomputed) — only the non-`random` augments are affected
#   (full_ft/random + frozen/random survive on D: and are already filled). Matrix
#   (45 runs):
#     full_ft × {none, plus_original} × {T1,T1b,T1c,T1d,T2} × seeds{0,1,2} = 30
#     frozen  × {plus_original}       × {T1,T1b,T1c,T1d,T2} × seeds{0,1,2} = 15
#   Each run saves BOTH val_metrics and test_metrics (trainer final eval).
#
#   Output tree: vit_outputs_hi_lr/aug_<aug>/ViT_B_mae75/<task>/seed_<s>/<strategy>/
#   (already globbed by 05.MODEL_TREES). DIAGNOSTIC tasks -> baseline split (no T4).
#
#   PREREQ: delete the stale (no-val, no-ckpt) cells first so skip-if-exists re-runs:
#     rm -rf $HOME/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/aug_none/ViT_B_mae75
#     rm -rf $HOME/rds/hpc-work/ADNI_MRI/vit_outputs_hi_lr/aug_plus_original/ViT_B_mae75
#   (Do NOT touch the unsuffixed vit_outputs_hi_lr/ViT_B_mae75 = full_ft/random,
#    nor vit_outputs_debug = frozen/random — those are filled.)
#
# Run (from the output folder):
#   cd /home/ec474/rds/hpc-work/ADNI_MRI
#   sbatch /home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/04j_finetune_ViT_diag_redo_submit_csd3.sh
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda run -n mri python -c "import torch, monai, sklearn, pandas" \
    || { echo "[ERROR] mri env missing deps (torch/monai/sklearn/pandas)"; exit 1; }

ROOT="/home/ec474/rds/hpc-work"
SCRIPT_DIR="${ROOT}/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="${ROOT}/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS="${ROOT}/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED="${ROOT}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DATA_DIR="${ROOT}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
mkdir -p "${ROOT}/ADNI_MRI/vit_outputs_hi_lr/slurm_logs"

# ── (strategy, augment, task, seed) job list ────────────────────────────────
TASKS=(T1_binary T1b_binary T1c_binary T1d_binary T2_multiclass)
SEEDS=(0 1 2)
CONFIGS=("full_ft none" "full_ft plus_original" "frozen plus_original")
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

OUT_DIR="${ROOT}/ADNI_MRI/vit_outputs_hi_lr/aug_${AUGMENT}"
RUN_DIR="${OUT_DIR}/ViT_B_mae75/${TASK}/seed_${SEED}/${STRATEGY}"
export WANDB_MODE=offline
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-MAE diag redo | array ${SLURM_ARRAY_TASK_ID}"
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
