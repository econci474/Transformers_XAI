#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_scratch_T1d
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/slurm_logs/vit_scratch_T1d_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/slurm_logs/vit_scratch_T1d_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%4
# =============================================================================
# 04b_finetune_ViT_scratch_T1d_backfill_csd3.sh
# =============================================================================
# Three-cell backfill: ViT-B from-scratch for T1d_binary (pMCI vs sMCI) at
# seeds 0/1/2. The main 04b sweep
# (mri_pipeline/04b_finetune_ViT_scratch_submit_csd3.sh) has TASKS without
# T1d, so the cross-model table currently shows ViT-scratch/baseline empty
# at T1d. This patches the gap.
#
# Strategy: identical to the main 04b sweep -- same trainer, same path
# tree, same wandb project (vit_baseline). Output lands at:
#   /home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline/ViT_B_scratch/T1d_binary/seed_<n>/scratch/
# matching the aggregator glob in 05_aggregate_mri_results.py.
#
# Expected outcome: same chance-floor as the other ViT-B-scratch tasks
# (T1, T1b, T1c, T2 all hit 0.500 val_bacc with 3 degenerate predictions).
# Running anyway so the cross-model table has T1d filled symmetrically.
#
# val_bacc, test_bacc, per-epoch loss/bacc curves, confusion matrices, and
# subject-level metrics are all logged automatically by the trainer:
#   - val_bacc + val_loss per epoch  -> W&B + metrics.json
#   - test_bacc + test_auc + test CM -> metrics.json + W&B summary
#   - test_subject_* metrics         -> metrics.json
#
# Submit:  sbatch mri_pipeline/baseline_vit/04b_finetune_ViT_scratch_T1d_backfill_csd3.sh
# Resume:  resubmitting auto-skips any cell whose metrics.json already exists.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

export WANDB_MODE=offline

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_baseline"

PY_SESSION_FLAGS="--long all"

# Three-cell sweep: TASK is fixed, only SEED varies (array IDs 0..2).
TASK="T1d_binary"
SEED="${SLURM_ARRAY_TASK_ID}"
STRATEGY="scratch"

RUN_DIR="${OUT_DIR}/ViT_B_scratch/${TASK}/seed_${SEED}/${STRATEGY}"
METRICS="${RUN_DIR}/metrics.json"

mkdir -p "${OUT_DIR}/slurm_logs"
mkdir -p "${RUN_DIR}"
export WANDB_DIR="${RUN_DIR}"

echo "============================================================"
echo "  ViT-from-scratch T1d backfill"
echo "  Job ID         : $SLURM_JOB_ID"
echo "  Array ID       : $SLURM_ARRAY_TASK_ID  (= seed)"
echo "  Node           : $SLURMD_NODENAME"
echo "  Task           : ${TASK}"
echo "  Seed           : ${SEED}"
echo "  Strategy       : ${STRATEGY} (random init, no MAE checkpoint)"
echo "  Output dir     : ${RUN_DIR}"
echo "  wandb          : offline -> project vit_baseline"
echo "  Started        : $(date)"
echo "============================================================"

if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs directory not found: ${VIT_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"; exit 1
fi
if [ -f "${METRICS}" ]; then
    echo "  metrics.json already exists -- skipping."; exit 0
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    ${PY_SESSION_FLAGS} \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --grad_accum_steps    8 \
    --wandb \
    --wandb_project       vit_baseline \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished       : $(date)"
echo "  Exit code      : ${EXIT_CODE}"
echo "  Output         : ${RUN_DIR}/"
echo "============================================================"
exit ${EXIT_CODE}
