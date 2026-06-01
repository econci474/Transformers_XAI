#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_embed
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/slurm_logs/bmvp_embed_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/slurm_logs/bmvp_embed_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-17%4
# =============================================================================
# 05_extract_brainmvp_full_ft_submit_csd3.sh
# =============================================================================
# Extracts per-scan embeddings + per-class probabilities from BrainMVP full_ft
# checkpoints for the (task, augment, seed) combos the user specified for
# downstream multimodal modelling. Per-task scope:
#
#   T1_binary     : aug = stochastic, plus_original     (BrainMVP best variants)
#   T1b_binary    : aug = stochastic, none
#   T1d_binary    : aug = plus_original
#   T2_multiclass : aug = stochastic
#
# 6 (task, aug) combinations x 3 seeds = 18 array tasks. Each task is a single
# GPU forward pass over the post-exclusion MRI master cohort (~1456 scans),
# wall-clock ~5-10 min on A100. %4 throttle per project convention.
#
# Inference output (per array task):
#   /home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/<task>/aug_<aug>/seed_<n>/
#       embeddings_seed_<n>.npz   (~3 MB, includes IDs + embeddings + probs)
#       embeddings_seed_<n>.csv   (flat CSV without embedding columns)
#
# After the array completes, scp the output tree back to local D::
#   scp -r ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings \
#       /d/ADNI_BIDS_project/derivatives/mri_embeddings_brainmvp/
#
# Submit:  sbatch mri_pipeline/brain_mvp/05_extract_brainmvp_full_ft_submit_csd3.sh
# =============================================================================

set -euo pipefail

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Hardcoded paths (BASH_SOURCE unreliable: SLURM copies the .sh to a spool dir)
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
EXTRACTOR="${SCRIPT_DIR}/05_extract_brainmvp_full_ft_embeddings.py"

INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
CKPT_ROOT="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug/brainmvp_debug"
MRI_MASTER="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
SPLITS_ROOT="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings"

mkdir -p "${OUT_DIR}/slurm_logs"

# Decode SLURM_ARRAY_TASK_ID -> (task, aug, seed)
# CELLS list: 6 (task, aug) combos × 3 seeds = 18 entries
# Order: outer loop over (task, aug); inner loop over seed.
# Index TASK_ID 0..17 maps to entry [TASK_ID // 3] of CELLS, seed = TASK_ID % 3.

CELLS=(
    "T1_binary stochastic"
    "T1_binary plus_original"
    "T1b_binary stochastic"
    "T1b_binary none"
    "T1d_binary plus_original"
    "T2_multiclass stochastic"
)

CELL_IDX=$(( SLURM_ARRAY_TASK_ID / 3 ))
SEED=$(( SLURM_ARRAY_TASK_ID % 3 ))

if [ "${CELL_IDX}" -ge "${#CELLS[@]}" ]; then
    echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (max=$((${#CELLS[@]}*3 - 1)))"
    exit 1
fi

CELL="${CELLS[$CELL_IDX]}"
TASK="$(echo "${CELL}" | awk '{print $1}')"
AUG="$(echo "${CELL}" | awk '{print $2}')"

# Skip if output already exists (idempotent resume).
OUT_NPZ="${OUT_DIR}/${TASK}/aug_${AUG}/seed_${SEED}/embeddings_seed_${SEED}.npz"
if [ -f "${OUT_NPZ}" ]; then
    echo "[skip] Output already exists: ${OUT_NPZ}"
    exit 0
fi

echo "============================================================"
echo "  BrainMVP full_ft embedding extraction"
echo "  Job ID    : ${SLURM_JOB_ID}"
echo "  Array ID  : ${SLURM_ARRAY_TASK_ID}  (cell ${CELL_IDX}/${#CELLS[@]}, seed ${SEED})"
echo "  Node      : ${SLURMD_NODENAME}"
echo "  GPU       : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Task      : ${TASK}"
echo "  Augment   : ${AUG}"
echo "  Seed      : ${SEED}"
echo "  Output    : ${OUT_NPZ}"
echo "  Started   : $(date)"
echo "============================================================"

# Run extraction
python "${EXTRACTOR}" \
    --task "${TASK}" \
    --augment "${AUG}" \
    --seed "${SEED}" \
    --inputs-dir "${INPUTS_DIR}" \
    --ckpt-root "${CKPT_ROOT}" \
    --master-csv "${MRI_MASTER}" \
    --splits-root "${SPLITS_ROOT}" \
    --out-dir "${OUT_DIR}" \
    --batch-size 4 \
    --num-workers 4

EXIT_CODE=$?

echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
