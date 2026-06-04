#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmvp_embed_T3
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/slurm_logs/bmvp_embed_T3_%A_%a.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/slurm_logs/bmvp_embed_T3_%A_%a.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-5%4
# =============================================================================
# 05_extract_brainmvp_full_ft_T3_submit_csd3.sh
# =============================================================================
# Per-patient probability extraction for BrainMVP full_ft (plus_original) on the
# T3 binary conversion horizons that were full-fine-tuned (T3a, T3c) — runs the
# fine-tuned encoder over the MRI volumes (the cached frozen embeddings can't be
# reused for full_ft). Feeds the upgraded T3/T4 late-fusion.
#
# 2 (task, aug) combos × 3 seeds = 6 array tasks. Output (per array task):
#   /home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/<task>/aug_plus_original/seed_<n>/
#       embeddings_seed_<n>.{npz,csv}   (csv = Patient_ID, adni_viscode, split, label, prob_class_*)
#
# Submit:  sbatch mri_pipeline/brain_mvp/05_extract_brainmvp_full_ft_T3_submit_csd3.sh
# Then scp the CSVs to local:
#   scp -r ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings/T3a_conv3y \
#       /d/ADNI_BIDS_project/derivatives/brainmvp_embeddings/
#   (and T3c_conv7y)
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline/brain_mvp"
EXTRACTOR="${SCRIPT_DIR}/05_extract_brainmvp_full_ft_embeddings.py"

INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"
CKPT_ROOT="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_debug"
MRI_MASTER="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
SPLITS_ROOT="/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/brainmvp_embeddings"

mkdir -p "${OUT_DIR}/slurm_logs"

# CELLS: 2 (task, aug) combos × 3 seeds = 6 array tasks. seed = TASK_ID % 3.
CELLS=(
    "T3a_conv3y plus_original"
    "T3c_conv7y plus_original"
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

OUT_NPZ="${OUT_DIR}/${TASK}/aug_${AUG}/seed_${SEED}/embeddings_seed_${SEED}.npz"
if [ -f "${OUT_NPZ}" ]; then
    echo "[skip] Output already exists: ${OUT_NPZ}"; exit 0
fi

echo "============================================================"
echo "  BrainMVP full_ft embedding extraction -- T3"
echo "  Array ID  : ${SLURM_ARRAY_TASK_ID}  (cell ${CELL_IDX}/${#CELLS[@]}, seed ${SEED})"
echo "  Task: ${TASK}   Augment: ${AUG}   Seed: ${SEED}"
echo "  Output: ${OUT_NPZ}   Started: $(date)"
echo "============================================================"

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
echo "  Finished  : $(date)   Exit code : ${EXIT_CODE}"
echo "============================================================"
exit ${EXIT_CODE}
