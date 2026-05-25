#!/bin/bash
#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vitmae_extract
#SBATCH --output=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/vit_mae/slurm_logs/vitmae_extract_%j.log
#SBATCH --error=/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/vit_mae/slurm_logs/vitmae_extract_%j.err
#SBATCH -p ampere
#SBATCH --exclude=gpu-q-[15,64,68-75,81-82,86-90]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:30:00
# =============================================================================
# 03c_extract_vit_mae_embeddings_submit_csd3.sh
# =============================================================================
# Single non-array GPU job. Forwards every preprocessed ViT input scan
# through the frozen MAE-pretrained Vision_Transformer3D encoder, takes
# the post-norm CLS token (head replaced with nn.Identity), saves paired
# (parquet + .pt) keyed by (bids_sub, bids_ses).
#
# Expected wall-clock: ~15 min for ~3000 scans at batch=4 on A100.
# =============================================================================

module purge
module load cuda/12.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/cached_embeddings/vit_mae"

mkdir -p "${OUT_DIR}/slurm_logs"

echo "============================================================"
echo "  ViT-MAE75 frozen-encoder embedding extraction"
echo "  Job ID      : $SLURM_JOB_ID"
echo "  Node        : $SLURMD_NODENAME"
echo "  GPU         : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  Out dir     : ${OUT_DIR}"
echo "  Started     : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"; exit 1
fi
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs dir not found: ${VIT_INPUTS_DIR}"; exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"; exit 1
fi

python "${SCRIPT_DIR}/03c_extract_vit_mae_embeddings.py" \
    --pretrained_ckpt        "${PRETRAINED_CKPT}" \
    --vit_inputs_dir         "${VIT_INPUTS_DIR}" \
    --matched_labels_csv     "${MATCHED_LABELS_CSV}" \
    --out_dir                "${OUT_DIR}" \
    --batch_size             4 \
    --num_workers            "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?
echo "============================================================"
echo "  Finished    : $(date)"
echo "  Exit code   : ${EXIT_CODE}"
echo "  Output      : ${OUT_DIR}/{vit_mae_pooled.parquet,vit_mae_pooled.pt,manifest.csv}"
echo "============================================================"
exit ${EXIT_CODE}
