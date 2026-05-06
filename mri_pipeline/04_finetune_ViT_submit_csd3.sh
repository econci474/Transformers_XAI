#!/bin/bash
# =============================================================================
# 04_finetune_ViT_submit_csd3.sh
# =============================================================================
# Single-config smoke test for the supervised fine-tuning of the MAE-pretrained
# ViT-B/3D on ADNI baseline T1w MRIs (ses-bl).
#
#   Task     : T1_binary (CN vs MCI+AD)
#   Seed     : 0
#   Strategy : full_ft
#
# Once green, expand to a 4 x 3 x 2 = 24 array (mirrors the clinical pipeline).
#
# Pre-flight on a CSD3 login node:
#   1. mri conda env created with: torch, monai, nibabel, pandas, scikit-learn,
#      einops, tqdm  (pip install einops scikit-learn pyyaml tqdm if needed).
#   2. scp inputs to RDS:
#        scp ViT_B_pretrained_..._077000.pth.tar ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/pretrained/
#        scp -r D:/.../vit_inputs/                ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/
#        scp -r D:/.../mri_clinical_matched/      ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/
#        scp -r D:/.../verbose/baseline/          ec474@login-cpu.hpc.cam.ac.uk:/home/ec474/rds/hpc-work/ADNI_CL/
#      (or run 03b on CSD3 directly to regenerate the matched CSV in-place).
#   3. mkdir -p logs
#   4. sbatch mri_pipeline/04_finetune_ViT_submit_csd3.sh
# =============================================================================

#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=vit_ft
#SBATCH --output=logs/vit_ft_%j.log
#SBATCH --error=logs/vit_ft_%j.err
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# ── Environment ───────────────────────────────────────────────────────────────
module purge
module load cuda/12.1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# Hardcoded paths (BASH_SOURCE doesn't work — SLURM copies .sh to spool dir)
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/mri_pipeline"
PRETRAINED_CKPT="/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
DATA_DIR="/home/ec474/rds/hpc-work/ADNI_CL/baseline"
VIT_INPUTS_DIR="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
MATCHED_LABELS_CSV="/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/mri_clinical_matched/master_mri_clinical_matched_labels.csv"
OUT_DIR="/home/ec474/rds/hpc-work/ADNI_MRI/vit_outputs"

# ── Smoke-test config ─────────────────────────────────────────────────────────
TASK="T1_binary"
SEED=0
STRATEGY="full_ft"
SESSION="bl"

mkdir -p "${OUT_DIR}"
mkdir -p logs

echo "============================================================"
echo "  ViT supervised fine-tuning — smoke test"
echo "  Job ID     : $SLURM_JOB_ID"
echo "  Node       : $SLURMD_NODENAME"
echo "  GPU        : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  CUDA_DEV   : $CUDA_VISIBLE_DEVICES"
echo "  Task       : ${TASK}"
echo "  Seed       : ${SEED}"
echo "  Strategy   : ${STRATEGY}"
echo "  Session    : ${SESSION}"
echo "  Pretrained : ${PRETRAINED_CKPT}"
echo "  Started    : $(date)"
echo "============================================================"

if [ ! -f "${PRETRAINED_CKPT}" ]; then
    echo "[ERROR] Pretrained checkpoint not found: ${PRETRAINED_CKPT}"
    exit 1
fi
if [ ! -d "${VIT_INPUTS_DIR}" ]; then
    echo "[ERROR] ViT inputs directory not found: ${VIT_INPUTS_DIR}"
    echo "        Run 03_prepare_ViT_submit_csd3.sh first."
    exit 1
fi
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV not found: ${MATCHED_LABELS_CSV}"
    echo "        Run 03b_match_mri_to_clinical.py first."
    exit 1
fi
if [ ! -d "${DATA_DIR}/seed_${SEED}" ]; then
    echo "[ERROR] Clinical splits dir not found: ${DATA_DIR}/seed_${SEED}"
    exit 1
fi

python "${SCRIPT_DIR}/04_supervised_finetuning_ViT.py" \
    --task                "${TASK}" \
    --seed                "${SEED}" \
    --strategy            "${STRATEGY}" \
    --session             "${SESSION}" \
    --pretrained_ckpt     "${PRETRAINED_CKPT}" \
    --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
    --data_dir            "${DATA_DIR}" \
    --vit_inputs_dir      "${VIT_INPUTS_DIR}" \
    --out_dir             "${OUT_DIR}" \
    --num_workers         "${SLURM_CPUS_PER_TASK}"

EXIT_CODE=$?

echo "============================================================"
echo "  Finished  : $(date)"
echo "  Exit code : ${EXIT_CODE}"
echo "  Output    : ${OUT_DIR}/ViT_B_mae75/${TASK}/seed_${SEED}/${STRATEGY}/"
echo "============================================================"

exit ${EXIT_CODE}
