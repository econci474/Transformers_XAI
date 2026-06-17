#!/bin/bash
# =============================================================================
# 02_finetune.sh (gwas_classification)
# =============================================================================
# SLURM job script for BMFM-DNA-SNP fine-tuning on binned z-score
# classification (7 classes: {-3, -2, -1, 0, +1, +2, +3}).
#
# Model  : ibm-research/biomed.dna.snp.modernbert.113m.v1  (ModernBERT, 113M)
# Task   : Classification — predict binned z-score from allele-specific DNA
# Data   : ~350 GW-sig × 2 alleles + N null SNPs  (from 08h)
# Input  : bmfm_gwas_classification_without_ukb/ (train/dev/test.csv)
# Config : 02_finetune.yaml
#
# GPU estimate:
#   ModernBERT 113M @ bf16, 2048 tokens, batch=4 → ~6 GB VRAM.
#   With ~280k+ rows, training is longer than regression (hours).
#
# Usage:
#   sbatch 02_finetune.sh
#
# Transfer dataset to CSD3 first:
#   scp -r D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_classification_without_ukb \
#       ec474@login-cpu.hpc.cam.ac.uk:~/ADNI_SNP/
# =============================================================================


#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmfm_gwas_cls
#SBATCH --output=logs/bmfm_gwas_cls_%j.out
#SBATCH --error=logs/bmfm_gwas_cls_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1             # 1 A100 40 GB
#SBATCH --cpus-per-task=3        # CSD3 max: 3 CPUs per GPU
#SBATCH -p ampere                # Cambridge CSD3 GPU partition (A100s)

# =============================================================================
# ── User-configurable paths ───────────────────────────────────────────────────
# =============================================================================

CONDA_ENV="bmfm"

# Classification dataset (08h, generated locally then scp'd to CSD3):
INPUT_DIRECTORY="/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_classification_without_ukb"

# Where to write checkpoints, predictions, logs:
OUTPUT_DIRECTORY="/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_classification_output/cls_snp_lr1e5_5ep"

# Directory containing 02_finetune.yaml (this script's directory).
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline/gwas_classification"
PARENT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline"

# =============================================================================
# ── Environment ───────────────────────────────────────────────────────────────
# =============================================================================

echo "============================================================"
echo " BMFM-DNA-SNP GWAS Classification Fine-tuning"
echo " Job ID      : $SLURM_JOB_ID"
echo " Node        : $SLURMD_NODENAME"
echo " GPU(s)      : $CUDA_VISIBLE_DEVICES"
echo " Input dir   : $INPUT_DIRECTORY"
echo " Output dir  : $OUTPUT_DIRECTORY"
echo " Start time  : $(date)"
echo "============================================================"

# Activate conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Flash Attention 2 glibc compat shim
export LD_PRELOAD="/home/ec474/glibc_compat.so"

# Verify bmfm is installed
if ! command -v bmfm-targets-run &> /dev/null; then
    echo "[ERROR] bmfm-targets-run not found in env '$CONDA_ENV'."
    echo "  Install with: pip install 'biomed-multi-omic[dna]'"
    exit 1
fi

mkdir -p "$OUTPUT_DIRECTORY"
mkdir -p logs

# =============================================================================
# ── Fine-tuning run ───────────────────────────────────────────────────────────
# Hydra config: 02_finetune.yaml (in SCRIPT_DIR)
# This extends dna_finetune_train_and_test_config with:
#   - classification label (z_bin, 7 classes)
#   - cross-entropy loss
#   - 2048 token max_length (8kb bp windows ≈ 1600–2000 BPE tokens)
#   - max_epochs=5, val_check_interval=0.5
#
# A100 80GB + bf16-mixed: batch_size=4, accumulate=16 → effective batch=64.
# =============================================================================

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCHDYNAMO_DISABLE=1   # Triton/icx compiler crashes on CSD3

python "$PARENT_DIR/force_sdpa_wrapper.py" \
    --config-path "$SCRIPT_DIR" \
    -cn 02_finetune \
    input_directory="$INPUT_DIRECTORY" \
    working_dir="$OUTPUT_DIRECTORY" \
    "checkpoint='ibm-research/biomed.dna.snp.modernbert.113m.v1'"

EXIT_CODE=$?

echo "============================================================"
echo " Finished: $(date)"
echo " Exit code: $EXIT_CODE"
echo "============================================================"

if [ $EXIT_CODE -ne 0 ]; then
    echo "[ERROR] Fine-tuning failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

echo ""
echo "Outputs written to: $OUTPUT_DIRECTORY"
echo "  predictions.csv   — model predictions on test set"
echo "  checkpoints/      — saved model weights"
