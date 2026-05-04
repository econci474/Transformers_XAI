#!/bin/bash
# =============================================================================
# bmfm_gwas_regression_finetuning.sh
# =============================================================================
# SLURM job script for BMFM-DNA-SNP fine-tuning on signed GWAS z-score
# regression.
#
# Model  : ibm-research/biomed.dna.snp.modernbert.113m.v1  (ModernBERT, 113M)
# Task   : Regression — predict GWAS z-score from allele-specific DNA sequence
# Data   : 430 SNPs × 2 alleles = 860 rows  (688 train / 86 dev / 86 test)
# Input  : bmfm_gwas_signed_regression_without_ukb/ (train/dev/test.csv)
# Config : bmfm_gwas_regression_finetuning.yaml
#
# GPU estimate:
#   ModernBERT 113M @ fp16, 1001 bp sequences, batch=16 → ~4 GB VRAM.
#   Full run (50 epochs, 43 steps/epoch) completes in ~10–20 min on 1× A100.
#   A V100 (32 GB) or A30 (24 GB) are both sufficient; RTX 3090 also works.
#
# Usage:
#   sbatch bmfm_gwas_regression_finetuning.sh
#
# Adjust CONDA_ENV, INPUT_DIRECTORY, OUTPUT_DIRECTORY and SCRIPT_DIR below.
# =============================================================================


#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmfm_gwas_combos
#SBATCH --output=logs/bmfm_gwas_combos_%j.out
#SBATCH --error=logs/bmfm_gwas_combos_%j.err
#SBATCH --time=03:00:00          # 2 epochs on A100 with 08f (~11k train seqs) ≈ 1–2 h
#SBATCH --nodes=1                # Cambridge HPC policy: must be explicit for GPU jobs
#SBATCH --gres=gpu:1             # 1 A100 40 GB
#SBATCH --cpus-per-task=3        # CSD3 max: 3 CPUs per GPU
#SBATCH -p ampere       # Cambridge CSD3 GPU partition (A100s)

# =============================================================================
# ── User-configurable paths ───────────────────────────────────────────────────
# =============================================================================

CONDA_ENV="bmfm"       # conda env with bmfm_targets installed
                        # install: pip install biomed-multi-omic[dna]

# Combinatorial EA/OA by-chrom dataset (08f, generated locally then scp'd to CSD3):
# scp -r D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos \
#     ec474@login-cpu.hpc.cam.ac.uk:~/ADNI_SNP/
INPUT_DIRECTORY="/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos"

# Where to write checkpoints, predictions, logs:
OUTPUT_DIRECTORY="/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_regression_output/bmfm_gwas_combos_lr1e6_2ep"

# Directory containing 09_bmfm_gwas_regression_finetuning.yaml
# NOTE: Cannot use BASH_SOURCE[0] here — SLURM copies the .sh to a spool
# directory (/var/spool/slurm/...) so dirname resolves to the wrong path.
# Hardcode to the actual git clone location on CSD3:
SCRIPT_DIR="/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline"

# =============================================================================
# ── Environment ───────────────────────────────────────────────────────────────
# =============================================================================

echo "============================================================"
echo " BMFM-DNA-SNP GWAS Regression Fine-tuning"
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

# Verify bmfm is installed
if ! command -v bmfm-targets-run &> /dev/null; then
    echo "[ERROR] bmfm-targets-run not found in env '$CONDA_ENV'."
    echo "  Install with: pip install 'biomed-multi-omic[dna]'"
    exit 1
fi

mkdir -p "$OUTPUT_DIRECTORY"
mkdir -p logs

# =============================================================================
# =============================================================================
# ── Fine-tuning run ───────────────────────────────────────────────────────────
# Hydra config: 09_bmfm_gwas_regression_finetuning.yaml (in SCRIPT_DIR)
# This extends dna_finetune_train_and_test_config with:
#   - regression label (z_score, is_regression_label=true)
#   - MSE loss
#   - 8192 bp max_length (chromosome-level 8kb windows)
#   - max_epochs=2, val_check_interval=0.5 (BMFM fine-tuning recommendation)
#
# A100 override: batch_size=4, accumulate=16 → effective batch=64 (same as L4
#   run at batch_size=1, accumulate=64, but ~4× faster per epoch on A100).
# =============================================================================

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCHDYNAMO_DISABLE=1   # Triton/icx compiler crashes on CSD3 AMD nodes

bmfm-targets-run \
    --config-path "$SCRIPT_DIR" \
    -cn 09_bmfm_gwas_regression_finetuning \
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
echo ""
echo "To run inference on new sequences:"
echo "  bmfm-targets-run -cn dna_predict \\"
echo "      input_directory=\$INPUT_DIRECTORY \\"
echo "      working_dir=\$OUTPUT_DIRECTORY \\"
echo "      checkpoint=\$OUTPUT_DIRECTORY/checkpoints/<best>.ckpt"
