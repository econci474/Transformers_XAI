#!/bin/bash
# =============================================================================
# 05_debug_csd3.sh (gwas_signed_regression)
# =============================================================================
# Run 05_debug_dataloader.py on an A100 to reproduce the bf16 MSE
# collapse seen on the production combos run (job 29285536). The probe reports
# logits stats + MSE under fp32, CUDA bf16-autocast, and CUDA fp16-autocast.
# On L4 (Ada) all three give MSE≈151 — non-zero. On A100 (Ampere) we expect
# bf16 to produce NaN (silently swallowed by calculate_losses → loss=0).
#
# Submit from any directory; output lands in ./logs/.
# =============================================================================

#SBATCH -A LIO-CHARM-SL2-GPU
#SBATCH --job-name=bmfm_probe
#SBATCH --output=logs/bmfm_probe_%j.out
#SBATCH --error=logs/bmfm_probe_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH -p ampere

SCRIPT_DIR=/home/ec474/rds/hpc-work/Transformers_XAI/snp_pipeline/gwas_signed_regression
INPUT_DIRECTORY=/home/ec474/rds/hpc-work/ADNI_SNP/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/forward_and_reverse
RUNDIR=$(pwd)

echo "============================================================"
echo " BMFM bf16 collapse probe (A100)"
echo " Job ID    : $SLURM_JOB_ID"
echo " Node      : $SLURMD_NODENAME"
echo " GPU       : $CUDA_VISIBLE_DEVICES"
echo " Started   : $(date)"
echo "============================================================"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bmfm

export LD_PRELOAD=/home/ec474/glibc_compat.so
export TORCHDYNAMO_DISABLE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export PYTHONIOENCODING=utf-8

mkdir -p logs

# The yaml's hardcoded hydra.searchpath is the correct CSD3 path. The probe
# auto-detects CUDA and runs CPU fp32 + CUDA fp32/bf16/fp16 paths.
python "$SCRIPT_DIR/05_debug_dataloader.py" \
    --config-path "$SCRIPT_DIR" \
    -cn 04_finetune \
    input_directory="$INPUT_DIRECTORY" \
    working_dir="$RUNDIR" \
    "checkpoint=ibm-research/biomed.dna.ref.modernbert.113m.v1" \
    accelerator=cpu \
    precision=32-true \
    batch_size=1

EXIT=$?
echo "============================================================"
echo " Finished: $(date)  exit=$EXIT"
echo "============================================================"
exit $EXIT
