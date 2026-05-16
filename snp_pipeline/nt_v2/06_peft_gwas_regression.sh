#!/bin/bash
#SBATCH -J ntv2_gwas_reg
#SBATCH -A BETHLEHEM-SL3-GPU
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ntv2_gwas_reg_%j.out
#SBATCH --error=logs/ntv2_gwas_reg_%j.err

set -euo pipefail
mkdir -p logs

# ── Environment ──────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate ntv2
export TORCHDYNAMO_DISABLE=1
export CUDA_LAUNCH_BLOCKING=0

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR="/rds/user/ec474/hpc-work/git/Transformers_XAI/snp_pipeline/nt_v2"
DATA_DIR="/rds/user/ec474/hpc-work/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom_combos/forward_and_reverse"
OUTDIR="/rds/user/ec474/hpc-work/results/ntv2_peft_gwas_regression"

echo "=== NT v2 GWAS Signed Z-Score Regression ==="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

python "$SCRIPT_DIR/06_peft_gwas_regression.py" \
    --data-dir "$DATA_DIR" \
    --outdir "$OUTDIR" \
    --epochs 5 \
    --lr 1e-5 \
    --batch-size 2 \
    --grad-accum 32 \
    --max-length 2048 \
    --device cuda

echo "=== Done ==="
