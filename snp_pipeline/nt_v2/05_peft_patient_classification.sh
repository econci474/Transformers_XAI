#!/bin/bash
#SBATCH -J ntv2_patient_clf
#SBATCH -A BETHLEHEM-SL3-GPU
#SBATCH -p ampere
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/ntv2_patient_clf_%j.out
#SBATCH --error=logs/ntv2_patient_clf_%j.err

set -euo pipefail
mkdir -p logs

# ── Environment ──────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate ntv2
export TORCHDYNAMO_DISABLE=1
export CUDA_LAUNCH_BLOCKING=0

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR="/rds/user/ec474/hpc-work/git/Transformers_XAI/snp_pipeline/nt_v2"
SEQ_CSV="/rds/user/ec474/hpc-work/bmfm_inputs/patient_window_sequences/all_patients.csv"
LABELS_DIR="/rds/user/ec474/hpc-work/clinical/no_cdr_stratified_ever_convert/tabular/baseline"
OUTDIR="/rds/user/ec474/hpc-work/results/ntv2_peft_patient_clf"

SEED=${1:-0}

echo "=== NT v2 Patient Classification (seed=$SEED) ==="
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

python "$SCRIPT_DIR/05_peft_patient_classification.py" \
    --sequences "$SEQ_CSV" \
    --labels-dir "$LABELS_DIR" \
    --seed "$SEED" \
    --outdir "$OUTDIR" \
    --epochs 10 \
    --lr 1e-4 \
    --batch-size 2 \
    --grad-accum 16 \
    --max-length 2048 \
    --device cuda

echo "=== Done ==="
