#!/usr/bin/env bash
# =============================================================================
# train_agms3d_l4_loop.sh   --   AG-MS3D vanilla rescue, sequential on L4
# =============================================================================
# Runs the 15 (task, seed) AG-MS3D cells one at a time on the L4 cluster
# (dev-gpu-acs), where there's no SLURM and we can't parallelise across cells.
# Identical bundled rescue config to the CSD3 submit (commit eede09b):
#   --backbone vanilla --lr 1e-3 --label_smoothing 0.1
#
# Total wall-clock estimate: ~1-1.5h per cell on an L4 with batch=4 +
# 60 epochs, ES patience 15 → realistically 30-60 min per cell after ES
# kicks in. 15 cells x ~40 min = ~10h sequential.
#
# Pre-flight (once):
#   conda env list | grep mri      # confirm the `mri` env is present
#   ls $HOME/ADNI_SMRIPREP/cnn_inputs    # confirm rsync done
#   ls $HOME/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv
#   ls $HOME/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline
#
# Run (detached so SSH disconnect doesn't kill it):
#   mkdir -p $HOME/ADNI_MRI/agms3d_outputs
#   cd $HOME/ADNI_MRI/agms3d_outputs
#   nohup bash $HOME/Transformers_XAI/mri_pipeline/3d_cnn_vit/train_agms3d_l4_loop.sh \
#       > l4_loop_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#   tail -f l4_loop_*.log    # watch progress
#
# Resume:
#   Re-running the script auto-skips any cell that already has metrics.json
#   (same skip logic as the CSD3 submit). So if a cell crashes mid-run,
#   re-running picks up from the next cell. If a cell got partway through
#   training (best_model.pt exists but no metrics.json), the trainer's
#   --no_resume default = false means it'll AUTO-RESUME from last_checkpoint.pt.
#
# Single-cell smoke (test paths + env first, before kicking off all 15):
#   bash $HOME/Transformers_XAI/mri_pipeline/3d_cnn_vit/train_agms3d_l4_loop.sh \
#       --single T1c_binary 0     # T1c, seed 0 only
# =============================================================================

set -u
# Note: NO `set -e` -- we want one failed cell to NOT abort the loop.

# -- Env --
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mri

# -- L4 paths (shallow layout: data dirs live directly under $HOME) -----------
# Override any of these via env var on the sbatch command line if your
# layout differs, e.g.:
#   CNN_INPUTS_DIR=/scratch/cnn_inputs bash train_agms3d_l4_loop.sh
SCRIPT_DIR="${SCRIPT_DIR:-${HOME}/Transformers_XAI/mri_pipeline/3d_cnn_vit}"
CNN_INPUTS_DIR="${CNN_INPUTS_DIR:-${HOME}/ADNI_SMRIPREP/cnn_inputs}"
MATCHED_LABELS_CSV="${MATCHED_LABELS_CSV:-${HOME}/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv}"
SPLITS_DIR="${SPLITS_DIR:-${HOME}/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline}"
OUT_DIR="${OUT_DIR:-${HOME}/ADNI_MRI/agms3d_outputs}"

# -- Hyperparams (locked to the CSD3 rescue config) ---------------------------
BACKBONE="vanilla"
LR="1e-3"
LABEL_SMOOTHING="0.1"
NUM_WORKERS="${NUM_WORKERS:-4}"
WANDB_PROJECT="${WANDB_PROJECT:-agms3d_vanilla_rescue}"
export WANDB_MODE="${WANDB_MODE:-offline}"

# -- Sweep grid ---------------------------------------------------------------
TASKS=("T1_binary" "T1b_binary" "T1c_binary" "T1d_binary" "T2_multiclass")
SEEDS=(0 1 2)

# -- --single TASK SEED (smoke) overrides the loop ----------------------------
if [ "${1:-}" = "--single" ] && [ -n "${2:-}" ] && [ -n "${3:-}" ]; then
    TASKS=("$2")
    SEEDS=("$3")
    echo "[INFO] --single mode: running only TASK=$2 SEED=$3"
fi

mkdir -p "${OUT_DIR}"

# -- Pre-flight checks --------------------------------------------------------
fail_count=0
for path in "${CNN_INPUTS_DIR}" "${SPLITS_DIR}"; do
    if [ ! -d "${path}" ]; then
        echo "[ERROR] Required directory missing: ${path}"
        fail_count=$((fail_count + 1))
    fi
done
if [ ! -f "${MATCHED_LABELS_CSV}" ]; then
    echo "[ERROR] Matched labels CSV missing: ${MATCHED_LABELS_CSV}"
    fail_count=$((fail_count + 1))
fi
if [ ! -f "${SCRIPT_DIR}/train_agms3d.py" ]; then
    echo "[ERROR] Trainer script missing: ${SCRIPT_DIR}/train_agms3d.py"
    fail_count=$((fail_count + 1))
fi
if [ "${fail_count}" -gt 0 ]; then
    echo "Aborting pre-flight (${fail_count} missing inputs)."
    exit 1
fi

# -- GPU sanity ---------------------------------------------------------------
echo "============================================================"
echo "  AG-MS3D vanilla rescue -- sequential L4 loop"
echo "  Host       : $(hostname)"
echo "  Started    : $(date)"
echo "  CUDA       : $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo 'no GPU detected')"
echo "  conda env  : $(conda info --envs | awk '/\*/{print $1}')"
echo "  Out dir    : ${OUT_DIR}"
echo "  N cells    : $((${#TASKS[@]} * ${#SEEDS[@]}))  (${#TASKS[@]} tasks x ${#SEEDS[@]} seeds)"
echo "============================================================"

# -- Run cells sequentially ---------------------------------------------------
cell_idx=0
n_cells=$((${#TASKS[@]} * ${#SEEDS[@]}))
loop_t0=$(date +%s)

for TASK in "${TASKS[@]}"; do
    for SEED in "${SEEDS[@]}"; do
        cell_idx=$((cell_idx + 1))
        RUN_DIR="${OUT_DIR}/AGMS3DCNN_vanilla/${TASK}/seed_${SEED}"
        METRICS="${RUN_DIR}/metrics.json"

        echo
        echo "[$cell_idx/$n_cells] -------------------------------------------"
        echo "  TASK=${TASK}  SEED=${SEED}"
        echo "  Started: $(date)"
        echo "  Run dir: ${RUN_DIR}"

        # Skip if already done
        if [ -f "${METRICS}" ]; then
            echo "  metrics.json already exists -- skipping."
            continue
        fi

        # Resume hint
        if [ -f "${RUN_DIR}/last_checkpoint.pt" ]; then
            echo "  last_checkpoint.pt found -- will AUTO-RESUME."
        fi

        mkdir -p "${RUN_DIR}"
        # Offline wandb data lands inside the run dir alongside metrics
        export WANDB_DIR="${RUN_DIR}"

        cell_t0=$(date +%s)
        python "${SCRIPT_DIR}/train_agms3d.py" \
            --task                "${TASK}" \
            --seed                "${SEED}" \
            --backbone            "${BACKBONE}" \
            --lr                  "${LR}" \
            --label_smoothing     "${LABEL_SMOOTHING}" \
            --cnn_inputs_dir      "${CNN_INPUTS_DIR}" \
            --matched_labels_csv  "${MATCHED_LABELS_CSV}" \
            --data_dir            "${SPLITS_DIR}" \
            --out_dir             "${OUT_DIR}" \
            --wandb \
            --wandb_project       "${WANDB_PROJECT}" \
            --num_workers         "${NUM_WORKERS}"
        EXIT_CODE=$?
        cell_t1=$(date +%s)
        cell_min=$(( (cell_t1 - cell_t0) / 60 ))

        echo "  Cell finished in ${cell_min} min, exit=${EXIT_CODE}"
        if [ "${EXIT_CODE}" -ne 0 ]; then
            echo "  [WARN] non-zero exit -- continuing to next cell."
        fi
    done
done

loop_t1=$(date +%s)
total_min=$(( (loop_t1 - loop_t0) / 60 ))
echo
echo "============================================================"
echo "  All ${n_cells} cells done at $(date)"
echo "  Total wall-clock: ${total_min} min"
echo "  Outputs        : ${OUT_DIR}/AGMS3DCNN_vanilla/"
echo "  W&B sync       : wandb sync ${OUT_DIR}/AGMS3DCNN_vanilla/*/seed_*/wandb/offline-run-*"
echo "============================================================"
