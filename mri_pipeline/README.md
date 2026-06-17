# `mri_pipeline/` — structural-MRI models for Alzheimer's prediction (ADNI)

Fine-tunes and benchmarks several 3D imaging encoders on sMRIprep-preprocessed
T1w MNI volumes, for the shared diagnostic/prognostic tasks (below). Every script
is numbered and runs in sequence; later scripts read earlier scripts' outputs from
canonical paths. There is no orchestrator — re-running a step requires its
predecessors' outputs to exist.

All data, checkpoints and figures live **outside** the repo (`D:\ADNI_BIDS_project\`
locally, `/home/ec474/rds/hpc-work/` on CSD3). `mri_pipeline/outputs/` is gitignored.
`mri_pipeline/tools/` (gitignored) holds local-only provenance helpers (W&B sync,
cached query CSVs) that are not part of the pipeline examiners need to read.

## Tasks (shared across all modality pipelines)
- **T1d** — sMCI vs pMCI (MCI conversion to AD).
- **T1e** — sCN vs pCN (CN progression to MCI/AD).
- **T2**  — 3-way CN / MCI / AD.
- **T3a/b/c** — conversion within ≤3 / ≤5 / ≤7 years (baseline-anchored binaries).
- **T4**  — AD-conversion horizon (3-class; converter-only cohort).
  (T1/T1b/T1c are auxiliary diagnostic binaries used during development.)

Splits are 80/10/10 across seeds 0/1/2 with test subjects held consistent across
modalities. Subject/session exclusions come from `bidsification/exclusions.py`
(the single source of truth), imported by the trainers.

## Model arms
| Arm | Encoder | Where |
|---|---|---|
| ViT-MAE75 | MAE-pretrained 3D ViT-B/16 (`_vit_recipe/`) | root `03/04`, `03d` |
| BrainDINO | DINOv3 ViT-B/16 (per-slice) | `brain_dino/` |
| BrainMVP  | UniFormer-Small | `brain_mvp/` |
| AG-MS3D-CNN | attention-gated multi-scale 3D CNN | `3d_cnn_vit/` |
| Spasov 3D-CNN | vanilla / separable 3D CNN baselines | `3d_conv_net/` |
| Cached-head | frozen-encoder embeddings + linear/MLP head | `cached_head_sweep/` |

## Pipeline order

**Inputs & MRI↔clinical matching**
- `03_prepare_ViT.py` — sMRIprep MNI T1w → ViT inputs (1.75 mm, 128³, RAS, nonzero z-score, foreground crop).
- `03b_match_mri_to_clinical.py` — pair each scan with its clinical row via date-derived VISCODE (with a strategy audit). Feeds the clinical coverage tables.
- `03c_match_mri_to_clinical_viscode2.py` — **canonical** MRI↔clinical match on ADNI's native VISCODE2. Its output (extended + post-exclusion) is the training-label CSV every finetune/sweep reads.
- `03d_extract_vit_mae_embeddings.py` — one-time frozen ViT-MAE75 pooled embeddings (for the cached-head sweep).

**Per-encoder input prep / finetune / embeddings** (each subfolder is self-contained)
- `brain_dino/01_*` prep → `02_supervised_finetuning_BrainDINO.py` (frozen / full_ft / LoRA) → `03_extract_braindino_embeddings.py`.
- `brain_mvp/03_prepare_BrainMVP.py` → `04_supervised_finetuning_BrainMVP.py` → `05_extract_brainmvp_full_ft_embeddings.py`; `06_gradcam_T1d.py` + `07_render_gradcam_T1d.py` (T1d saliency on the real MNI T1w).
- `3d_conv_net/00_prepare_CNN_inputs.py` → `train_3dcnn.py` (Spasov vanilla/separable).
- `3d_cnn_vit/train_agms3d.py` (AG-MS3D-CNN).
- `04_supervised_finetuning_ViT.py` — the ViT trainer **and shared library**: its task config and label/metric helpers are imported (by exact path, via importlib) by the BrainDINO / BrainMVP / 3D-CNN / cached-head trainers, so its filename is intentionally stable.

**Cached-head HP sweep** (`cached_head_sweep/`)
- `04_head_finetune_from_embeddings.py` — train a head on cached embeddings over the HP grid.
- `05b_select_best_hp_per_task.py` — pick the HP winner per (model, task) by mean val balanced-accuracy.
- `04_finalize_winner.py` — retrain the winner to materialise its checkpoint.
- `04c_render_sweep_table.py` — render sweep result tables.

**Hyper-parameter search & metric collection**
- `08_hp_optuna_full_ft.py` — Optuna full-finetune HP search (BrainMVP / ViT-MAE).
- `09_collect_val_test_metrics.py` — recompute val + test metrics from saved checkpoints.

**Aggregation & thesis tables**
- `05_aggregate_mri_results.py` — aggregate every run's `metrics.json` (mean ± std over seeds; incl. true val macro-F1).
- `05_extract_cached_head_probs.py` / `05_extract_embeddings_braindino_frozen.py` — per-scan probabilities/embeddings for late fusion.
- `06_render_cross_model_table.py`, `06b_widen_cross_model_tables.py`, `06b_export_hp_provenance.py`, `06c_render_styled_cross_model.py`, `07_render_prognostic_conversion_tables.py`, `_render_mri_t3_val_macrof1.py` — the cross-model comparison tables (diagnostic + prognostic) used in the thesis.

**Shared utilities** — `_vit_recipe/` (3D ViT architecture + MAE checkpoint loading), `tests/test_vit_preprocessing.py`.

## Compute & conventions
`*_submit_csd3.sh` are CSD3 SLURM wrappers around the matching `.py`; see the root
`tool.md` for account/partition/wall-time specifics. Scripts carry a LOCAL/HPC
path block at the top — switch it when porting. Run from the output folder so
relative SLURM logs land beside the checkpoints.
