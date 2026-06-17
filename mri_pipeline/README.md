# `mri_pipeline/` — structural-MRI models for Alzheimer's prediction (ADNI)

Fine-tunes and benchmarks several 3D imaging encoders on sMRIprep-preprocessed T1w MNI volumes for the
shared diagnostic/prognostic tasks (below). Scripts are numbered and run in sequence; later steps read
earlier steps' outputs from canonical paths. There is no orchestrator — re-running a step requires its
predecessors' outputs to exist.

All data, checkpoints and figures live **outside** the repo (`D:\ADNI_BIDS_project\` locally,
`/home/ec474/rds/hpc-work/` on CSD3). `mri_pipeline/outputs/` and `mri_pipeline/tools/` (local-only
provenance helpers: W&B sync scripts, cached query CSVs) are gitignored.

## Tasks (shared across all modality pipelines)
- **T1d** — sMCI vs pMCI (MCI conversion to AD).
- **T1e** — sCN vs pCN (CN progression to MCI/AD).
- **T2**  — 3-way CN / MCI / AD.
- **T3a/b/c** — conversion within ≤3 / ≤5 / ≤7 years (baseline-anchored binaries).
- **T4**  — AD-conversion horizon (3-class; converter-only cohort).
  (T1/T1b/T1c are auxiliary diagnostic binaries used during development.)

Splits are 80/10/10 across seeds 0/1/2 with test subjects held consistent across modalities.
Subject/session exclusions come from `bidsification/exclusions.py` (the single source of truth).

## Layout
```
mri_pipeline/
  03b_match_mri_to_clinical.py           MRI<->clinical match (date-derived VISCODE; clinical-coverage)
  03c_match_mri_to_clinical_viscode2.py  MRI<->clinical match (ADNI VISCODE2) — CANONICAL training labels
  04_supervised_finetuning.py            shared 3D-ViT trainer + library (task config / label & metric
                                         helpers) — imported by exact path by the other arms' trainers
  05_extract_cached_head_probs.py        per-scan probabilities for late fusion (cross-model)
  vit_mae/      ViT-MAE arm — 03_prepare_ViT.py, 03d_extract_vit_mae_embeddings.py, ViT 04* submit scripts
  brain_dino/   BrainDINO arm — 01 prep, 02 finetune (frozen/full_ft/LoRA), 03 extract embeddings,
                05_extract_embeddings_braindino_frozen.py (head -> probs+embeddings for fusion)
  brain_mvp/    BrainMVP arm — 03 prep, 04 finetune, 05 extract; 06/07 T1d Grad-CAM (real-MNI overlay)
  agms3d/       AG-MS3D-CNN arm — AGMS3DCNN.py + train_agms3d.py
  3d_conv_net/  Spasov 3D-CNN baselines — 00 prep, 3DCNN/3DSCNN, train_3dcnn.py
  cached_head_sweep/  frozen-encoder head HP sweep — 04 train, 05b select winner, 04_finalize_winner,
                      04c_render_sweep_table
  tables/       cross-model result tables — 05_aggregate_mri_results + 06/06b/06c/07 renderers + macro-F1
  optuna_sweep/ exploratory Optuna full-finetune HP search (ViT-MAE / BrainMVP) — not used for final models
  baseline_vit/ from-scratch ViT control submits (ViT-tiny scratch, ViT high-LR)
  _vit_recipe/  shared 3D ViT-B/16 architecture + MAE checkpoint loading
  tests/        test_vit_preprocessing.py (ViT preprocessing + shared-trainer input-gate smoke)
```

## Pipeline order
1. **Inputs & matching** — `vit_mae/03_prepare_ViT.py` builds ViT inputs (1.75 mm, 128³, RAS, nonzero
   z-score, foreground crop); `03b`/`03c` pair scans to clinical rows (the `03c` VISCODE2 output, extended
   + post-exclusion, is the training-label CSV every finetune/sweep reads); `vit_mae/03d_extract` caches
   frozen ViT-MAE embeddings. Each other arm has its own prep in its folder.
2. **Fine-tune** — `04_supervised_finetuning.py` (ViT-MAE) is the shared trainer; the BrainDINO / BrainMVP
   / AG-MS3D / 3D-CNN / cached-head trainers import its task config + label/metric helpers via importlib.
3. **Cached-head sweep** — `cached_head_sweep/`: train heads on cached embeddings → `05b_select_best_hp_per_task`
   → `04_finalize_winner` (materialise the winner checkpoint).
4. **HP search** — `08_hp_optuna_full_ft.py` (Optuna full-finetune search).
5. **Aggregate & render** — `tables/05_aggregate_mri_results.py` aggregates every `metrics.json`
   (mean ± std over seeds, true val macro-F1); `tables/06*/07*` render the thesis cross-model tables.
6. **Explainability** — `brain_mvp/06_gradcam_T1d.py` + `07_render_gradcam_T1d.py` (T1d saliency
   back-projected onto the real MNI T1w).

## Conventions
`*_submit_csd3.sh` are CSD3 SLURM wrappers kept beside the `.py` they launch (the ViT submits in `vit_mae/`
point `SCRIPT_DIR` at `mri_pipeline/` since they invoke the shared trainer there). Scripts carry a LOCAL/HPC
path block at the top — switch it when porting. See the root `tool.md` for account/partition/wall-time
specifics. Run a submit from the output folder so relative SLURM logs land beside the checkpoints.
