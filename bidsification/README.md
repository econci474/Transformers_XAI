# ADNI BIDSification Pipeline

Scripts to convert raw ADNI NIfTI data into a [BIDS 1.7.0](https://bids-specification.readthedocs.io) compliant dataset.

## Overview

This pipeline takes pre-processed N3m-corrected NIfTI files from `sourcedata/ADNI/` and organises them into `bids/` following the Brain Imaging Data Structure specification.

## Pipeline Scripts

| Script | Description |
|--------|-------------|
| `01_build_session_map.py` | Parses `MRIQUALITY.csv` + `image_series_mapping.csv` to generate a comprehensive `metadata/session_map.csv` and `metadata/ses_to_visit_code.csv` |
| `02_build_qc_selection.py` | Ranks T1w scans by QC score (MAYOADIRL) + series preference (GradWarp/B1/N3) → `metadata/scan_selection.csv` |
| `03_copy_niftis_to_bids.py` | Copies + gzip-compresses NIfTIs into `bids/sub-*/ses-*/anat/` (with `run-1`/`run-2` labels when multiple scans exist) |
| `04_generate_json_sidecars.py` | Generates BIDS T1w JSON sidecars from MRIQUALITY scanner metadata |
| `05_build_participants_tsv.py` | Builds `bids/participants.tsv` + `participants.json` from ADNIMERGE + APOERES |
| `06_populate_phenotype.py` | Copies key clinical CSVs into `bids/phenotype/` as TSVs + JSON sidecars |
| `07_generate_scans_tsv.py` | Generates per-session `*_scans.tsv` files with acquisition dates |

## Running the Pipeline

Run scripts **in order** from the project root:

```bash
cd D:\ADNI_BIDS_project

python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\01_build_session_map.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\02_build_qc_selection.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\03_copy_niftis_to_bids.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\04_generate_json_sidecars.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\05_build_participants_tsv.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\06_populate_phenotype.py
python C:\Users\elena\iCloudDrive\Desktop\ACS_MPhil\Thesis\git\Transformers_XAI\bidsification\07_generate_scans_tsv.py
```

## Key Design Decisions

### Session Labels
BIDS sessions use **date-based labels** (`ses-YYYYMMDD`), e.g. `ses-20110602`.  
The correspondence to ADNI visit codes (`bl`, `m06`, etc.) is preserved in:
```
D:\ADNI_BIDS_project\metadata\ses_to_visit_code.csv
```

### Multi-Scan Handling
When a subject has >1 T1w scan per session, **both** are kept:
- `*_run-1_T1w.nii.gz` — preferred scan (highest QC score + best series type)
- `*_run-2_T1w.nii.gz` — secondary scan

The **recommended** scan per session is documented in:
```
D:\ADNI_BIDS_project\metadata\scan_selection.csv
```

### Series Preference Ranking (within a session)
Lower rank = more preferred:

| Rank | Series Type |
|------|-------------|
| 0 | GradWarp + B1 + N3m corrected |
| 1 | GradWarp + N3m |
| 2 | N3m only (`MT1__N3m`) |
| 3 | GradWarp only |
| 4 | Accelerated MPRAGE (ADNI3) |
| 5 | Standard MPRAGE / IR-FSPGR |

### NIfTI Source
All NIfTIs come from `sourcedata/ADNI/<SubjectID>/MT1__N3m/`.  
These are ADNI's N3-corrected images (no re-conversion from DICOM needed).

## Output Structure

```
bids/
├── dataset_description.json
├── participants.tsv
├── participants.json
├── phenotype/
│   ├── adnimerge.tsv / .json
│   ├── cdr.tsv / .json
│   ├── adas.tsv / .json
│   ├── mmse.tsv / .json
│   ├── apoe.tsv / .json
│   └── study_arm.tsv / .json
└── sub-<label>/
    └── ses-<YYYYMMDD>/
        ├── sub-<label>_ses-<YYYYMMDD>_scans.tsv
        └── anat/
            ├── sub-<label>_ses-<YYYYMMDD>[_run-N]_T1w.nii.gz
            └── sub-<label>_ses-<YYYYMMDD>[_run-N]_T1w.json
```

## Intermediate Metadata Files

| File | Description |
|------|-------------|
| `metadata/session_map.csv` | Full record of all T1w images (source → BIDS mapping) |
| `metadata/ses_to_visit_code.csv` | Correspondence: `bids_ses` ↔ ADNI `VISCODE2` |
| `metadata/scan_selection.csv` | Per-session recommended scan + all alternatives |
| `metadata/03_nifti_copy.log` | Log from script 03 (copy results) |

## Requirements

```
pip install pandas numpy
```

Python ≥ 3.8. No DICOM tools needed (NIfTIs already available in sourcedata/).
