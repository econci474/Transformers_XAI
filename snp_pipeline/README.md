# SNP Pipeline — README

## Overview

This pipeline takes raw ADNI genotype data (GRCh37/hg19), performs quality
control, runs GWAS, and prepares fine-tuning datasets for the
**BMFM-DNA-SNP** foundation model
(`ibm-research/biomed.dna.snp.modernbert.113m.v1`) to predict Alzheimer's
disease GWAS z-scores from allele-specific DNA sequences.

---

## Directory Layout

```
snp_pipeline/
│
├── 01_qc.sh                                  PLINK QC (MAF, missingness, HWE)
├── 02_pca.py                                 Principal component analysis
├── 03_prepare_pheno_cn_vs_ad.py              Phenotype file: CN vs AD
├── 03_prepare_pheno_cn_vs_emci.py            Phenotype file: CN vs eMCI
│
├── 04a_run_gwas_cn_vs_ad.sh                  GWAS logistic regression (CN vs AD)
├── 04a_run_gwas_cn_vs_emci.sh                GWAS logistic regression (CN vs eMCI)
├── 04b_parse_gwas_results.py                 Parse GWAS results; compute FDR thresholds
├── 04c_gwas_plots.py                         Manhattan + QQ plots
├── 04d_gwas_annotate.py                      Annotate GWAS hits
│
├── 05a_update_snp_ids.py                     Convert kgp*/probe IDs → dbSNP rsIDs
├── 05b_resolve_preexisting_dups.py           Resolve duplicate variant IDs
├── 05c_validate_rsid.py                      Validate rsID format
├── 05d_inspect_duplicate_rsid.py             Inspect remaining duplicates
│
├── 06_liftover_hg19_to_hg38.py              Liftover ADNI genotypes GRCh37 → GRCh38
├── 07_ld_07_comparison_script.py             LD comparison / clumping
│
├── 08a_liftover_hg19_to_hg38_wightman_2021.py  Liftover Wightman 2021 GWAS → GRCh38
├── 08b_prepare_bmfm_labels_from_external_gwas.py  Merge external GWAS → label TSV
├── 08c_prepare_bmfm_gwas_regression_inputs.py     1,001 bp window inputs (EA/OA pairs)
├── 08d_augment_bmfm_gwas_regression.py            Reverse complement augmentation (×2)
├── 08e_augment_bmfm_gwas_regression_by_chrom.py   8 kb chromosome-context windows
│
├── 09_bmfm_gwas_regression_finetuning.sh    SLURM job submission script (CSD3 ampere)
├── 09_bmfm_gwas_regression_finetuning.yaml  Hydra fine-tuning configuration
│
└── 10_prepare_bmfm_dna_snp_inputs.py        BMFM-DNA inference inputs (classification)
```

---

## Pipeline Flow

```
Raw ADNI genotypes (hg19 PLINK binary)
         │
    01_qc.sh          — MAF ≥ 0.01, missingness ≤ 0.05, HWE p > 1e-6
         │
    05a–05d           — rsID update, duplicate resolution (PLINK --update-name)
         │
    06_liftover       — PLINK binary files lifted to GRCh38
         │
    02_pca.py         — PCA on QC'd GRCh38 genotypes (for covariate file)
         │
    03_prepare_pheno  — PLINK phenotype + covariate files (age, sex, APOE, PCs)
         │
    04a_run_gwas      — PLINK2 logistic regression (CN vs AD / CN vs eMCI)
         │
    04b_parse_gwas    — Extract GW-sig hits (p < 5×10⁻⁸); compute FDR
         │
    08a_liftover_wightman — Lift published Wightman 2021 summary stats to GRCh38
         │
    08b_prepare_labels    — Merge ADNI GWAS + Wightman + Bellenguez 2022
                            → external_gwas_labels.tsv (label=1 GW-sig, label=0 null)
         │
    08c_prepare_inputs    — Extract 1,001 bp EA/OA sequences from GRCh38 FASTA
                            → train/dev/test CSVs (sequence, z_score, seq_id)
         │
    08d_augment           — Add reverse complement strands (×2 sequences per SNP)
         │
    08e_augment_by_chrom  — Build 8 kb multi-SNP chromosome-context windows
                            → forward/ and forward_and_reverse/ CSVs
         │
    09_bmfm_finetuning.sh — Fine-tune BMFM-DNA-SNP on Cambridge CSD3 (A100 GPU)
```

---

## Script Usage

### Conda environments

- **`snp`** — data preparation (local Windows machine)
- **`bmfm`** — model fine-tuning (CSD3 HPC, Linux)

Activate before running each group of scripts:
```bash
conda activate snp    # for 08a–08e, local
conda activate bmfm   # for 09, on HPC
```

---

### 08b — Build external GWAS label TSV

```bash
python 08b_prepare_bmfm_labels_from_external_gwas.py
```

**Inputs:**
- `D:/ADNI_SNP_Omni2.5M_20140220/gwas_results/` (ADNI GWAS output)
- Wightman 2021 summary stats (GRCh38 lifted)
- Bellenguez 2022 summary stats

**Output:**
- `D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_labels.tsv`
  Columns: `SNP, CHR, BP, EA, OA, REF, z_score, label, source`

---

### 08c — Prepare 1,001 bp regression inputs

```bash
# Default (500 bp flanks, 80/10/10 split, seed 42)
python 08c_prepare_bmfm_gwas_regression_inputs.py

# Custom flank size
python 08c_prepare_bmfm_gwas_regression_inputs.py --flank 200

# Quick test with 50 SNPs
python 08c_prepare_bmfm_gwas_regression_inputs.py --max-snps 50
```

**Inputs:**
- `external_gwas_labels.tsv` (label=1 SNPs only)
- `Homo_sapiens.GRCh38.dna.primary_assembly.fa`

**Outputs** → `bmfm_gwas_signed_regression_without_ukb/`:
- `train.csv`, `dev.csv`, `test.csv` — columns: `sequence, z_score, seq_id`
- `metadata.tsv` — full record per sequence
- `dataset_summary.txt`
- `finetune_config.yaml`

**Split logic:** Random shuffle at the SNP level (seed=42).
EA and OA pairs for the same SNP always land in the same split.
Default: 80% train / 10% dev / 10% test.

**Z-score convention:**
- EA sequence → `+z` (risk allele)
- OA sequence → `−z` (other allele)

---

### 08d — Reverse complement augmentation (×2)

```bash
# Dry run (print counts only)
python 08d_augment_bmfm_gwas_regression.py --check

# Full run
python 08d_augment_bmfm_gwas_regression.py
```

**Input:** `bmfm_gwas_signed_regression_without_ukb/`

**Output** → `bmfm_gwas_signed_regression_without_ukb_augmented/`:
- 4 rows per SNP: `rs*_EA_f`, `rs*_EA_r`, `rs*_OA_f`, `rs*_OA_r`
- RC gets the **same z-score** as the forward sequence (same biological effect)

---

### 08e — Chromosome-context 8 kb windows

```bash
# Dry run
python 08e_augment_bmfm_gwas_regression_by_chrom.py --check

# Full run (default: 8,192 bp windows, 50 bp min flank)
python 08e_augment_bmfm_gwas_regression_by_chrom.py

# Smaller windows
python 08e_augment_bmfm_gwas_regression_by_chrom.py --max-len 4096
```

**Design:**
- GW-sig SNPs grouped by chromosome, packed greedily into ≤8,192 bp windows
- For each window: one sequence with **all SNPs → EA** (`z = Σ z_k`),
  one with **all SNPs → OA** (`z = −Σ z_k`)
- All other genomic positions remain GRCh38 reference
- Minimum flank: `--min-flank 50` (default) ensures ≥50 bp on each side

**Outputs** → `bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom/`:
- `forward/train.csv`, `forward/dev.csv`, `forward/test.csv`
- `forward_and_reverse/train.csv`, `forward_and_reverse/dev.csv`, `forward_and_reverse/test.csv`

**Split logic:** Inherited from 08c. For multi-SNP windows with mixed
splits, majority vote determines the assigned split.

**Metadata columns** (right of `seq_id`):
`window, chrom, n_snps_in_window, win_len_bp,`
`rsid_snp1, position_hg38_snp1, allele_type_snp1, z_score_snp1,`
`rsid_snp2, position_hg38_snp2, allele_type_snp2, z_score_snp2, …`

---

### 09 — Fine-tune on CSD3 (HPC)

```bash
# On the HPC login node — transfer data first
scp -r D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb \
    ec474@login-cpu.hpc.cam.ac.uk:/rds/hpc-work/ADNI_SNP/

# Create log directory
mkdir -p /rds/hpc-work/ADNI_SNP/logs

# Submit job
sbatch /rds/hpc-work/Transformers_XAI/snp_pipeline/09_bmfm_gwas_regression_finetuning.sh

# Monitor
squeue -u ec474
tail -f /rds/hpc-work/ADNI_SNP/logs/bmfm_gwas_regression_<JOBID>.out
```

**Key SLURM settings** (`09_bmfm_gwas_regression_finetuning.sh`):
- Partition: `ampere` (A100 GPUs)
- Account: `COMPUTERLAB-SL2-GPU`
- Memory: 32 GB RAM, 1 GPU
- Wall time: 2 hours

**Hydra configuration** (`09_bmfm_gwas_regression_finetuning.yaml`):
- Task: signed z-score regression (MSE loss)
- Tokenizer: `ref2vec` (character-level, compatible with pre-trained weights)
- Max sequence length: 1,001 bp (08c) or 8,192 bp (08e)
- Checkpoint: `ibm-research/biomed.dna.snp.modernbert.113m.v1`

---

## Data Paths (Local Windows)

| Resource | Path |
|---|---|
| Raw genotypes (hg19) | `D:/ADNI_SNP_Omni2.5M_20140220/` |
| GRCh38 reference FASTA | `D:/ADNI_SNP_Omni2.5M_20140220/Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| GWAS label TSV | `D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/external_gwas_labels.tsv` |
| 08c output (1,001 bp) | `D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb/` |
| 08d output (+ RC) | `D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented/` |
| 08e output (8 kb) | `D:/ADNI_SNP_Omni2.5M_20140220/bmfm_inputs/bmfm_gwas_signed_regression_without_ukb_augmented_by_chrom/` |

## Data Paths (HPC RDS)

| Resource | Path |
|---|---|
| Input CSVs | `/rds/hpc-work/ADNI_SNP/bmfm_gwas_signed_regression_without_ukb/` |
| Model output | `/rds/hpc-work/ADNI_SNP/bmfm_gwas_regression_output/` |
| SLURM logs | `/rds/hpc-work/ADNI_SNP/logs/` |
| Pipeline scripts | `/rds/hpc-work/Transformers_XAI/snp_pipeline/` |

---

## Environment Requirements

### `snp` — Data preparation (local, Python 3.10+)

Used for: scripts 01–08e (all data preparation steps).

```bash
conda create -n snp python=3.10
conda activate snp
```

**Core packages:**

```
pandas>=2.0
numpy>=1.24
scipy
matplotlib
seaborn
pyfaidx          # FASTA indexing and random access
biopython        # sequence utilities (optional)
pyranges         # genomic interval operations (optional, used in liftover)
pyliftover       # coordinate liftover (hg19 → hg38)
requests         # downloading chain files
tqdm
```

Install:
```bash
pip install pandas numpy scipy matplotlib seaborn pyfaidx biopython \
            pyranges pyliftover requests tqdm
```

**External tools** (must be on PATH):
- `plink` / `plink2` — QC, GWAS, liftover
- `bcftools` (optional, for VCF operations)

---

### `bmfm` — Model fine-tuning (HPC, Python 3.10/3.12)

Used for: script 09 (SLURM job on CSD3 ampere partition).

```bash
conda create -n bmfm python=3.10
conda activate bmfm
```

**Core packages:**

```
torch>=2.1           # PyTorch with CUDA support
transformers>=4.37   # HuggingFace (ModernBERT backbone)
biomed-multi-omic    # IBM BMFM framework (bmfm-targets-run CLI)
hydra-core>=1.3      # configuration management
omegaconf
lightning>=2.1       # PyTorch Lightning (training loop)
huggingface-hub      # model download
pybind11             # required build dependency
```

Install on HPC (within Conda environment, after loading CUDA modules):
```bash
module load cuda/12.4 cudnn/9.2_cuda-12.4

# pybind11 must be installed before biomed-multi-omic
CFLAGS="-I$CONDA_PREFIX/include/python3.10" \
CXXFLAGS="-I$CONDA_PREFIX/include/python3.10" \
pip install pybind11

pip install torch --index-url https://download.pytorch.org/whl/cu124
pip install biomed-multi-omic
pip install hydra-core omegaconf lightning huggingface-hub
```

**Pre-cache the model weights** (run once on a login node before submitting):
```python
from huggingface_hub import snapshot_download
snapshot_download("ibm-research/biomed.dna.snp.modernbert.113m.v1")
```

---

## Key Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Tokenizer | `ref2vec` | Character-level ACGT; compatible with pre-trained ModernBERT weights |
| Regression target | Signed z-score (+z for EA, −z for OA) | Encodes effect direction; MSE loss |
| Window size (08c) | 1,001 bp (500 bp each side) | Standard flanking for TF motif context |
| Window size (08e) | 8,192 bp (ModernBERT max) | Maximum genomic context |
| Z-score (08e multi-SNP) | Σ z_k for all SNPs in window | Additive polygenic assumption |
| Train/dev/test split | 80 / 10 / 10, seed=42, SNP-level | EA+OA pairs always in same split |
| RC augmentation | Same z-score as forward | Same biological effect on complementary strand |
