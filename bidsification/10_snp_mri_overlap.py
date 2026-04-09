"""
Script 10 — SNP + MRI Subject Overlap
======================================
Identifies subjects in the BIDS dataset (3T T1w MRI) who also have
SNP genotype data in PLINK binary format, and saves:

  bids/genotype/subjects_with_snp_and_mri.tsv   — full metadata per subject
  bids/genotype/snp_subjects_with_mri.fam        — PLINK keep-list for genotype subsetting

Inputs:
  D:/ADNI_SNP_Omni2.5M_20140220/WGS_Omni25_BIN_wo_ConsentsIssues.fam
  D:/ADNI_BIDS_project/bids/participants.tsv
"""

import pandas as pd
import os

# ── Paths ─────────────────────────────────────────────────────────────────────
FAM_FILE  = r'D:\ADNI_SNP_Omni2.5M_20140220\WGS_Omni25_BIN_wo_ConsentsIssues.fam'
BIDS_DIR  = r'D:\ADNI_BIDS_project\bids'
OUT_DIR   = os.path.join(BIDS_DIR, 'genotype')
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Load SNP subject IDs from PLINK .fam ───────────────────────────────────
# .fam columns: FamilyID, SubjectID (IID), PaternalID, MaternalID, Sex, Phenotype
fam = pd.read_csv(FAM_FILE, sep=r'\s+', header=None,
                  names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'],
                  engine='python')
snp_ids = set(fam['IID'].str.strip())
print(f'SNP subjects in .fam:        {len(fam)}')

# ── 2. Load BIDS participants ──────────────────────────────────────────────────
pt = pd.read_csv(os.path.join(BIDS_DIR, 'participants.tsv'), sep='\t')
pt['adni_subject_id'] = pt['adni_subject_id'].str.strip()
print(f'BIDS subjects with 3T T1w:   {len(pt)}')

# ── 3. Find overlap ────────────────────────────────────────────────────────────
pt['has_snp'] = pt['adni_subject_id'].isin(snp_ids)
overlap = pt[pt['has_snp']].copy()

print(f'\nSubjects with BOTH MRI+SNP:  {len(overlap)}')
print(f'SNP-only (no MRI):           {len(snp_ids - set(pt["adni_subject_id"]))}')
print(f'MRI-only (no SNP):           {len(pt) - len(overlap)}')
print(f'\nDiagnosis breakdown:')
print(overlap['diagnosis_bl'].value_counts().to_string())
print(f'\nSex breakdown:')
print(overlap['sex'].value_counts().to_string())

# ── 4. Save outputs ────────────────────────────────────────────────────────────
out_cols = ['participant_id', 'adni_subject_id', 'diagnosis_bl', 'sex',
            'age', 'education_years', 'apoe4_dosage', 'apoe_genotype', 'site']
out_cols = [c for c in out_cols if c in overlap.columns]
overlap[out_cols].to_csv(os.path.join(OUT_DIR, 'subjects_with_snp_and_mri.tsv'),
                         sep='\t', index=False)

fam_overlap = fam[fam['IID'].isin(set(overlap['adni_subject_id']))].copy()
fam_overlap.to_csv(os.path.join(OUT_DIR, 'snp_subjects_with_mri.fam'),
                   sep=' ', index=False, header=False)

print(f'\nOutputs written to {OUT_DIR}/')
print('  subjects_with_snp_and_mri.tsv')
print('  snp_subjects_with_mri.fam')
