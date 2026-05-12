import pandas as pd
df = pd.read_csv(r"D:\ADNI_BIDS_project\metadata\ses_to_visit_code.csv")
print("Columns:", df.columns.tolist())
print(df.head(15).to_string())
print(f"\nRows: {len(df)}")

# Also check clinical VISCODE_2
clin = pd.read_csv(r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\longitudinal\master_clinical_tabular.csv",
                   usecols=["Patient_ID", "VISCODE_long", "VISCODE_2"], low_memory=False)
print(f"\nClinical master cols: VISCODE_long unique = {clin['VISCODE_long'].nunique()}, VISCODE_2 unique = {clin['VISCODE_2'].nunique()}")
print("VISCODE_2 samples:", sorted(clin["VISCODE_2"].dropna().unique())[:20])
print("VISCODE_long samples:", sorted(clin["VISCODE_long"].dropna().unique(), key=lambda s: int(s[1:]) if s.startswith("m") and s[1:].isdigit() else -1)[:20])

# Check MRI master bids_ses
mri = pd.read_csv(r"D:\ADNI_BIDS_project\derivatives\mri_clinical_matched\master_mri_clinical_matched_labels.csv",
                  usecols=["bids_sub", "bids_ses"], low_memory=False)
print(f"\nMRI master bids_ses unique: {sorted(mri['bids_ses'].unique())}")
