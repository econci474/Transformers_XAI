"""
Quick diagnostic: how many subjects have sc-only (no bl row)?
Do any subjects have NaN Label_bl_multi?

Run:  python clinical_pipeline/_check_bl_sc_unknown.py
"""
import pandas as pd

MASTER = r"D:\ADNI_BIDS_project\derivatives\clinical\no_cdr_stratified\tabular\longitudinal\master_clinical_tabular.csv"

df = pd.read_csv(MASTER, usecols=["Patient_ID", "VISCODE_long", "Label_bl_multi"],
                 low_memory=False)

all_pids = df["Patient_ID"].unique()
bl_pids  = set(df[df["VISCODE_long"] == "bl"]["Patient_ID"].unique())
sc_pids  = set(df[df["VISCODE_long"] == "sc"]["Patient_ID"].unique())

sc_only = sc_pids - bl_pids  # have sc but no bl

print(f"Total subjects          : {len(all_pids)}")
print(f"Subjects with bl row    : {len(bl_pids)}")
print(f"Subjects with sc row    : {len(sc_pids)}")
print(f"Subjects with sc but NO bl: {len(sc_only)}")
if sc_only:
    print(f"  IDs: {sorted(sc_only)}")
else:
    print("  (none)")

# Check for NaN Label_bl_multi among baseline-like rows
bl_rows = df[df["VISCODE_long"] == "bl"]
sc_rows = df[df["VISCODE_long"] == "sc"]
bl_set  = set(bl_rows["Patient_ID"])
baseline = pd.concat([bl_rows, sc_rows[~sc_rows["Patient_ID"].isin(bl_set)]]).drop_duplicates("Patient_ID")

nan_dx = baseline[baseline["Label_bl_multi"].isna()]
print(f"\nBaseline-like rows total: {len(baseline)}")
print(f"Rows with NaN Label_bl_multi: {len(nan_dx)}")
if len(nan_dx) > 0:
    print(f"  IDs: {sorted(nan_dx['Patient_ID'].tolist())}")
else:
    print("  (none — all subjects have a valid diagnosis)")
