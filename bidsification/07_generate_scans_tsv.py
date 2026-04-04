"""
Script 07 — Generate scans.tsv Files
======================================
For each BIDS subject/session, generates a scans.tsv listing all
T1w NIfTI files and their acquisition dates.

BIDS scans.tsv format:
    filename                              acq_time
    anat/sub-002S0295_ses-20110602_T1w.nii.gz   2011-06-02T00:00:00

The scans.tsv is placed at:
    bids/sub-<label>/ses-<label>/sub-<label>_ses-<label>_scans.tsv
"""

import pandas as pd
import os
import glob

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT   = r"D:\ADNI_BIDS_project"
BIDS_DIR       = os.path.join(PROJECT_ROOT, "bids")
SESSION_MAP    = os.path.join(PROJECT_ROOT, "metadata", "session_map.csv")
SES_VISIT      = os.path.join(PROJECT_ROOT, "metadata", "ses_to_visit_code.csv")

# Load session date lookup (bids_ses = YYYYMMDD → ISO date for scans.tsv)
ses_visit = pd.read_csv(SES_VISIT, low_memory=False)
ses_date_lookup = {
    f"{row['bids_sub']}|{row['bids_ses']}": str(row["StudyDate"])[:10]
    for _, row in ses_visit.iterrows()
}

print(f"Loaded {len(ses_date_lookup):,} session date entries")

# ── Walk BIDS directory by session ────────────────────────────────────────────
session_dirs = glob.glob(os.path.join(BIDS_DIR, "sub-*", "ses-*"))
created = 0

for ses_dir in sorted(session_dirs):
    # Parse sub and ses labels
    parts = ses_dir.replace("\\", "/").split("/")
    sub_label = parts[-2].replace("sub-", "")
    ses_label = parts[-1].replace("ses-", "")
    
    anat_dir = os.path.join(ses_dir, "anat")
    if not os.path.isdir(anat_dir):
        continue
    
    nii_files = sorted(glob.glob(os.path.join(anat_dir, "*.nii.gz")))
    if not nii_files:
        continue
    
    # Get acquisition date
    key = f"{sub_label}|{ses_label}"
    study_date = ses_date_lookup.get(key, "")
    if study_date and study_date != "nan":
        # Format as ISO 8601 datetime (time unknown → T00:00:00)
        acq_time = f"{study_date}T00:00:00"
    else:
        acq_time = "n/a"
    
    # Build scans.tsv rows
    rows = []
    for nii_path in nii_files:
        rel_path = f"anat/{os.path.basename(nii_path)}"
        rows.append({"filename": rel_path, "acq_time": acq_time})
    
    scans_df = pd.DataFrame(rows)
    out_tsv = os.path.join(ses_dir, f"sub-{sub_label}_ses-{ses_label}_scans.tsv")
    scans_df.to_csv(out_tsv, sep="\t", index=False)
    created += 1

print(f"Created {created:,} scans.tsv files")
