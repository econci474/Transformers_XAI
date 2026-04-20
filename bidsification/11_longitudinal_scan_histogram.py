"""
Script 11 — Longitudinal MRI Visit Table and Histogram
========================================================
For subjects with both SNP and 3T T1w MRI data, produces:

  bids/imaging/mri_visit_dates_snp_subjects.csv
      Pivot table: rows = subjects, columns = BIDS visit codes (bl, m12, …),
      values = acquisition date (YYYY-MM-DD) or 'n/a'.

  bids/imaging/mri_scan_count_histogram.png
      Stacked bar chart of scan counts per subject, colour-coded by
      baseline diagnosis.

Requires: bids/genotype/subjects_with_snp_and_mri.tsv (from script 10)
"""

import pandas as pd
import os
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
BIDS_DIR  = r'D:\ADNI_BIDS_project\bids'
GENO_DIR  = os.path.join(BIDS_DIR, 'genotype')
OUT_DIR   = os.path.join(BIDS_DIR, 'imaging')
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Load overlapping subjects ──────────────────────────────────────────────
overlap  = pd.read_csv(os.path.join(GENO_DIR, 'subjects_with_snp_and_mri.tsv'), sep='\t')
snp_pids = set(overlap['participant_id'])
print(f'Overlapping subjects: {len(snp_pids)}')

# ── 2. Collect session dates from scans.tsv files ─────────────────────────────
records = []
for st_path in glob.glob(os.path.join(BIDS_DIR, 'sub-*', 'ses-*', 'scans.tsv')):
    norm   = st_path.replace('\\', '/')
    tokens = norm.split('/')
    sub = next((t for t in tokens if t.startswith('sub-')), None)
    ses = next((t for t in tokens if t.startswith('ses-')), None)
    if sub not in snp_pids or sub is None or ses is None:
        continue
    df  = pd.read_csv(st_path, sep='\t')
    acq = df['acq_time'].iloc[0] if 'acq_time' in df.columns else 'n/a'
    date = str(acq)[:10] if acq not in ('n/a', None, 'nan') else 'n/a'
    records.append({'participant_id': sub, 'ses': ses.replace('ses-', ''), 'date': date})

long = pd.DataFrame(records)
print(f'Session records: {len(long)} across {long["participant_id"].nunique()} subjects')

# ── 3. Pivot to subject × visit-code ─────────────────────────────────────────
visit_order = ['sc', 'bl', 'm06', 'm12', 'm18', 'm24', 'm36',
               'm48', 'm60', 'm72', 'm84', 'm96', 'm108', 'm120']
all_visits  = sorted(long['ses'].unique(),
                     key=lambda v: visit_order.index(v) if v in visit_order else 999)

pivot = (long.pivot_table(index='participant_id', columns='ses', values='date', aggfunc='first')
         .reindex(columns=all_visits)
         .fillna('n/a'))
pivot.index.name = 'participant_id'
pivot = pivot.sort_index()

csv_path = os.path.join(OUT_DIR, 'mri_visit_dates_snp_subjects.csv')
pivot.to_csv(csv_path)
print(f'Saved: {csv_path}  ({pivot.shape[0]} rows × {pivot.shape[1]} visits)')

# ── 4. Scan counts per subject ────────────────────────────────────────────────
scan_counts = (pd.DataFrame(records)
               .groupby('participant_id')
               .size()
               .reset_index(name='n_scans'))
df = scan_counts.merge(overlap[['participant_id', 'diagnosis_bl']], on='participant_id', how='left')
df['diagnosis_bl'] = df['diagnosis_bl'].fillna('Unknown')

print(f'\nScan count distribution:\n{scan_counts["n_scans"].value_counts().sort_index().to_string()}')
print(f'Mean scans: {scan_counts["n_scans"].mean():.1f}  '
      f'(range {scan_counts["n_scans"].min()}–{scan_counts["n_scans"].max()})')

# ── 5. Stacked histogram ──────────────────────────────────────────────────────
diag_order  = ['CognitivelyNormal', 'EarlyMCI', 'LateMCI', 'AlzheimersDisease']
diag_labels = {'CognitivelyNormal': 'CN',
               'EarlyMCI':          'Early MCI',
               'LateMCI':           'Late MCI',
               'AlzheimersDisease': 'AD'}
palette     = {'CognitivelyNormal': '#4DAF4A',
               'EarlyMCI':          '#377EB8',
               'LateMCI':           '#FF7F00',
               'AlzheimersDisease': '#E41A1C'}

n_scans_range = sorted(df['n_scans'].unique())
stacked = pd.DataFrame(0, index=n_scans_range, columns=diag_order)
for (n, dx), grp in df.groupby(['n_scans', 'diagnosis_bl']):
    if dx in diag_order:
        stacked.loc[n, dx] = len(grp)

fig, ax = plt.subplots(figsize=(10, 6))
bottoms = np.zeros(len(n_scans_range))
diag_counts = df['diagnosis_bl'].value_counts()
for dx in diag_order:
    vals = stacked[dx].values.astype(int)
    n = int(diag_counts.get(dx, 0))
    ax.bar(n_scans_range, vals, bottom=bottoms,
           color=palette[dx], edgecolor='white', linewidth=0.6,
           width=0.7, label=f"{diag_labels[dx]}  (n={n})")
    bottoms += vals

totals = stacked.sum(axis=1).values.astype(int)
for x, total in zip(n_scans_range, totals):
    ax.text(x, total + 0.5, str(total), ha='center', va='bottom',
            fontsize=10, fontweight='bold', color='#222222')

ax.set_xlabel('Number of 3T T1w MRI Scans', fontsize=13)
ax.set_ylabel('Number of Subjects', fontsize=13)
ax.set_title('Histogram of Longitudinal 3T T1w MRI Scans for Subjects with SNP Data',
             fontsize=13, fontweight='bold', pad=14)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.tick_params(axis='both', labelsize=11)
ax.set_xlim(0.5, max(n_scans_range) + 0.5)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(False)
ax.legend(title='Baseline Diagnosis', title_fontsize=11, fontsize=10,
          loc='upper right', framealpha=0.9)

plt.tight_layout()
hist_path = os.path.join(OUT_DIR, 'mri_scan_count_baseline_histogram.png')
plt.savefig(hist_path, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {hist_path}')
