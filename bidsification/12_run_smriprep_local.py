"""
Script 12 — Run sMRIprep Locally via Docker (Parallel)
=======================================================
Iterates over all subjects in subjects_with_snp_and_mri.tsv and runs
sMRIprep in Docker with N_PARALLEL simultaneous containers.

Output spaces (matching the HPC apptainer command):
  --output-spaces T1w MNI152NLin2009cAsym:res-1

Recommended settings for 8-CPU machine:
  N_PARALLEL   = 2   (2 subjects at once)
  NPROCS       = 4   (4 CPUs per subject)
  OMP_NTHREADS = 2
  MEM_GB       = 8   (adjust to your available RAM / N_PARALLEL)

Time estimate (printed at start):
  ~1598 sessions x ~45 min / N_PARALLEL jobs ~ days (see output)

Usage:
  python 12_run_smriprep_local.py
  python 12_run_smriprep_local.py --dry-run    # preview commands only
"""

import subprocess, os, sys, time, datetime, argparse, threading
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── Configuration ──────────────────────────────────────────────────────────────
BIDS_DIR      = r'D:\ADNI_BIDS_project\bids'
OUT_DIR       = r'D:\ADNI_BIDS_project\derivatives\smriprep'
WORK_DIR      = r'D:\ADNI_BIDS_project\work'
FS_LICENSE    = r'D:\FreeSurfer_licence\license.txt'
SNP_TSV       = r'D:\ADNI_BIDS_project\bids\genotype\subjects_with_snp_and_mri.tsv'
LOG_FILE      = r'D:\ADNI_BIDS_project\smriprep_run_log.tsv'
SMRIPREP_IMG  = 'nipreps/smriprep:0.19.2'
# Full path to Docker on Windows (not always in PATH)
DOCKER_CMD    = r'C:\Program Files\Docker\Docker\resources\bin\docker.exe'

# ── Resource allocation per Docker container ───────────────────────────────────
# 8-CPU machine: 2 parallel x 4 CPUs = 8 total, 2 x 8 GB = 16 GB RAM
N_PARALLEL    = 2   # simultaneous Docker containers
NPROCS        = 4   # CPUs per container  (N_PARALLEL x NPROCS <= total CPUs)
OMP_NTHREADS  = 2   # threads per process (<= NPROCS)
MEM_GB        = 8   # RAM per container   (N_PARALLEL x MEM_GB <= total RAM)


OUTPUT_SPACES = 'T1w MNI152NLin2009cAsym:res-1'

# ── Parse args ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--dry-run', action='store_true', help='Print commands without executing')
args = parser.parse_args()
DRY_RUN = args.dry_run

# ── Load subjects ──────────────────────────────────────────────────────────────
tsv = pd.read_csv(SNP_TSV, sep='\t')
subjects = tsv['participant_id'].str.replace('sub-', '', regex=False).tolist()

# ── Time estimate ──────────────────────────────────────────────────────────────
N_SESSIONS    = 1598
MINS_PER_SES  = 45
wall_h = N_SESSIONS * MINS_PER_SES / 60 / N_PARALLEL
print('=' * 65)
print(f'sMRIprep Docker run -- {len(subjects)} subjects, ~{N_SESSIONS} sessions')
print(f'Parallelism : {N_PARALLEL} subjects x {NPROCS} CPUs  |  {MEM_GB} GB RAM each')
print(f'Estimate    : {N_SESSIONS} x {MINS_PER_SES} min / {N_PARALLEL} = '
      f'{wall_h:.0f} h  ({wall_h/24:.1f} days wall-clock)')
print()
print('WARNING: Even with parallelism this is impractical for all 1,598 sessions.')
print('  For the full dataset, use Cambridge HPC (Singularity, SLURM array jobs).')
print('  Consider running baseline-only sessions locally to test pipeline.')
print('=' * 65)

os.makedirs(OUT_DIR,  exist_ok=True)
os.makedirs(WORK_DIR, exist_ok=True)

# ── Logging ────────────────────────────────────────────────────────────────────
log_lock = threading.Lock()
if not os.path.exists(LOG_FILE):
    with open(LOG_FILE, 'w') as f:
        f.write('participant_id\tstart_time\tend_time\tduration_min\treturn_code\n')

def log_result(pid, t_start, t_end, rc):
    dur = round((t_end - t_start) / 60, 1)
    with log_lock:
        with open(LOG_FILE, 'a') as f:
            f.write(f'{pid}\t'
                    f'{datetime.datetime.fromtimestamp(t_start).isoformat()}\t'
                    f'{datetime.datetime.fromtimestamp(t_end).isoformat()}\t'
                    f'{dur}\t{rc}\n')

# ── Docker path helpers ────────────────────────────────────────────────────────
def win_to_docker(path):
    """Convert Windows path D:\\foo\\bar -> /d/foo/bar for Docker on Windows."""
    p = path.replace('\\', '/')
    if len(p) >= 2 and p[1] == ':':
        p = '/' + p[0].lower() + p[2:]
    return p

BIDS_D = win_to_docker(BIDS_DIR)
OUT_D  = win_to_docker(OUT_DIR)
LIC_D  = win_to_docker(FS_LICENSE)

# ── Per-subject function ───────────────────────────────────────────────────────
def run_subject(pid, index, total):
    sub_out = os.path.join(OUT_DIR, 'smriprep', f'sub-{pid}')
    if os.path.isdir(sub_out):
        print(f'[{index}/{total}] sub-{pid}: already done, skipping.')
        return pid, 'skipped', 0

    # Each subject gets its own scratch directory to avoid conflicts
    work_sub = win_to_docker(os.path.join(WORK_DIR, pid))

    cmd = [
        DOCKER_CMD, 'run', '--rm',
        '-v', f'{BIDS_D}:/data:ro',
        '-v', f'{OUT_D}:/out',
        '-v', f'{work_sub}:/scratch',
        '-v', f'{LIC_D}:/opt/freesurfer/license.txt:ro',
        SMRIPREP_IMG,
        '/data', '/out', 'participant',
        '--participant-label', pid,
        '--nprocs',          str(NPROCS),
        '--omp-nthreads',    str(OMP_NTHREADS),
        '--mem-gb',          str(MEM_GB),
        '--output-spaces',   'T1w', 'MNI152NLin2009cAsym:res-1',
        '--fs-no-reconall',
        '-w', '/scratch',
        '--fs-license-file', '/opt/freesurfer/license.txt',
    ]

    print(f'[{index}/{total}] sub-{pid} START -- {datetime.datetime.now():%H:%M:%S}')
    if DRY_RUN:
        print('  CMD: ' + ' '.join(cmd))
        print('  DRY RUN -- skipping.')
        return pid, 'dry-run', 0

    os.makedirs(os.path.join(WORK_DIR, pid), exist_ok=True)

    t0 = time.time()
    result = subprocess.run(cmd, capture_output=False)
    t1 = time.time()
    rc = result.returncode
    log_result(f'sub-{pid}', t0, t1, rc)

    status = 'ok' if rc == 0 else f'FAILED({rc})'
    print(f'[{index}/{total}] sub-{pid} {status} -- {(t1-t0)/60:.1f} min')
    return pid, status, rc

# ── Main parallel loop ─────────────────────────────────────────────────────────
results = {'ok': [], 'skipped': [], 'failed': []}

with ThreadPoolExecutor(max_workers=N_PARALLEL) as pool:
    futures = {pool.submit(run_subject, pid, i+1, len(subjects)): pid
               for i, pid in enumerate(subjects)}
    for future in as_completed(futures):
        pid, status, rc = future.result()
        if status == 'skipped':
            results['skipped'].append(pid)
        elif rc == 0:
            results['ok'].append(pid)
        else:
            results['failed'].append(pid)

# ── Summary ────────────────────────────────────────────────────────────────────
print('\n' + '='*65)
print(f'Run complete.')
print(f'  Skipped (already done): {len(results["skipped"])}')
print(f'  Newly processed:        {len(results["ok"])}')
print(f'  Failed:                 {len(results["failed"])}')
if results['failed']:
    print(f'  Failed IDs: {results["failed"]}')
print(f'  Log: {LOG_FILE}')
