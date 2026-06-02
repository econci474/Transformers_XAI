# env_logs/ — dated snapshots of each conda environment

**Why this exists:** a silent `torch 2.3.1 → 2.4.1` drift in the `clinical` env (the live
env diverged from `clinical_pipeline/environment.yml`, which pins 2.3.1) broke ModernBERT
**full fine-tune** with step-1 `grad_norm=NaN` and cost real compute to diagnose. (The same
env had also drifted python 3.11 → 3.10.) These dated, version-controlled snapshots make env
drift visible and reversible.

## Layout

```
env_logs/<env>/<env>_<YYYY-MM-DD>_<host>.txt
```

| env | used for | spec file |
|---|---|---|
| `clinical` | ModernBERT / BioClinical-ModernBERT fine-tuning (clinical_pipeline) | `clinical_pipeline/environment.yml` |
| `snp` | PLINK QC + GWAS data prep (snp_pipeline 01–08e) | — |
| `bmfm` | BMFM-DNA-SNP fine-tuning (snp_pipeline 09, CSD3) | — |
| `survml` | survival analysis (clinical_pipeline 02k / 02l) | — |

Each file records: python version, a short key-library version block, the full `pip freeze`,
and `conda list`.

## How to snapshot

Run on the machine whose env you want to record (envs differ across CSD3 / dept-L4 / local —
the `<host>` suffix keeps them distinct):

```bash
conda activate <env>
bash snapshot_env.sh            # infers env from the active one
# or:  bash snapshot_env.sh <env>
git add env_logs/ && git commit -m "env(<env>): snapshot $(date +%F)" && git push
```

Take a fresh snapshot **whenever you (re)build or `pip install` into an env**, and **always
after confirming an env is in a known-good state**, so there is a point to roll back to.

## Known-good references

- **clinical** (CSD3, 2026-06-02): `torch==2.3.1+cu121`, `torchvision==0.18.1`,
  `transformers 4.57.6`, `tokenizers 0.22.2`, `huggingface_hub 0.36.2`, `numpy 1.26.4`,
  python 3.10. **`torch 2.4.x` NaNs ModernBERT full fine-tune — do NOT let it drift up.**
  `flash-attn` is absent on the CSD3 clinical env (SDPA is used); that is fine.
