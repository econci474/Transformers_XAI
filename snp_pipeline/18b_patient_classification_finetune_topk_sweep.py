"""
18b_patient_classification_finetune_topk_sweep.py
==================================================
Resumable hyperparameter-sweep driver for the BMFM-DNA-REF + LoRA patient
classifier (`18_patient_classification_finetune.py`), restricted to the
existing windows that contain the **Top-K SNPs by |GWAS z-score|**.

Design (keeps script 18 UNMODIFIED for the data path):
  * Rank SNPs that appear in ≥1 existing window by |z_score|, take the top-K,
    keep the already-built windows that contain any of them.
  * Write a FILTERED windows.tsv (same columns, kept rows, original order) and
    pass it to script 18 as `--windows`. script 18's `build_chrom_map()`
    already yields contiguous chrom indices from whatever windows.tsv it is
    given, so NO change to script 18 is needed — the dataset only tokenises
    the kept windows ⇒ ≈ (K/183)× the cost (~7× faster at K=25).
  * Per HP point, pass a UNIQUE `--output-root <base>/<hp_tag>` and
    `--wandb-run-name topk<K>__<hp_tag>` so out_dir / resumable checkpoint /
    W&B runs never collide; skip a point whose final metrics.json exists
    (resumable across Colab sessions).

The cosine-LR schedule is the one (default-guarded) addition inside script 18
(`--lr-scheduler cosine --warmup-frac …`); everything else is its existing CLI.

Usage (Colab)
-------------
    python snp_pipeline/18b_patient_classification_finetune_topk_sweep.py \\
        --base-root /content/drive/MyDrive/ADNI_SNP/classification_finetune_topk \\
        --tier 1 --seeds 0
    # add --dry-run to print the selection + commands without running
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT18 = Path(__file__).parent / "18_patient_classification_finetune.py"
DRIVE = "/content/drive/MyDrive/ADNI_SNP"

# Fixed = the 183-window anchor config (for a like-for-like comparison).
FIXED = dict(label_mode="ever_convert", aggregation="chrom_mean_then_attn",
             per_window_pool="mean", loss="weighted_bce", precision="bf16",
             micro_batch=16, lora_dropout=0.1, weight_decay=1e-3,
             warmup_frac=0.02)
# Tier-1 = LR sweep at the anchor rank; Tier-2 adds a rank axis at anchor LR.
TIER1 = [dict(lr=3e-4, lora_r=16, lora_alpha=32),
         dict(lr=1e-4, lora_r=16, lora_alpha=32),   # = anchor cfg @ Top-K
         dict(lr=3e-5, lora_r=16, lora_alpha=32)]
TIER2_EXTRA = [dict(lr=1e-4, lora_r=8,  lora_alpha=16),
               dict(lr=1e-4, lora_r=32, lora_alpha=64)]


def _fs(x) -> str:
    return str(x).replace("-", "m").replace(".", "p").replace("+", "")


def hp_tag(p: dict, k: int, seed: int) -> str:
    return (f"k{k}_lr{_fs(p['lr'])}_r{p['lora_r']}_a{p['lora_alpha']}"
            f"_wd{_fs(FIXED['weight_decay'])}_s{seed}")


def select_windows_for_top_snps(wins_df: pd.DataFrame, gwas_labels_path: Path,
                                 k: int) -> tuple[list[str], pd.DataFrame]:
    """Return (top_snps, kept_rows). Rank SNPs that appear in ≥1 window by
    |z_score| (deterministic tie-break by SNP id), take top-k, keep the
    windows containing any of them. No re-windowing."""
    g = pd.read_csv(gwas_labels_path, sep="\t", low_memory=False)
    zmap = {str(s): abs(float(z))
            for s, z in zip(g["SNP"], g["z_score"]) if pd.notna(z)}
    snp_to_windows: dict[str, list[str]] = {}
    for r in wins_df.itertuples(index=False):
        for s in str(r.snp_ids).split("|"):
            snp_to_windows.setdefault(s, []).append(str(r.window_id))
    ranked = sorted((s for s in snp_to_windows if s in zmap),
                    key=lambda s: (-zmap[s], s))
    top_snps = ranked[:k]
    keep: set[str] = set()
    for s in top_snps:
        keep.update(snp_to_windows[s])
    kept = wins_df[wins_df["window_id"].astype(str).isin(keep)].copy()
    # provenance: per-kept-window max|z| over its SNPs
    kept["max_abs_z"] = kept["snp_ids"].map(
        lambda s: max((zmap.get(x, 0.0) for x in str(s).split("|")), default=0.0))
    return top_snps, kept


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-root", type=Path, required=True,
                    help="Sweep root; each HP point goes in <base-root>/<hp_tag>/.")
    ap.add_argument("--sequences", default=f"{DRIVE}/inputs/patient_sequences/all_patients.csv")
    ap.add_argument("--windows",   default=f"{DRIVE}/inputs/windows/windows.tsv",
                    help="FULL windows.tsv; the driver writes a Top-K-filtered copy.")
    ap.add_argument("--splits-root", default=f"{DRIVE}/inputs/no_cdr_stratified_ever_convert/tabular/baseline")
    ap.add_argument("--gwas-labels", default=f"{DRIVE}/inputs/external_gwas_labels.tsv",
                    help="external_gwas_labels.tsv (cols SNP, z_score).")
    ap.add_argument("--top-k", type=int, default=25)
    ap.add_argument("--tier", type=int, choices=[1, 2], default=1)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0])
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--wandb-project", default="bmfm_classification_lora_ft_topk")
    ap.add_argument("--device", default=None, help="pass-through to script 18")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    args.base_root.mkdir(parents=True, exist_ok=True)

    # ── Build the Top-K-SNP filtered windows.tsv (once, shared by all points) ──
    wins_df = pd.read_csv(args.windows, sep="\t")
    top_snps, kept = select_windows_for_top_snps(
        wins_df, Path(args.gwas_labels), args.top_k)
    kept_in_order = wins_df[wins_df["window_id"].astype(str).isin(
        set(kept["window_id"].astype(str)))].copy()           # original order
    chroms = sorted(kept_in_order["chrom"].astype(str).unique(),
                    key=lambda c: int(c) if c.isdigit() else 99)
    filt_windows = args.base_root / f"topk{args.top_k}_windows.tsv"
    prov = args.base_root / f"topk{args.top_k}_selected.tsv"
    print(f"[topk] {len(top_snps)} top SNPs → {len(kept_in_order)}/{len(wins_df)} "
          f"windows over chroms {chroms}")
    print(f"[topk] top SNPs (by |z|): {', '.join(top_snps[:10])}"
          f"{' …' if len(top_snps) > 10 else ''}")
    if args.dry_run:
        print(f"[dry-run] would write {filt_windows} and {prov}")
    else:
        kept_in_order[wins_df.columns].to_csv(filt_windows, sep="\t", index=False)
        kept.sort_values("max_abs_z", ascending=False)[
            ["window_id", "chrom", "max_abs_z", "snp_ids"]
        ].to_csv(prov, sep="\t", index=False)
        print(f"[topk] wrote {filt_windows}\n[topk] wrote {prov}")

    points = TIER1 + (TIER2_EXTRA if args.tier == 2 else [])
    grid = [(p, s) for p in points for s in args.seeds]
    print(f"\n[sweep] tier {args.tier}: {len(grid)} run(s) "
          f"(K={args.top_k}, epochs≤{args.epochs}, patience {args.patience}, "
          f"cosine LR)\n")

    done = skipped = failed = 0
    for i, (p, seed) in enumerate(grid, 1):
        tag = hp_tag(p, args.top_k, seed)
        out_root = args.base_root / tag
        metrics = (out_root / FIXED["label_mode"] / "raw_bmfm_ref" / "lora"
                   / FIXED["aggregation"] / f"seed_{seed}" / "metrics.json")
        if metrics.exists():
            skipped += 1
            print(f"[{i}/{len(grid)}] {tag}  SKIP (metrics.json exists)")
            continue
        cmd = [sys.executable, "-u", str(SCRIPT18),
               "--sequences", str(args.sequences),
               "--windows", str(filt_windows),
               "--splits-root", str(args.splits_root),
               "--label-mode", FIXED["label_mode"],
               "--aggregation", FIXED["aggregation"],
               "--per-window-pool", FIXED["per_window_pool"],
               "--loss", FIXED["loss"], "--precision", FIXED["precision"],
               "--micro-batch", str(FIXED["micro_batch"]),
               "--grad-checkpoint",
               "--lora-dropout", str(FIXED["lora_dropout"]),
               "--weight-decay", str(FIXED["weight_decay"]),
               "--epochs", str(args.epochs), "--patience", str(args.patience),
               "--lr-scheduler", "cosine",
               "--warmup-frac", str(FIXED["warmup_frac"]),
               "--lr", str(p["lr"]),
               "--lora-r", str(p["lora_r"]),
               "--lora-alpha", str(p["lora_alpha"]),
               "--seed", str(seed),
               "--output-root", str(out_root),
               "--wandb", "--wandb-project", args.wandb_project,
               "--wandb-run-name", f"topk{args.top_k}__{tag}",
               "--log-every", "25"]
        if args.device:
            cmd += ["--device", args.device]
        if args.dry_run:
            print(f"[{i}/{len(grid)}] {tag}  DRY: {' '.join(cmd)}")
            continue
        print(f"[{i}/{len(grid)}] {tag}  RUN")
        r = subprocess.run(cmd)                       # stream live (no capture)
        if r.returncode == 0 and metrics.exists():
            done += 1
            print(f"[{i}/{len(grid)}] {tag}  OK")
        else:
            failed += 1
            print(f"[{i}/{len(grid)}] {tag}  FAIL rc={r.returncode}")

    print(f"\n[sweep] done={done} skipped={skipped} failed={failed} "
          f"→ {args.base_root}")


if __name__ == "__main__":
    main()
