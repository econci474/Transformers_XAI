"""
08_hp_optuna_full_ft.py — Optuna Bayesian full_ft HP search (BrainMVP / ViT-MAE)
=================================================================================
BrainMVP + ViT-MAE are strongest FULL fine-tuned (not frozen), and full_ft is
expensive (A100, ~hours/run), so a Bayesian (TPE) search beats a grid. A
wandb-native sweep would need the controller to reach the wandb server, but
CSD3 compute nodes are offline — so we use **Optuna** (local Bayesian, no
server) with a file-based **JournalStorage** study shared by parallel array
agents on Lustre, and log each trial to **wandb offline** (sync later).

Selection (user decision): SEARCH at **seed 0**, optimise **val balanced
accuracy** (`config.best_val_balanced_acc`); then REFIT the **top-k** configs on
**3 seeds** and pick the winner by **mean val-bACC**. macro-AUC is recorded for
the paper comparison but not optimised.

It subprocess-calls the EXISTING trainers (no edits to them); each search trial
writes to a unique trial dir, we read its metrics.json, then DELETE its *.pt
(checkpoint-free search → "metrics_only" without touching the trainers).

Modes:
  search  : run --n-trials trials against the shared study (array → 1 trial/agent).
  report  : print the study leaderboard + write <study-dir>/top_configs.json (top-k).
  refit   : --refit-index decodes (config_k, seed) from top_configs.json → one
            full_ft run KEEPING the checkpoint, under --winner-dir.
  winner  : aggregate the refit metrics.json → mean±std per config, pick the
            winner by mean val-bACC, write winner_summary.csv.

Examples (CSD3):
  # search agent (one trial) — array driver passes this:
  python 08_hp_optuna_full_ft.py --mode search --arch brainmvp --n-trials 1 \
      --study-dir .../optuna/brainmvp_T2
  python 08_hp_optuna_full_ft.py --mode report --arch brainmvp --study-dir .../optuna/brainmvp_T2 --topk 3
  python 08_hp_optuna_full_ft.py --mode refit  --arch brainmvp --study-dir .../optuna/brainmvp_T2 \
      --refit-index $SLURM_ARRAY_TASK_ID --winner-dir .../brainmvp_T2_winner
  python 08_hp_optuna_full_ft.py --mode winner --arch brainmvp --winner-dir .../brainmvp_T2_winner
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

MRI = Path(__file__).resolve().parent

# ── CSD3 post-exclusion paths (same master/splits as the cached-head + fusion) ──
DEF_MASTER = "/home/ec474/rds/hpc-work/ADNI_MRI/master_mri_clinical_matched_viscode2_extended_post_exclusion.csv"
DEF_DATA   = "/home/ec474/rds/hpc-work/ADNI_CL/no_cdr_stratified_post_exclusion/tabular/baseline"
VIT_PRE    = "/home/ec474/rds/hpc-work/ViT_pretrained/ViT_B_pretrained_noaug_mae75_BRATS2023_IXI_OASIS3_seed_8456_999_077000.pth.tar"
VIT_IN     = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/vit_inputs"
BMVP_PRE   = "/home/ec474/rds/hpc-work/ViT_pretrained/BrainMVP_uniformer.pt"
BMVP_IN    = "/home/ec474/rds/hpc-work/ADNI_SMRIPREP/derivatives/brainmvp_inputs"


def vit_space(trial):
    return {"--lr":             trial.suggest_float("lr", 1e-5, 1e-3, log=True),
            "--weight_decay":   trial.suggest_float("weight_decay", 1e-6, 1e-2, log=True),
            "--drop_path_rate": trial.suggest_float("drop_path_rate", 0.0, 0.3),
            "--attn_dropout":   trial.suggest_float("attn_dropout", 0.0, 0.3),
            "--label_smoothing": trial.suggest_categorical("label_smoothing", [0.0, 0.1]),
            "--augment":        trial.suggest_categorical("augment", ["none", "random", "plus_original"])}


def bmvp_space(trial):
    return {"--lr":             trial.suggest_float("lr", 1e-5, 1e-3, log=True),
            "--weight_decay":   trial.suggest_float("weight_decay", 1e-6, 1e-2, log=True),
            "--drop_rate":      trial.suggest_float("drop_rate", 0.0, 0.3),
            "--label_smoothing": trial.suggest_categorical("label_smoothing", [0.0, 0.1]),
            "--augment":        trial.suggest_categorical("augment", ["none", "stochastic", "plus_original"])}


ARCHS = {
    "vit_mae":  dict(script=MRI / "04_supervised_finetuning_ViT.py", model_slug="ViT_B_mae75",
                     pretrained=VIT_PRE, inputs_flag="--vit_inputs_dir", inputs=VIT_IN,
                     space=vit_space, wandb_project="vit_mae_T2_optuna"),
    "brainmvp": dict(script=MRI / "brain_mvp" / "04_supervised_finetuning_BrainMVP.py",
                     model_slug="BrainMVP_uniformer", pretrained=BMVP_PRE,
                     inputs_flag="--brainmvp_inputs_dir", inputs=BMVP_IN,
                     space=bmvp_space, wandb_project="brainmvp_T2_optuna"),
}


def make_storage(journal_path: str):
    """File-based JournalStorage (NFS/Lustre-safe for parallel agents)."""
    import optuna
    try:                                            # optuna >= 3.4 layout
        from optuna.storages.journal import JournalFileBackend
        return optuna.storages.JournalStorage(JournalFileBackend(journal_path))
    except Exception:                               # older layout
        from optuna.storages import JournalFileStorage
        return optuna.storages.JournalStorage(JournalFileStorage(journal_path))


def _read_metrics(run_dir: Path):
    """Find the trainer's metrics.json under run_dir; return (val_bacc, test_bacc,
    test_auc) or None if absent."""
    hits = glob.glob(str(run_dir / "**" / "metrics.json"), recursive=True)
    if not hits:
        return None
    d = json.load(open(hits[0]))
    cfg, tm = d.get("config", {}), d.get("test_metrics", {})
    auc = next((tm[k] for k in ("auc_roc_ovr", "auc_roc") if k in tm), None)
    return (cfg.get("best_val_balanced_acc"), tm.get("balanced_acc"), auc)


def _run_trainer(arch, a, hp_flags, seed, out_dir, keep_ckpt, extra_env):
    acfg = ARCHS[arch]
    cmd = [sys.executable, str(acfg["script"]),
           "--task", a.task, "--seed", str(seed), "--strategy", "full_ft",
           "--long", a.long,
           "--matched_labels_csv", a.matched_labels_csv, "--data_dir", a.data_dir,
           acfg["inputs_flag"], a.inputs_dir or acfg["inputs"],
           "--pretrained_ckpt", a.pretrained or acfg["pretrained"],
           "--out_dir", str(out_dir), "--num_workers", str(a.num_workers)]
    if a.search_epochs and not keep_ckpt:
        cmd += ["--epochs", str(a.search_epochs)]
    for k, v in hp_flags.items():
        cmd += [k, str(v)]
    if a.wandb:
        cmd += ["--wandb", "--wandb_project", acfg["wandb_project"]]
    env = dict(os.environ, **extra_env)
    print("  +", " ".join(cmd), flush=True)
    return subprocess.run(cmd, env=env).returncode


def mode_search(a):
    import optuna
    acfg = ARCHS[a.arch]
    study_dir = Path(a.study_dir); (study_dir / "trials").mkdir(parents=True, exist_ok=True)
    storage = make_storage(str(study_dir / "optuna_journal.log"))
    study = optuna.create_study(direction="maximize", study_name=f"{a.arch}_{a.task}",
                                storage=storage, load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=0, n_startup_trials=10))

    def objective(trial):
        hp = acfg["space"](trial)
        run_dir = study_dir / "trials" / f"trial_{trial.number}"
        rc = _run_trainer(a.arch, a, hp, seed=0, out_dir=run_dir, keep_ckpt=False,
                          extra_env={"WANDB_MODE": "offline", "WANDB_DIR": str(run_dir)})
        if rc != 0:
            raise optuna.TrialPruned()             # crash/OOM -> don't poison TPE model
        m = _read_metrics(run_dir)
        if m is None or m[0] is None:
            raise optuna.TrialPruned()
        val_bacc, test_bacc, test_auc = m
        trial.set_user_attr("test_bacc", test_bacc)
        trial.set_user_attr("test_auc", test_auc)
        trial.set_user_attr("hp", {k.lstrip("-"): v for k, v in hp.items()})
        # checkpoint-free search: drop the trial's *.pt (keep metrics.json).
        for p in glob.glob(str(run_dir / "**" / "*.pt"), recursive=True):
            try: os.remove(p)
            except OSError: pass
        print(f"  trial {trial.number}: val_bACC={val_bacc:.4f} test_bACC={test_bacc} "
              f"test_AUC={test_auc}", flush=True)
        return float(val_bacc)

    study.optimize(objective, n_trials=a.n_trials)
    print(f"[search] {a.arch} study now has {len(study.trials)} trials; "
          f"best val_bACC={study.best_value:.4f}")


def mode_report(a):
    import optuna
    storage = make_storage(str(Path(a.study_dir) / "optuna_journal.log"))
    study = optuna.load_study(study_name=f"{a.arch}_{a.task}", storage=storage)
    done = [t for t in study.trials if t.value is not None]
    done.sort(key=lambda t: t.value, reverse=True)
    print(f"\n[report] {a.arch} {a.task} — {len(done)} completed trials (top {a.topk}):")
    top = []
    for rank, t in enumerate(done[:a.topk]):
        hp = t.user_attrs.get("hp", t.params)
        print(f"  #{rank} val_bACC={t.value:.4f}  test_bACC={t.user_attrs.get('test_bacc')}  "
              f"test_AUC={t.user_attrs.get('test_auc')}  {hp}")
        top.append({"rank": rank, "trial": t.number, "val_bacc": t.value,
                    "test_bacc": t.user_attrs.get("test_bacc"),
                    "test_auc": t.user_attrs.get("test_auc"), "params": hp})
    out = Path(a.study_dir) / "top_configs.json"
    json.dump(top, open(out, "w"), indent=2)
    print(f"[report] wrote {out}  ({len(top)} configs × {len(a.seeds_list)} seeds to refit)")


def _params_to_flags(arch, params):
    """Map an Optuna params dict back to trainer CLI flags for the given arch."""
    keys = (["lr", "weight_decay", "drop_path_rate", "attn_dropout", "label_smoothing", "augment"]
            if arch == "vit_mae" else
            ["lr", "weight_decay", "drop_rate", "label_smoothing", "augment"])
    return {f"--{k}": params[k] for k in keys if k in params}


def mode_refit(a):
    top = json.load(open(Path(a.study_dir) / "top_configs.json"))
    seeds = a.seeds_list
    k, s_idx = divmod(a.refit_index, len(seeds))
    if k >= len(top):
        print(f"[refit] index {a.refit_index} beyond topk×seeds — nothing to do."); return
    cfg, seed = top[k], seeds[s_idx]
    flags = _params_to_flags(a.arch, cfg["params"])
    out_dir = Path(a.winner_dir) / f"config_{k}"
    print(f"[refit] config #{k} (val_bACC {cfg['val_bacc']:.4f}) seed {seed} -> {out_dir}")
    rc = _run_trainer(a.arch, a, flags, seed=seed, out_dir=out_dir, keep_ckpt=True,
                      extra_env={"WANDB_MODE": "offline", "WANDB_DIR": str(out_dir)})
    sys.exit(rc)


def mode_winner(a):
    import numpy as np
    rows = []
    for mp in glob.glob(str(Path(a.winner_dir) / "config_*" / "**" / "metrics.json"), recursive=True):
        ck = mp.split(os.sep)
        cfg_k = next((p for p in ck if p.startswith("config_")), "config_?")
        d = json.load(open(mp)); c, tm = d.get("config", {}), d.get("test_metrics", {})
        auc = next((tm[k] for k in ("auc_roc_ovr", "auc_roc") if k in tm), None)
        rows.append(dict(config=cfg_k, seed=c.get("seed"),
                         val_bacc=c.get("best_val_balanced_acc"), val_auc=c.get("best_val_auc"),
                         test_bacc=tm.get("balanced_acc"), test_auc=auc, test_f1=tm.get("macro_f1", tm.get("f1"))))
    if not rows:
        print("[winner] no refit metrics.json found."); return
    import pandas as pd
    df = pd.DataFrame(rows)
    agg = df.groupby("config").agg(
        val_bacc_mean=("val_bacc", "mean"), val_bacc_std=("val_bacc", "std"),
        test_bacc_mean=("test_bacc", "mean"), test_bacc_std=("test_bacc", "std"),
        test_auc_mean=("test_auc", "mean"), test_auc_std=("test_auc", "std"),
        n=("seed", "count")).reset_index().sort_values("val_bacc_mean", ascending=False)
    out = Path(a.winner_dir) / "winner_summary.csv"
    df.sort_values(["config", "seed"]).to_csv(Path(a.winner_dir) / "refit_per_seed.csv", index=False)
    agg.to_csv(out, index=False)
    win = agg.iloc[0]
    print(f"\n[winner] {a.arch} {a.task} — WINNER {win['config']} "
          f"(mean val-bACC {win['val_bacc_mean']:.4f}):")
    print(f"    TEST bACC = {win['test_bacc_mean']:.4f} ± {win['test_bacc_std']:.4f}")
    print(f"    TEST AUC  = {win['test_auc_mean']:.4f} ± {win['test_auc_std']:.4f}   (paper 0.954)")
    print(f"[winner] per-config table -> {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", required=True, choices=["search", "report", "refit", "winner"])
    p.add_argument("--arch", required=True, choices=list(ARCHS))
    p.add_argument("--task", default="T2_multiclass")
    p.add_argument("--study-dir", default=None, help="Optuna study dir (journal + trials).")
    p.add_argument("--winner-dir", default=None, help="Refit winner-tree dir.")
    p.add_argument("--n-trials", type=int, default=1, help="Trials this agent runs (array → 1).")
    p.add_argument("--topk", type=int, default=3)
    p.add_argument("--refit-index", type=int, default=0)
    p.add_argument("--seeds", default="0,1,2")
    p.add_argument("--long", default="all")
    p.add_argument("--search-epochs", type=int, default=None,
                   help="Cap epochs during SEARCH only (speed). Refit uses the trainer default.")
    p.add_argument("--num_workers", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    p.add_argument("--matched_labels_csv", default=DEF_MASTER)
    p.add_argument("--data_dir", default=DEF_DATA)
    p.add_argument("--pretrained", default=None, help="Override arch default pretrained ckpt.")
    p.add_argument("--inputs-dir", default=None, help="Override arch default inputs dir.")
    p.add_argument("--wandb", action="store_true", default=True)
    p.add_argument("--no-wandb", dest="wandb", action="store_false")
    a = p.parse_args()
    a.seeds_list = [int(s) for s in a.seeds.split(",")]
    if a.mode in ("search", "report", "refit") and not a.study_dir:
        p.error("--study-dir required for search/report/refit")
    if a.mode in ("refit", "winner") and not a.winner_dir:
        p.error("--winner-dir required for refit/winner")
    {"search": mode_search, "report": mode_report,
     "refit": mode_refit, "winner": mode_winner}[a.mode](a)


if __name__ == "__main__":
    main()
