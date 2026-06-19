r"""
23c_cox_survival_resolved.py   (LOCAL, env: snp)
=================================================
Cox proportional-hazards survival models for the 8 MAF-resolved PRS
combinations exported by `25d_export_maf_resolved_sets.py` to
`D:\…\GWAS_filtered_maf_resolved\<COMBO>\`. Parallels script 23 but
consumes the new self-contained `(snps.tsv + patient_dosage.tsv)` trio
per combo.

Three transitions (same as script 23):
  CN_to_MCI    at-risk bl_dx==CN     event ever_conversion_MCI
  MCI_to_AD    at-risk bl_dx==MCI    event progressive_MCI
  CNMCI_to_AD  at-risk bl_dx∈{CN,MCI} event ever_conversion_AD

Per combo, two PRS/dosage atoms are added to the joint atom set:
  betaPRS_<combo>   weighted PRS  Σ β·dosage (1 column)
  dosage_<combo>    raw per-SNP dosage block (collapsed in fit_eval)

Atom registry: age, sex, APOE4, APOE2, betaPRS_<combo>, dosage_<combo>
(per combo). `covsets` enforces ≤1 PRS/dosage signature per model so
combos are never mixed. APOE4/APOE2 + a non-APOE PRS is allowed (all
combos here are noAPOE by construction post-MAF-resolution).

Output schema matches script 23 so script 23b can render it:
  <out>/<transition>/<covset>/seed_<n>/metrics.json

Usage:
  conda run -n snp python snp_pipeline/23c_cox_survival_resolved.py \
    --out "D:/ADNI_SNP_Omni2.5M_20140220/cox_survival_resolved"
"""
from __future__ import annotations

import argparse
import importlib.util as _il
import itertools
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── importlib reuse of scripts 20 + 23 ───────────────────────────────────────
_S20 = Path(__file__).parent / "20_genetic_baselines.py"
_spec = _il.spec_from_file_location("s20", _S20)
_s20 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s20)
_covars = _s20._covars
SPLITS = _s20.SPLITS

_S23 = Path(__file__).parent / "23_cox_survival_models.py"
_spec = _il.spec_from_file_location("s23", _S23)
_s23 = _il.module_from_spec(_spec)
_spec.loader.exec_module(_s23)
build_frame = _s23.build_frame
fit_eval = _s23.fit_eval
split_ids = _s23.split_ids
TRANSITIONS = _s23.TRANSITIONS
LABELS_TSV = _s23.LABELS_TSV
SPLIT_ROOT = _s23.SPLIT_ROOT

BASE = Path("D:/ADNI_SNP_Omni2.5M_20140220")
COMBO_ROOT_DEFAULT = BASE / "GWAS_filtered_maf_resolved"
COMBOS = ["sourceB", "sourceW", "sourceD", "sourceK",
          "5W26B", "5W26B14D", "5W26B13D", "5W26B14D2K"]

CLIN_ATOMS = ["age", "sex", "APOE4", "APOE2"]
APOE_GENE = {"APOE4", "APOE2"}


def covsets_resolved(combos: list[str]) -> list[list[str]]:
    """Non-empty atom subsets with ≤1 PRS/dosage signature per model.
    No APOE-PRS exclusion is needed because all resolved combos are
    noAPOE-equivalent (the published APOE-region SNPs are inert post-MAF)."""
    sig_atoms = [f"betaPRS_{c}" for c in combos] + [f"dosage_{c}" for c in combos]
    atoms = CLIN_ATOMS + sig_atoms
    sig = set(sig_atoms)
    out = []
    for k in range(1, len(atoms) + 1):
        for c in itertools.combinations(atoms, k):
            s = set(c)
            if len(s & sig) > 1:
                continue
            out.append(list(c))
    return out


def build_cov_matrix(combos: list[str], combo_root: Path, seed: int) -> pd.DataFrame:
    """Build a per-PTID cov matrix with clinical atoms + one betaPRS_<combo>
    column + one dos__<combo>__<rsID> column per SNP per combo."""
    parts = [pd.read_csv(SPLIT_ROOT / f"seed_{seed}" / f"{sp}.csv")
             for sp in ("train", "val", "test")]
    raw = pd.concat(parts, ignore_index=True).drop_duplicates("Patient_ID")
    cv = _covars(raw)
    cv.index = cv.index.astype(str)

    m = pd.DataFrame(index=cv.index)
    m["age"] = cv["AGE"]
    m["sex"] = cv["SEX"]
    m["APOE4"] = cv["APOE4"]
    m["APOE2"] = cv["APOE2"]

    for combo in combos:
        d = combo_root / combo
        snp = pd.read_csv(d / f"{combo}_snps.tsv", sep="\t").set_index("rsID")
        dos = pd.read_csv(d / f"{combo}_patient_dosage.tsv", sep="\t")
        dos["PTID"] = dos["PTID"].astype(str)
        dos = dos.set_index("PTID")
        dos = dos.apply(lambda c: c.fillna(c.mean()))
        cols = [s for s in snp.index if s in dos.columns]
        w = snp.loc[cols, "beta_A1"].astype(float).values
        prs = dos[cols].astype(float).mul(w, axis=1).sum(axis=1)
        prs.index = prs.index.astype(str)
        m[f"betaPRS_{combo}"] = prs.reindex(m.index)
        # raw-dosage block
        dblk = dos[cols].copy()
        dblk.index = dblk.index.astype(str)
        dblk = dblk.reindex(m.index)
        dblk.columns = [f"dos__{combo}__{c}" for c in cols]
        m = pd.concat([m, dblk], axis=1)
    return m


def covset_name(atoms: list[str]) -> str:
    """Stable directory name for an atom subset."""
    return "+".join(atoms)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combo-root", type=Path, default=COMBO_ROOT_DEFAULT)
    ap.add_argument("--combos", nargs="+", default=COMBOS)
    ap.add_argument("--out", type=Path,
                     default=BASE / "cox_survival_resolved")
    ap.add_argument("--labels-tsv", type=Path, default=LABELS_TSV)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--skip-if-exists", action="store_true")
    ap.add_argument("--clinical-only-once", action="store_true", default=True,
                     help="emit clinical-only covsets (age/sex/APOE4/APOE2 "
                          "subsets) once per transition, not per combo.")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(args.labels_tsv, sep="\t")
    lab["Patient_ID"] = lab["Patient_ID"].astype(str)

    print(f"[setup] out={args.out}  combos={args.combos}  seeds={args.seeds}")
    cov = {sd: build_cov_matrix(args.combos, args.combo_root, sd)
            for sd in args.seeds}
    cs = covsets_resolved(args.combos)
    # Deduplicate clinical-only covsets (avoid emitting them once per combo) —
    # they don't depend on any combo. Keep all combos × signatures.
    print(f"[grid] {len(TRANSITIONS)} transitions × {len(cs)} covsets × "
          f"{len(args.seeds)} seeds = "
          f"{len(TRANSITIONS) * len(cs) * len(args.seeds)} fits")

    done = skipped = failed = 0
    for tn, (ar, ev, tc, fb) in TRANSITIONS.items():
        fr = build_frame(lab, ar, ev, tc, fb)
        print(f"\n=== {tn}: n={len(fr)}  events={int(fr['event'].sum())} ===")
        for atoms in cs:
            for sd in args.seeds:
                name = covset_name(atoms)
                out_dir = args.out / tn / name / f"seed_{sd}"
                mfile = out_dir / "metrics.json"
                if args.skip_if_exists and mfile.exists():
                    skipped += 1
                    continue
                try:
                    r = fit_eval(fr, cov[sd], atoms, split_ids(sd))
                    out_dir.mkdir(parents=True, exist_ok=True)
                    json.dump({"transition": tn, "covset": name,
                                "atoms": atoms, "seed": sd,
                                "n_atrisk": int(len(fr)),
                                "n_events": int(fr["event"].sum()),
                                **r}, open(mfile, "w"), indent=2)
                    if r.get("status") == "ok":
                        done += 1
                    else:
                        skipped += 1
                except Exception as e:
                    failed += 1
                    print(f"  [FAIL] {tn}/{name}/seed{sd}: "
                          f"{type(e).__name__}: {e}")

    print(f"\n[done] ok={done}  skipped/insufficient={skipped}  failed={failed} "
          f"→ {args.out}")


if __name__ == "__main__":
    main()
