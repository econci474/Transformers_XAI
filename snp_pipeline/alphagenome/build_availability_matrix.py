r"""
build_availability_matrix.py
============================
Audit script. Queries AlphaGenome metadata and emits TWO parallel views
of what tracks would be included in each extraction run:

  Version A (inclusive)  — DLPFC + hippocampus + parietal lobe (bulk adult)
                            + all 5 astrocyte CL codes (any life stage)
                            + 5 scorer modalities + CHIP_HISTONE + CHIP_TF
  Version B (strict)     — DLPFC + hippocampus (bulk adult)
                            + only hippocampal astrocyte (CL:0002604)
                            + same scorer set

Outputs (under snp_pipeline/alphagenome/):
  alphagenome_track_availability__A.csv  + .png
  alphagenome_track_availability__B.csv  + .png
  alphagenome_track_availability__diff.csv  (A \ B)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).parent
API_KEY = (HERE / "alphagenome_api.txt").read_text(encoding="utf-8").strip()

from alphagenome.models import dna_client
md = dna_client.create(API_KEY).output_metadata(
    organism=dna_client.Organism.HOMO_SAPIENS).concatenate()
md = md.assign(
    gtex_tissue=md["gtex_tissue"].fillna(""),
    biosample_name=md["biosample_name"].fillna(""),
    ontology_curie=md["ontology_curie"].fillna(""),
    biosample_life_stage=md["biosample_life_stage"].fillna(""),
)
md["output_short"] = md["output_type"].astype(str).str.replace("OutputType.", "", regex=False)

# ─── Region + cell-type rosters per version ─────────────────────────────
VERSIONS = {
    "A": {
        "regions": {
            "UBERON:0009834": "DLPFC",
            "UBERON:0001954": "Hippocampus",
            "UBERON:0001872": "Parietal lobe",
        },
        "celltypes": {
            "CL:0000127": "astrocyte (generic)",
            "CL:0002603": "astrocyte of the cerebellum",
            "CL:0002604": "astrocyte of the hippocampus",
            "CL:0002605": "astrocyte of the cerebral cortex",
            "CL:0002606": "astrocyte of the spinal cord",
        },
    },
    "B": {
        "regions": {
            "UBERON:0009834": "DLPFC",
            "UBERON:0001954": "Hippocampus",
            "UBERON:0002305": "Layer of hippocampus",
        },
        "celltypes": {
            "CL:0002604": "astrocyte of the hippocampus",
        },
    },
}

OUTPUT_COLS = ["RNA_SEQ", "DNASE", "ATAC", "CAGE", "PROCAP",
                "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
                "CHIP_HISTONE", "CHIP_TF"]


def select_tracks(version: dict) -> pd.DataFrame:
    """Region rows must be adult bulk tissue; cell-type rows take any life stage."""
    keep_region = (md["ontology_curie"].isin(version["regions"]) &
                    (md["biosample_type"] == "tissue") &
                    md["biosample_life_stage"].isin(["adult", ""]))
    keep_cells = md["ontology_curie"].isin(version["celltypes"])
    sel = md[keep_region | keep_cells].copy()

    def label(row):
        if row["ontology_curie"] in version["regions"]:
            return f"{version['regions'][row['ontology_curie']]} (bulk tissue)"
        return version["celltypes"].get(row["ontology_curie"], row["biosample_name"])

    sel["row_label"] = sel.apply(label, axis=1)
    return sel


def build_matrix(sel: pd.DataFrame, version: dict) -> tuple[np.ndarray, list[str]]:
    # Fix row order: regions first (in their dict order), then cell types
    row_labels = ([f"{lab} (bulk tissue)" for lab in version["regions"].values()] +
                  list(version["celltypes"].values()))
    M = np.zeros((len(row_labels), len(OUTPUT_COLS)), dtype=int)
    for i, lab in enumerate(row_labels):
        sub = sel[sel["row_label"] == lab]
        for j, out in enumerate(OUTPUT_COLS):
            M[i, j] = int((sub["output_short"] == out).sum())
    return M, row_labels


def render_png(M: np.ndarray, row_labels: list[str], title: str, path: Path):
    fig, ax = plt.subplots(figsize=(12.0, 0.55 * len(row_labels) + 2.4))
    cmap = LinearSegmentedColormap.from_list("avail",
        ["#f7f7f7", "#fde9c8", "#f6b26b", "#cc4125"], N=256)
    vmax = max(M.max(), 1)
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v == 0:
                ax.text(j, i, "—", ha="center", va="center",
                        color="#aaaaaa", fontsize=9)
            else:
                ax.text(j, i, str(v), ha="center", va="center",
                        color=("white" if v > vmax * 0.55 else "#222"),
                        fontsize=9, weight="bold")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9.5)
    ax.set_xticks(range(len(OUTPUT_COLS)))
    ax.set_xticklabels(OUTPUT_COLS, rotation=35, ha="right", fontsize=9)
    ax.set_title(title, fontsize=10.5)
    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("n tracks", fontsize=9)
    # Total per row on the right
    for i, lab in enumerate(row_labels):
        total = int(M[i].sum())
        ax.text(len(OUTPUT_COLS) - 0.4, i, f"  Σ={total}",
                ha="left", va="center", fontsize=9, color="#333", weight="bold")
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


# ─── Build and save both versions ──────────────────────────────────────
all_sel = {}
for v_name, v_cfg in VERSIONS.items():
    sel = select_tracks(v_cfg)
    sel = sel[["ontology_curie", "biosample_name", "biosample_type",
                "biosample_life_stage", "row_label", "output_short",
                "data_source", "histone_mark", "transcription_factor", "name"]]
    csv_path = HERE / f"alphagenome_track_availability__{v_name}.csv"
    sel.to_csv(csv_path, index=False)
    M, row_labels = build_matrix(sel, v_cfg)
    title = (f"AlphaGenome track availability — Version {v_name}  "
              f"(total {len(sel)} tracks)\n"
              + ("DLPFC + Hippocampus + Parietal lobe bulk + all 5 astrocyte CL codes + CHIP scorers added"
                  if v_name == "A" else
                  "DLPFC + Hippocampus bulk + ONLY hippocampal astrocyte + CHIP scorers added"))
    png_path = HERE / f"alphagenome_track_availability__{v_name}.png"
    render_png(M, row_labels, title, png_path)
    print(f"[{v_name}] {len(sel)} tracks  →  {csv_path.name}, {png_path.name}")
    all_sel[v_name] = sel

# Diff A \ B
diff = all_sel["A"].merge(all_sel["B"]["name"], on="name", how="left", indicator=True)
diff = diff[diff["_merge"] == "left_only"].drop(columns=["_merge"])
diff.to_csv(HERE / "alphagenome_track_availability__diff_A_minus_B.csv", index=False)
print(f"[A\\B] {len(diff)} extra tracks in A  →  alphagenome_track_availability__diff_A_minus_B.csv")

# Per-row totals printed to stdout
print("\n=== Version A row totals ===")
for v_name, sel in all_sel.items():
    print(f"  [{v_name}]")
    for lab, n in sel.groupby("row_label").size().items():
        print(f"      {lab:<45s}  {n:3d} tracks")


# ─── Narrowed-AD-spec audit figure (2026-06-07 refresh) ─────────────────
# Three row colour classes:
#   GREEN — in this run (strict allow-list, region-anchored bulk-adult or
#           CL:0002604 hippocampal astrocyte).
#   RED   — user-requested but the (region, cell-type) combination has
#           ZERO tracks anywhere in AlphaGenome's metadata.
#   GREY  — present in AlphaGenome but not allow-listed (no region anchor /
#           in-vitro / embryonic / not adult bulk).
#
# Column groups (user spec 2026-06-07):
SPLICE_KEYS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
SCORER_GROUPS = {
    "RNA\nexpression": ("RNA_SEQ",),
    "Accessibility":   ("ATAC", "DNASE"),
    "TSS":             ("CAGE", "PROCAP"),
    "Splicing":        SPLICE_KEYS,
    "Histone\nmarks":  ("CHIP_HISTONE",),
    "TF binding":      ("CHIP_TF",),
}

# Region-grouped row roster (DLPFC block, Hippocampus block, then
# hippo-astro green row, then grey).
# Tuple format: (label, colour, curie, bulk_filter_dict_or_None)
DLPFC_BULK_FILTER = {"biosample_type": "tissue", "life_stage_in": ("adult", "")}
HIPPO_BULK_FILTER = {"biosample_type": "tissue", "life_stage_in": ("adult", "")}

REGION_ROWS = [
    # ── DLPFC region block ──
    ("DLPFC · bulk tissue (all cells)",            "green", "UBERON:0009834", DLPFC_BULK_FILTER),
    ("DLPFC · microglia",                          "red",   None, None),
    ("DLPFC · astrocyte",                          "red",   None, None),
    ("DLPFC · neuron",                             "red",   None, None),
    ("DLPFC · oligodendrocyte",                    "red",   None, None),
    ("DLPFC · oligodendrocyte precursor cell (OPC)", "red", None, None),
    # ── Hippocampus region block ──
    ("Hippocampus · bulk tissue (all cells)",      "green", "UBERON:0001954", HIPPO_BULK_FILTER),
    ("Layer of hippocampus · bulk tissue",         "green", "UBERON:0002305", HIPPO_BULK_FILTER),
    ("Hippocampus · microglia",                    "red",   None, None),
    ("Hippocampus · astrocyte",                    "red",   None, None),
    ("Hippocampus · neuron",                       "red",   None, None),
    ("Hippocampus · oligodendrocyte",              "red",   None, None),
    ("Hippocampus · oligodendrocyte precursor cell (OPC)", "red", None, None),
    # ── Sorted-cell green row ──
    ("Hippocampus · astrocyte of the hippocampus", "green", "CL:0002604",     None),
]

GREY_ROWS = [
    ("astrocyte (generic, no region)",            "CL:0000127"),
    ("glutamatergic neuron (in vitro)",           "CL:0000679"),
    ("neuronal stem cell (in vitro)",             "CL:0000047"),
    ("motor neuron",                              "CL:0000100"),
    ("OPC (generic, no region)",                  "CL:0002453"),
    ("parietal lobe (embryonic)",                 "UBERON:0001872"),
]

# Indices where horizontal separator lines should fall (after the named row)
_REGION_SEPARATOR_AFTER = [5, 12, 13]  # after DLPFC block, after Hippo block, after hippo-astro


def _count_group_tracks_for_curie(curie: str, group_keys: tuple,
                                    bulk_filter: dict | None = None) -> int:
    sub = md[md["ontology_curie"] == curie]
    if bulk_filter:
        if "biosample_type" in bulk_filter:
            sub = sub[sub["biosample_type"] == bulk_filter["biosample_type"]]
        if "life_stage_in" in bulk_filter:
            sub = sub[sub["biosample_life_stage"].isin(bulk_filter["life_stage_in"])]
    return int(sub["output_short"].isin(group_keys).sum())


def _render_narrowed_ad_spec_png(out_path):
    # Build (rows × scorer-groups) count matrix in region-grouped order
    rows = []   # (label, colour, count_row)
    group_keys_list = list(SCORER_GROUPS.values())
    group_names = list(SCORER_GROUPS.keys())

    for lab, colour, curie, bulk_filter in REGION_ROWS:
        if colour == "red" or curie is None:
            counts = [0] * len(group_keys_list)
        else:
            counts = [_count_group_tracks_for_curie(curie, g, bulk_filter)
                      for g in group_keys_list]
        rows.append((lab, colour, counts))

    for lab, curie in GREY_ROWS:
        counts = [_count_group_tracks_for_curie(curie, g) for g in group_keys_list]
        rows.append((lab, "grey", counts))

    M = np.array([r[2] for r in rows], dtype=int)
    row_labels = [r[0] for r in rows]
    row_colours = [r[1] for r in rows]

    colour_hex = {
        "green": "#2a8a3e",
        "red":   "#b03030",
        "grey":  "#5a5a5a",
    }

    fig, ax = plt.subplots(figsize=(11.5, 0.42 * len(rows) + 2.6))
    cmap = LinearSegmentedColormap.from_list("avail",
        ["#ffffff", "#fde9c8", "#f6b26b", "#cc4125"], N=256)
    vmax = max(M.max(), 1)
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v == 0:
                ax.text(j, i, "—", ha="center", va="center",
                        color="#bbbbbb", fontsize=10)
            else:
                ax.text(j, i, str(v), ha="center", va="center",
                        color=("white" if v > vmax * 0.55 else "#222"),
                        fontsize=10, weight="bold")

    # Row labels with class colour
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(row_labels, fontsize=10)
    for tick, c in zip(ax.get_yticklabels(), row_colours):
        tick.set_color(colour_hex[c])
        if c == "green":
            tick.set_weight("bold")

    ax.set_xticks(range(len(group_names)))
    ax.set_xticklabels(group_names, rotation=0, ha="center", fontsize=10.5)
    ax.tick_params(axis="x", which="both", pad=4)

    # Per-row totals — show as an annotation OUTSIDE the axes box on the
    # right, so the heatmap width stays unchanged.
    for i in range(M.shape[0]):
        total = int(M[i].sum())
        ax.annotate(f"Σ={total}",
                     xy=(1.0, i), xycoords=("axes fraction", "data"),
                     xytext=(6, 0), textcoords="offset points",
                     ha="left", va="center", fontsize=9.5,
                     color=colour_hex[row_colours[i]], weight="bold",
                     annotation_clip=False)

    # Section separators: after each region block + before grey block
    n_region = len(REGION_ROWS)
    for sep_after in _REGION_SEPARATOR_AFTER + [n_region - 1]:
        ax.axhline(sep_after + 0.5, color="#999999", linewidth=0.8)

    ax.set_title("AlphaGenome track availability", fontsize=12)

    # Restore black border around the heatmap (the previous render had this
    # automatically via the colorbar's parent axes; now we set it explicitly).
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(1.0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


_NARROWED_PNG = HERE / "alphagenome_track_availability.png"
_render_narrowed_ad_spec_png(_NARROWED_PNG)
_n_green = sum(1 for _, c, _, _ in REGION_ROWS if c == "green")
_n_red   = sum(1 for _, c, _, _ in REGION_ROWS if c == "red")
print(f"\n[narrowed-AD-spec figure] {_NARROWED_PNG.name}  "
      f"({_n_green} green / {_n_red} red / {len(GREY_ROWS)} grey rows, "
      f"{len(SCORER_GROUPS)} scorer-group columns)")
