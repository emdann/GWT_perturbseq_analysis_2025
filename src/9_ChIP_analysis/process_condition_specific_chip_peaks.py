#!/usr/bin/env python
"""
Process condition-specific ChIP Atlas peak data for TRIM28 and KDM6B:
  1. Group tracks by TF + condition using the experiment manifest
  2. Compute a union consensus peak set per (TF, condition)
       - TRIM28: "resting" vs "stimulated"  (anti-CD3/CD28 + IL-2, 3 days; GSE81872)
       - KDM6B:  "resting" (T0) vs "stimulated_8h" (T8; CD3/CD28 beads; GSE122219)
         Note: SRX4985784 is excluded (H3K27me3 ChIP, mislabeled as KDM6B in ChIP Atlas)
  3. For each condition, overlap consensus peaks with gene promoters (TSS ± 2 kb)
  4. Save one binary TF × gene matrix per condition

Usage:
    source activate gwt-env
    python process_condition_specific_chip_peaks.py

Outputs (in OUTDIR):
    tf_promoter_binary_matrix.<condition>.csv / .parquet
"""

import glob
import os

import bioframe as bf
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths / settings
# ---------------------------------------------------------------------------
CHIP_DIR = "/mnt/oak/users/emma/data/ChIP_atlas_data"
MANIFEST = (
    "/mnt/oak/users/emma/bin/GWT_perturbseq_analysis"
    "/src/9_ChIP_analysis/chipatlas_cd4_tf_experiments.tsv"
)
ANNOTATION_FILE = (
    "/mnt/oak/users/emma/bin/GWT_perturbseq_analysis"
    "/src/5_sgRNA_annotation/genome/gene_annotations_hg38.parquet"
)
OUTDIR = (
    "/mnt/oak/users/emma/bin/GWT_perturbseq_analysis"
    "/src/9_ChIP_analysis/results"
)
THRESHOLD = "20"        # q-value cutoff matching downloaded BED filenames
PROMOTER_WINDOW = 2000  # bp upstream/downstream of TSS

# TFs to process with condition-specific consensus; maps TF → conditions to include
TF_CONDITIONS = {
    "TRIM28": ["resting", "stimulated"],
    "KDM6B":  ["resting", "stimulated_8h"],
}

# SRX IDs to exclude despite appearing in the manifest under a TF name
EXCLUDE_SRX = {
    "SRX4985784",   # H3K27me3 ChIP mislabeled as KDM6B in ChIP Atlas
}

os.makedirs(OUTDIR, exist_ok=True)

BED_COLS = ["chrom", "start", "end"]


# ---------------------------------------------------------------------------
# Helper functions (same as process_chip_peaks.py)
# ---------------------------------------------------------------------------

def load_bed(filepath: str) -> pd.DataFrame:
    """Load first 3 columns of a BED file as a bioframe-compatible DataFrame."""
    return pd.read_csv(filepath, sep="\t", usecols=[0, 1, 2], names=BED_COLS,
                       dtype={"chrom": str, "start": int, "end": int})


def union_peaks(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge all peak DataFrames into a single union/consensus set."""
    combined = pd.concat([df[BED_COLS] for df in dfs], ignore_index=True)
    return bf.merge(combined)[BED_COLS]


# ---------------------------------------------------------------------------
# Load manifest and filter to relevant TFs/conditions
# ---------------------------------------------------------------------------
manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)
manifest = manifest[~manifest["srx_id"].isin(EXCLUDE_SRX)]
manifest = manifest[manifest["tf"].isin(TF_CONDITIONS)]
manifest = manifest[manifest["condition"].notna() & (manifest["condition"] != "NA")]

# Keep only the conditions we want per TF
rows_to_keep = []
for tf, conditions in TF_CONDITIONS.items():
    rows_to_keep.append(
        manifest[(manifest["tf"] == tf) & (manifest["condition"].isin(conditions))]
    )
manifest = pd.concat(rows_to_keep, ignore_index=True)

print(f"Manifest entries to process: {len(manifest)}")
print(manifest[["srx_id", "tf", "condition"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Step 1: Consensus peaks per (TF, condition)
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Step 1: Computing condition-specific consensus peaks")
print("=" * 60)

# { (tf, condition): DataFrame of consensus peaks }
condition_consensus: dict[tuple[str, str], pd.DataFrame] = {}

for (tf, condition), group in manifest.groupby(["tf", "condition"]):
    cache_file = os.path.join(
        CHIP_DIR, tf, f"{tf}_{condition}_consensus_peaks.bed"
    )

    if os.path.exists(cache_file):
        print(f"  {tf} [{condition}]: loading cached consensus peaks")
        consensus = pd.read_csv(cache_file, sep="\t", usecols=[0, 1, 2],
                                names=BED_COLS,
                                dtype={"chrom": str, "start": int, "end": int})
        condition_consensus[(tf, condition)] = consensus
        continue

    srx_ids = group["srx_id"].tolist()
    bed_files = [
        os.path.join(CHIP_DIR, tf, f"{srx}.{THRESHOLD}.bed")
        for srx in srx_ids
    ]
    missing = [f for f in bed_files if not os.path.exists(f)]
    if missing:
        print(f"  {tf} [{condition}]: WARNING — missing BED files: {missing}")
        bed_files = [f for f in bed_files if os.path.exists(f)]

    if not bed_files:
        print(f"  {tf} [{condition}]: no BED files found, skipping")
        continue

    print(f"  {tf} [{condition}]: {len(bed_files)} tracks → union consensus")
    dfs = [bf.merge(load_bed(f)) for f in bed_files]
    consensus = union_peaks(dfs)

    consensus.to_csv(cache_file, sep="\t", header=False, index=False)
    condition_consensus[(tf, condition)] = consensus
    print(f"    -> {len(consensus):,} peaks saved to {os.path.basename(cache_file)}")


# ---------------------------------------------------------------------------
# Step 2: Gene promoters from GENCODE hg38 annotations
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Step 2: Defining promoters (TSS ± {:,} bp)".format(PROMOTER_WINDOW))
print("=" * 60)

genes = pd.read_parquet(ANNOTATION_FILE)
print(f"  Loaded {len(genes):,} genes from GENCODE v48 (hg38)")

genes["tss"] = np.where(genes["strand"] == "+", genes["start"], genes["end"])
genes["prom_start"] = (genes["tss"] - PROMOTER_WINDOW).clip(lower=0)
genes["prom_end"] = genes["tss"] + PROMOTER_WINDOW

promoters = (
    genes[["chrom", "prom_start", "prom_end", "gene_id", "gene_name"]]
    .rename(columns={"prom_start": "start", "prom_end": "end"})
    .copy()
)
print(f"  Promoters defined for {len(promoters):,} genes")


# ---------------------------------------------------------------------------
# Step 3: Build binary matrix per condition
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Step 3: Building condition-specific TF × gene binary matrices")
print("=" * 60)

# Collect all conditions present
all_conditions = sorted({cond for (_, cond) in condition_consensus})

for condition in all_conditions:
    tfs_in_condition = sorted(
        tf for (tf, cond) in condition_consensus if cond == condition
    )
    if not tfs_in_condition:
        continue

    print(f"\n  Condition: {condition}  |  TFs: {tfs_in_condition}")

    binary_matrix = pd.DataFrame(
        0,
        index=genes["gene_id"],
        columns=tfs_in_condition,
        dtype=np.int8,
    )
    binary_matrix.index.name = "gene_id"

    for tf in tfs_in_condition:
        consensus = condition_consensus[(tf, condition)]
        if consensus.empty:
            print(f"    {tf}: 0 consensus peaks — all zeros")
            continue

        ov = bf.overlap(promoters, consensus, how="inner", suffixes=("", "_chip"))
        hit_gene_ids = ov["gene_id"].unique()
        binary_matrix.loc[hit_gene_ids, tf] = 1

        n_hits = int(binary_matrix[tf].sum())
        print(f"    {tf}: {n_hits:,} genes with promoter peak")

    binary_matrix.insert(0, "gene_name", genes.set_index("gene_id")["gene_name"])

    out_csv = os.path.join(OUTDIR, f"tf_promoter_binary_matrix.{condition}.csv")
    out_parquet = os.path.join(OUTDIR, f"tf_promoter_binary_matrix.{condition}.parquet")
    binary_matrix.to_csv(out_csv)
    binary_matrix.to_parquet(out_parquet)

    print(f"    Saved: {os.path.basename(out_csv)} ({binary_matrix.shape[0]:,} genes × {len(tfs_in_condition)} TFs)")

print("\n" + "=" * 60)
print("Done.")
print("=" * 60)
