#!/usr/bin/env python
"""
Process ChIP Atlas peak data for CD4+ T cell TFs:
  1. Compute consensus peaks per TF (intersection across all tracks)
  2. Define gene promoters as TSS ± PROMOTER_WINDOW from Ensembl/GENCODE hg38 annotations
  3. Overlap TF consensus peaks with promoters
  4. Save binary matrix: rows = genes (Ensembl ID), columns = TFs

Usage:
    source activate gwt-env
    python process_chip_peaks.py

Both ChIP Atlas data (hg38) and gene annotations (GENCODE v48, hg38) use the
same genome build, so no liftover is needed.
"""

import glob
import os

import bioframe as bf
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CHIP_DIR = "/mnt/oak/users/emma/data/ChIP_atlas_data"
ANNOTATION_FILE = (
    "/mnt/oak/users/emma/bin/GWT_perturbseq_analysis"
    "/src/5_sgRNA_annotation/genome/gene_annotations_hg38.parquet"
)
OUTDIR = (
    "/mnt/oak/users/emma/bin/GWT_perturbseq_analysis"
    "/src/9_ChIP_analysis/results"
)
THRESHOLD = "20"          # matches the downloaded BED files (q < 1e-20)
PROMOTER_WINDOW = 2000    # bp upstream and downstream of TSS

os.makedirs(OUTDIR, exist_ok=True)

BED_COLS = ["chrom", "start", "end"]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_bed(filepath: str) -> pd.DataFrame:
    """Load first 3 columns of a BED file as a bioframe-compatible DataFrame."""
    return pd.read_csv(filepath, sep="\t", usecols=[0, 1, 2], names=BED_COLS,
                       dtype={"chrom": str, "start": int, "end": int})


def intersect_two(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Return intervals that are present in both df1 and df2, trimmed to the
    overlapping coordinates and then merged.
    """
    ov = bf.overlap(df1[BED_COLS], df2[BED_COLS], how="inner", suffixes=("", "_"))
    ov["start"] = ov[["start", "start_"]].max(axis=1)
    ov["end"] = ov[["end", "end_"]].min(axis=1)
    valid = ov[ov["end"] > ov["start"]][BED_COLS].copy()
    if valid.empty:
        return valid
    return bf.merge(valid)[BED_COLS]  # drop the n_intervals column bf.merge adds


def consensus_peaks_for_tf(tf_dir: str, tf: str) -> pd.DataFrame | None:
    """
    Compute the intersection of all peak BED files for a given TF.
    Each BED in tf_dir is first merged (self-overlapping peaks collapsed),
    then iteratively intersected with the running consensus.
    Returns None if no BED files are found.
    """
    bed_files = sorted(glob.glob(os.path.join(tf_dir, f"*.{THRESHOLD}.bed")))
    if not bed_files:
        return None

    print(f"  {tf}: {len(bed_files)} tracks")
    dfs = [bf.merge(load_bed(f)) for f in bed_files]

    result = dfs[0][BED_COLS]
    for df in dfs[1:]:
        result = intersect_two(result, df[BED_COLS])
        if result.empty:
            print(f"    WARNING: empty consensus after intersection — no peaks survive")
            return result

    return result


# ---------------------------------------------------------------------------
# Step 1: Consensus peaks per TF
# ---------------------------------------------------------------------------
print("=" * 60)
print("Step 1: Computing consensus peaks per TF")
print("=" * 60)

tf_dirs = sorted(
    d for d in glob.glob(os.path.join(CHIP_DIR, "*")) if os.path.isdir(d)
)

tf_consensus: dict[str, pd.DataFrame] = {}

for tf_dir in tf_dirs:
    tf = os.path.basename(tf_dir)
    consensus_file = os.path.join(tf_dir, f"{tf}_consensus_peaks.bed")

    if os.path.exists(consensus_file):
        print(f"  {tf}: loading existing consensus peaks")
        consensus = pd.read_csv(consensus_file, sep="\t", usecols=[0, 1, 2],
                                names=BED_COLS,
                                dtype={"chrom": str, "start": int, "end": int})
        tf_consensus[tf] = consensus
        continue

    consensus = consensus_peaks_for_tf(tf_dir, tf)
    if consensus is None:
        print(f"  {tf}: no BED files found, skipping")
        continue
    if consensus.empty:
        print(f"  {tf}: 0 consensus peaks")
        continue

    consensus.to_csv(consensus_file, sep="\t", header=False, index=False)
    tf_consensus[tf] = consensus
    print(f"    -> {len(consensus)} consensus peaks saved to {os.path.basename(consensus_file)}")

print(f"\n{len(tf_consensus)} TFs with consensus peaks")


# ---------------------------------------------------------------------------
# Step 2: Gene promoters from Ensembl/GENCODE hg38 annotations
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Step 2: Defining promoters (TSS ± {:,} bp)".format(PROMOTER_WINDOW))
print("=" * 60)

genes = pd.read_parquet(ANNOTATION_FILE)
print(f"  Loaded {len(genes):,} genes from GENCODE v48 (hg38)")
print(f"  Gene types: {genes['gene_type'].value_counts().head(5).to_dict()}")

# TSS: start coordinate for + strand genes, end coordinate for - strand genes
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
# Step 3: Overlap consensus peaks with promoters → binary matrix
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Step 3: Building TF × gene binary matrix")
print("=" * 60)

tf_list = sorted(tf_consensus.keys())

# Index: gene_id (unique Ensembl ID); include gene_name as a column
binary_matrix = pd.DataFrame(
    0,
    index=genes["gene_id"],
    columns=tf_list,
    dtype=np.int8,
)
binary_matrix.index.name = "gene_id"

for tf in tf_list:
    consensus = tf_consensus[tf]
    if consensus.empty:
        continue

    ov = bf.overlap(promoters, consensus, how="inner", suffixes=("", "_chip"))
    hit_gene_ids = ov["gene_id"].unique()
    binary_matrix.loc[hit_gene_ids, tf] = 1

    n_hits = int(binary_matrix[tf].sum())
    print(f"  {tf}: {n_hits:,} genes with promoter peak")

# Add gene_name as the first column for readability
binary_matrix.insert(0, "gene_name", genes.set_index("gene_id")["gene_name"])


# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
out_csv = os.path.join(OUTDIR, "tf_promoter_binary_matrix.csv")
out_parquet = os.path.join(OUTDIR, "tf_promoter_binary_matrix.parquet")

binary_matrix.to_csv(out_csv)
binary_matrix.to_parquet(out_parquet)

print("\n" + "=" * 60)
print("Done.")
print(f"  Matrix shape: {binary_matrix.shape[0]:,} genes × {len(tf_list)} TFs")
print(f"  Saved CSV    : {out_csv}")
print(f"  Saved Parquet: {out_parquet}")
print("=" * 60)
