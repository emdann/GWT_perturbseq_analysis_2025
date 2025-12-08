#!/usr/bin/env python3
"""
Merge technical replicates for a biological sample and save to upload directory.

Usage:
    python merge_sample_for_upload.py --donor D1 --condition Rest
"""

import argparse
import anndata
import pandas as pd
import yaml
import glob
import os
import sys
from collections import defaultdict

# Add path to import custom modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../1_preprocess/')))
from qc_samples import *
from preprocess import _convert_oak_path


def parse_sample_upload(f, datadir, sgrna_library_metadata, sample_id, guide_group_only=False):
    """Parse and process a single sample."""
    adata = anndata.experimental.read_lazy(f)
    if not isinstance(adata.obs, pd.DataFrame):
        adata.obs = adata.obs.to_dataframe()
    sgrna_assignments2adata(adata, datadir, sgrna_library_metadata=sgrna_library_metadata, sample_id=sample_id)

    # Get PuroR expression for all cells before filtering
    puro_expr = adata[:, 'CUSTOM001_PuroR'].X.compute()
    adata.obs['PuroR'] = puro_expr.toarray().flatten()

    # Calculate guide group stats (before filtering)
    adata.obs['guide_group'] = 'targeting single sgRNA'
    adata.obs['guide_group'] = np.where(adata.obs['guide_id'].isna(), 'no sgRNA', adata.obs['guide_group'])
    adata.obs['guide_group'] = np.where(adata.obs['guide_id'] == 'multi_sgRNA', 'multi sgRNA', adata.obs['guide_group'])

    qc_summary = get_qc_summary(adata)

    return adata


def extract_bio_sample(sample_id):
    """Extract biological sample prefix (e.g., 'D1_Rest', 'D2_Stim8hr')"""
    parts = sample_id.split('_')
    # Pattern is typically: CD4i_R1_D1_Rest_... or CD4i_R1_D2_Stim8hr_...
    # Find the donor (D1, D2) and condition (Rest, Stim8hr)
    donor = None
    condition = None
    for i, part in enumerate(parts):
        if part.startswith('D') and part[1:].isdigit():
            donor = part
        elif part in ['Rest', 'Stim8hr', 'Stim48hr']:
            condition = part
            break

    if donor and condition:
        return f"{donor}_{condition}"
    return None


def main():
    parser = argparse.ArgumentParser(description='Merge technical replicates for a biological sample')
    parser.add_argument('--donor', required=True, help='Donor ID (e.g., D1, D2, D3, D4)')
    parser.add_argument('--condition', required=True, help='Condition (e.g., Rest, Stim8hr, Stim48hr)')
    parser.add_argument('--outdir', default='/mnt/oak/users/emma/data/GWT/to_share/',
                        help='Output directory')
    parser.add_argument('--config', default='../../metadata/experiments_config.yaml',
                        help='Path to experiments config YAML')
    parser.add_argument('--sgrna-lib', default='../../metadata/sgRNA_library_curated.csv',
                        help='Path to sgRNA library metadata CSV')

    args = parser.parse_args()

    # Construct biological sample ID
    bio_sample = f"{args.donor}_{args.condition}"
    print(f"Processing biological sample: {bio_sample}")

    # Load configuration
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_file = os.path.join(script_dir, args.config)
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    # Load sgRNA library metadata
    sgrna_lib_path = os.path.join(script_dir, args.sgrna_lib)
    sgrna_library_metadata = pd.read_csv(sgrna_lib_path, index_col=0)

    # Step 1: Collect all samples matching the biological sample ID
    matching_samples = []

    for exp in ['CD4iR1_Psomagen', 'CD4iR2_Psomagen']:
        exp_config = config[exp]
        datadir = _convert_oak_path(exp_config['datadir'])
        sample_metadata_csv = _convert_oak_path(exp_config['sample_metadata'])

        sample_metadata = pd.read_csv(sample_metadata_csv, index_col=0)

        scRNA_files = glob.glob(os.path.join(datadir, 'tmp', '*.scRNA.h5ad'))

        for f in scRNA_files:
            sample_id = os.path.basename(f).replace('.scRNA.h5ad', '')
            sample_bio_id = extract_bio_sample(sample_id)

            if sample_bio_id == bio_sample:
                matching_samples.append({
                    'sample_id': sample_id,
                    'file_path': f,
                    'datadir': datadir,
                    'exp': exp
                })

    if not matching_samples:
        print(f"ERROR: No samples found for {bio_sample}")
        sys.exit(1)

    print(f"Found {len(matching_samples)} technical replicates for {bio_sample}")

    # Step 2: Load and process all matching samples
    adatas = []
    for i, sample_info in enumerate(matching_samples, 1):
        sample_file = sample_info['file_path']
        datadir = sample_info['datadir']
        sample_id = sample_info['sample_id']

        print(f"  [{i}/{len(matching_samples)}] Loading {sample_id}...")
        adata = parse_sample_upload(sample_file, datadir, sgrna_library_metadata, sample_id)

        # Convert to DataFrame if needed
        if not isinstance(adata.var, pd.DataFrame):
            adata.var = adata.var.to_dataframe()
        if not isinstance(adata.obs, pd.DataFrame):
            adata.obs = adata.obs.to_dataframe()

        # Move counts to X
        adata.X = adata.layers['counts'].copy()
        del adata.layers['counts']
        cols2drop = [
            'library_id', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes',
            'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes',
            'pct_counts_in_top_500_genes', 'total_counts_mt', 'sequence',
            'log1p_total_counts_mt', 'n_genes'
        ]
        adata.obs = adata.obs.drop(cols2drop, axis=1)

        adatas.append(adata)

    # Step 3: Concatenate all samples
    print(f"Concatenating {len(adatas)} samples...")
    merged_adata = anndata.concat(adatas, join='outer', merge='same')
    
    # # Sort by perturbation name and lane ID
    # sorted_obs_names = merged_adata.obs[['lane_id', 'perturbed_gene_name']].sort_values(by=['perturbed_gene_name', 'lane_id']).index.values

    # Step 4: Write to output directory
    os.makedirs(args.outdir, exist_ok=True)
    output_file = os.path.join(args.outdir, f'{bio_sample}.assigned_guide.h5ad')
    print(f"Saving to {output_file}...")
    merged_adata.write_h5ad(output_file)

    print(f"SUCCESS: Saved {bio_sample} with {merged_adata.n_obs} cells")


if __name__ == '__main__':
    main()
