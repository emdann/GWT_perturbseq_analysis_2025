'''
To run on comino:
for gem_group in {1..56}; do sbatch --gres=gpu:1 --mem=10G --wrap="PYTHONNOUSERSITE=1 /home/emmadann/miniforge3/envs/gwt_de/bin/python Jurkat_DE/jurkat_make_pseudobulk.py --gem_group ${gem_group}"; done
'''
import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import rapids_singlecell

import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate pseudobulk data for Jurkat gem_group')
    parser.add_argument('--gem_group', type=int, required=True, help='Gem group index to process')
    return parser.parse_args()

def main():
    args = parse_args()
    g = args.gem_group

    ## -- Load output of DE analysis in T cells -- ##
    print(f"Loading DE analysis data for T cells ...")
    datadir = '/mnt/oak/users/emma/data/GWT/CD4iR1_Psomagen/'
    experiment_name = 'CD4iR1_Psomagen'
    adata_de = sc.read_h5ad(datadir + f'/DE_results_all_confounders/{experiment_name}.merged_DE_results.h5ad')
    print(f"Loaded DE data with shape: {adata_de.shape}")

    print("Computing z-scores...")
    adata_de.layers['zscore'] = adata_de.layers['log_fc'] / adata_de.layers['lfcSE']
    adata_de.layers['zscore'][np.where(adata_de.layers['zscore'] > 50)] = 50
    adata_de.var_names = adata_de.var['gene_name'].values

    # Load summary stats
    DE_stats = pd.read_csv(datadir + f'/DE_results_all_confounders/DE_summary_stats_per_target.csv', index_col=0)
    # Keep targets with significant effects in any condition
    keep_targets = DE_stats[DE_stats['ontarget_significant'] & (DE_stats.n_total_de_genes > 1)].target_name.unique().tolist()
    print(len(keep_targets))

    ## -- Load Jurkat data -- ##
    print(f"Loading Jurkat data for gem_group {g}...")
    jurkat_adata = anndata.experimental.read_lazy('/mnt/oak/users/emma/data/wesvae-data/Jurkat_essential.h5ad', load_annotation_index=True)
    batch_adata = jurkat_adata[jurkat_adata.obs['gem_group'] == g]
    print(f"Loaded gem_group data with shape: {batch_adata.shape}")

    # Filter to common tests (targets and genes)
    print("Filtering to common genes and targets...")
    all_targets = batch_adata.obs['gene'].to_dataframe()['gene'].unique()
    # var_names are already Ensembl IDs
    jurkat_gs = pd.Index(batch_adata.var_names)
    common_gs = np.intersect1d(adata_de.var['gene_ids'], jurkat_gs)
    print(f"Found {len(common_gs)} common genes")

    control_level = 'non-targeting'
    common_targets = np.intersect1d(all_targets, keep_targets + [control_level]).tolist()
    print(f"Found {len(common_targets)} common targets")

    print("Subsetting data to common features...")
    batch_adata = batch_adata[:, jurkat_gs.isin(common_gs)]
    batch_adata = batch_adata[batch_adata.obs['gene'].isin(common_targets)].copy()
    print(f"Filtered data shape: {batch_adata.shape}")

    ## -- Subset and pseudobulk -- ##
    print("Moving data to GPU and computing pseudobulk...")
    SPARSE_CHUNK_SIZE = 100_000
    batch_adata.X = batch_adata.X.map_blocks(lambda x: x.toarray(), dtype=batch_adata.X.dtype, meta=np.array([]))
    batch_adata.X = batch_adata.X.rechunk((SPARSE_CHUNK_SIZE, batch_adata.shape[1]))
    batch_adata = batch_adata.to_memory()
    rapids_singlecell.get.anndata_to_GPU(batch_adata)
    pbulk_batch_adata = rapids_singlecell.get.aggregate(batch_adata, by=['gem_group', 'gene'], func='sum')
    print(f"Pseudobulk data shape: {pbulk_batch_adata.shape}")

    # Save pseudobulk data for this gem_group
    print(f"Saving pseudobulk data for gem_group {g}...")
    pbulk_batch_adata.write_h5ad(f'/mnt/oak/users/emma/data/GWT/Jurkat_DE_analysis/tmp/Jurkat_essential_filt_pseudobulk_gemgroup_{g}.h5ad')
    print("Done!")

if __name__ == '__main__':
    main()
