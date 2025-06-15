'''
To run on comino:
for batch in {1..273}; do sbatch --gres=gpu:1 --mem=10G --wrap="python K562_DE/k562_make_pseudobulk.py --batch ${batch}"; done
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
    parser = argparse.ArgumentParser(description='Generate pseudobulk data for K562 batch')
    parser.add_argument('--batch', type=int, required=True, help='Batch index to process')
    return parser.parse_args()

def main():
    args = parse_args()
    b = args.batch

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

    ## -- Load Replogle data -- ##
    print(f"Loading K562 data for batch {b}...")
    load_annotation_index = True
    k562_adata = anndata.experimental.read_lazy('/mnt/oak/users/emma/data/wesvae-data/K562_gwps_full.h5ad', load_annotation_index = True)
    batch_adata = k562_adata[k562_adata.obs['batch'] == b]
    print(f"Loaded batch data with shape: {batch_adata.shape}")

    # Filter to common tests (targets and genes)
    print("Filtering to common genes and targets...")
    all_targets = batch_adata.obs['gene'].to_dataframe()['gene'].unique()
    k562_gs = batch_adata.var['ensembl_id'].compute()
    common_gs = np.intersect1d(adata_de.var['gene_ids'], k562_gs)
    print(f"Found {len(common_gs)} common genes")

    control_level = 'non-targeting'
    common_targets = np.intersect1d(all_targets, keep_targets + [control_level]).tolist()
    print(f"Found {len(common_targets)} common targets")

    print("Subsetting data to common features...")
    batch_adata = batch_adata[:, k562_gs.isin(common_gs)]
    batch_adata = batch_adata[batch_adata.obs['gene'].isin(common_targets)].copy()
    print(f"Filtered data shape: {batch_adata.shape}")

    ## -- Subset and pseudobulk -- ##
    print("Moving data to GPU and computing pseudobulk...")
    SPARSE_CHUNK_SIZE = 100_000
    batch_adata.X = batch_adata.X.map_blocks(lambda x: x.toarray(), dtype=batch_adata.X.dtype, meta=np.array([]))
    batch_adata.X = batch_adata.X.rechunk((SPARSE_CHUNK_SIZE, batch_adata.shape[1]))
    batch_adata = batch_adata.to_memory()
    rapids_singlecell.get.anndata_to_GPU(batch_adata) 
    pbulk_batch_adata = rapids_singlecell.get.aggregate(batch_adata, by=['batch', 'gene'], func='sum')
    print(f"Pseudobulk data shape: {pbulk_batch_adata.shape}")
    
    # Save pseudobulk data for this batch
    print(f"Saving pseudobulk data for batch {b}...")
    pbulk_batch_adata.write_h5ad(f'/mnt/oak/users/emma/data/GWT/K562_gwps_filt_pseudobulk_batch_{b}.h5ad')
    print("Done!")

if __name__ == '__main__':
    main()