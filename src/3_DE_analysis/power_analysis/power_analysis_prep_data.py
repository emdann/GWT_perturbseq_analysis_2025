import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm.notebook import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from copy import deepcopy

# Add the parent directory to the path to import from sibling directory
sys.path.append(os.path.abspath('../../1_preprocess/'))
sys.path.append(os.path.abspath('../../'))
from preprocess import _convert_oak_path

from tqdm.notebook import tqdm

def merge_pseudobulks(h5ad_files, sample_cols = ['cell_sample_id', 'guide_id'], selected_perturbed_genes = None, NTC_control_label = 'NTC'):
    # Get var table for all genes
    all_vars = [anndata.experimental.read_lazy(a).var.to_dataframe().copy() for a in h5ad_files]
    var_merged = pd.concat(all_vars)[['gene_ids', 'gene_name']].drop_duplicates()

    all_obs = []
    keep_cols = ['10xrun_id', 'cell_sample_id', 'donor_id', 'culture_condition',
        'guide_id', 'perturbed_gene_name', 'perturbed_gene_id', 'guide_type']
    for a in h5ad_files:
        a_obs = anndata.experimental.read_lazy(a).obs.to_dataframe()[keep_cols]
        a_obs["sample_id"] = a_obs[sample_cols].apply(lambda x: "_".join(x), axis=1)
        all_obs.append(a_obs)
    obs_merged = pd.concat(all_obs).drop_duplicates()
    obs_merged = obs_merged.set_index('sample_id')
    if selected_perturbed_genes is not None:
        if NTC_control_label not in selected_perturbed_genes:
            selected_perturbed_genes.append(NTC_control_label)
        obs_merged = obs_merged[obs_merged.perturbed_gene_name.isin(selected_perturbed_genes)]

    # Read and sum pseudobulks
    merged_adata = sc.read_h5ad(h5ad_files[0])
    if selected_perturbed_genes is not None:
        merged_adata = merged_adata[merged_adata.obs.perturbed_gene_name.isin(selected_perturbed_genes)]
    merged_adata.obs["sample_id"] = merged_adata.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)
    n_cells_merged = merged_adata.obs[['sample_id', 'n_cells']].set_index('sample_id')

    for f in tqdm(h5ad_files[1:], desc="Merging pseudobulk files"):
        a = sc.read_h5ad(f)
        if selected_perturbed_genes is not None:
            a = a[a.obs.perturbed_gene_name.isin(selected_perturbed_genes)]
        a.obs["sample_id"] = a.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)
        # Keep track of n_cells
        n_cells_a = a.obs[['sample_id', 'n_cells']].set_index('sample_id')
        all_indices = n_cells_merged.index.union(n_cells_a.index)
        n_cells_merged = n_cells_merged.reindex(all_indices).fillna(0)
        n_cells_a = n_cells_a.reindex(all_indices).fillna(0)
        n_cells_merged = n_cells_merged + n_cells_a
        merged_adata = anndata.concat([merged_adata, a], join='outer')
        merged_adata = sc.get.aggregate(merged_adata, by='sample_id', func='sum', layer='sum')

    merged_adata.var = var_merged.loc[merged_adata.var_names]
    merged_adata.obs = obs_merged.loc[merged_adata.obs_names]
    merged_adata.obs['n_cells'] = n_cells_merged.loc[merged_adata.obs_names]
    return(merged_adata)


def merge_pseudobulk_from_list(adatas, sample_cols = ['cell_sample_id', 'guide_id']):
    for a in adatas:
        a.obs["sample_id"] = a.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)

    all_vars = [a.var for a in adatas]
    var_merged = pd.concat(all_vars)[['gene_ids', 'gene_name']].drop_duplicates()
    all_obs = [a.obs.drop('n_cells', axis=1) for a in adatas]
    obs_merged = pd.concat(all_obs).drop_duplicates()
    obs_merged = obs_merged.set_index('sample_id')
    
    # Keep track of n_cells
    n_cells_merged = pd.concat([a.obs[['n_cells']] for a in adatas], axis=1).fillna(0).sum(1)
    merged_adata = anndata.concat(adatas, join='outer')
    merged_adata = sc.get.aggregate(merged_adata, by='sample_id', func='sum', layer='sum')

    merged_adata.var = var_merged.loc[merged_adata.var_names]
    merged_adata.obs = obs_merged.loc[merged_adata.obs_names]
    merged_adata.obs['n_cells'] = n_cells_merged.loc[merged_adata.obs_names]
    return(merged_adata)

if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Prepare data for power analysis')
    parser.add_argument('--sample_id', type=str, required=True,
                       help='Sample ID to analyze')

    args = parser.parse_args()

    sample_id = args.sample_id

    # Read config
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)

    datadir = _convert_oak_path(config['sample_ids'][sample_id]['datadir'])
    random_seeds = config['random_seeds']
    selected_perturbed_genes = config['selected_perturbed_genes']
    subset_size = config['subset_size']

    # Get pseudobulk h5ad files
    tmp_dir = os.path.join(datadir, 'tmp')
    h5ad_files = [f'{tmp_dir}/{f}' for f in os.listdir(tmp_dir) if f.endswith('.postQC.DE_pseudobulk.h5ad')]
    sample_h5ad_files = [x for x in h5ad_files if sample_id in x]
    power_analysis_adatas = {}

    for s in random_seeds:
        np.random.seed(s)

        # Split files into non-overlapping sets of ~n lanes
        all_indices = np.arange(len(sample_h5ad_files))
        np.random.shuffle(all_indices)
        split_indices = [all_indices[i:i+subset_size] for i in range(0, len(all_indices), subset_size)]

        merged_adatas = []
        for ixs in split_indices:
            subset_sample_h5ad_files = [sample_h5ad_files[i] for i in ixs] # Merge set of 5 lanes
            merged_adata = merge_pseudobulks(subset_sample_h5ad_files, selected_perturbed_genes = selected_perturbed_genes)
            merged_adatas.append(merged_adata)

        # Sum different number of lanes
        power_analysis_adatas[f'valsplit_seed{s}'] = merge_pseudobulk_from_list(merged_adatas[0:1])
        for i in np.arange(2, len(merged_adatas)+1):
            power_analysis_adatas[f'split{i}_seed{s}'] = merge_pseudobulk_from_list(merged_adatas[1:i])

    # Store objects
    power_analysis_mdata = md.MuData(power_analysis_adatas)
    power_analysis_mdata.write_h5mu(f'./power_analysis.{sample_id}.h5mu')
