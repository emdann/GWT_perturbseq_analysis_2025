import os,sys
import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import glob
from tqdm import tqdm

from copy import deepcopy
import argparse

def parse_DE_results_2_adata(df):
    all_dfs = {}
    for stat in ['baseMean', 'log_fc', 'lfcSE', 'p_value','adj_p_value']:
        stat_df = df.pivot(values=stat, columns='variable', index='target_contrast')
        all_dfs[stat] = stat_df

    DE_anndata = anndata.AnnData(
        layers = all_dfs.copy()
    )

    DE_anndata.obs_names = all_dfs['log_fc'].index.tolist()
    DE_anndata.var_names = all_dfs['log_fc'].columns.tolist()
    DE_anndata.obs = df[['target_contrast', 'target_contrast_gene_name', 'culture_condition']].drop_duplicates().set_index('target_contrast').loc[DE_anndata.obs_names]
    DE_anndata.obs['target_contrast'] = DE_anndata.obs_names.values
    DE_anndata.obs_names = DE_anndata.obs['target_contrast'] + '_' + DE_anndata.obs['culture_condition']
    return(DE_anndata)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare data for differential expression analysis')
    parser.add_argument('--experiment_name', type=str, required=True, help='Name of the experiment')
    parser.add_argument('--datadir', type=str, default='/mnt/oak/users/emma/data/GWT/', help='Directory containing the GWT data files')
    args = parser.parse_args()

    # Extract parameters from command line arguments
    experiment_name = args.experiment_name
    datadir = args.datadir + experiment_name + '/'

    # Read cell-level metadata
    obs_df = pd.read_csv(f'{datadir}/{experiment_name}_merged.gex.lognorm.postQC_obs.csv', compression='gzip', index_col=0)
    gene_name_to_id = dict(zip(obs_df['perturbed_gene_id'], obs_df['perturbed_gene_name']))
    var_df = sc.read_h5ad(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad', backed=True).var.copy()

    de_results_dir = datadir + '/DE_results/tmp/'

    try:
        combined_de_adata = sc.read_h5ad(datadir + f'/DE_results/{experiment_name}.merged_DE_results.h5ad')
    except:
        # Read all csv.gz files from the DE results directory
        de_results_files = glob.glob(de_results_dir + '*.csv.gz')
        de_results_adatas = []

        for file in tqdm(de_results_files, desc="Processing DE result files"):
            df = pd.read_csv(file, compression='gzip', index_col=0)
            df = df.rename({'contrast': 'target_contrast'}, axis=1)
            df['target_contrast_gene_name'] = df['target_contrast'].map(lambda x: gene_name_to_id.get(x, x))
            de_results_adatas.append(parse_DE_results_2_adata(df))

        combined_de_adata = anndata.concat(de_results_adatas, label='chunk')
        combined_de_adata.obs_names = combined_de_adata.obs_names.str.split('-').str[0]
        assert combined_de_adata.obs_names.is_unique

    # Annotate number of cells per target gene
    guide_cell_counts = obs_df[['perturbed_gene_id', 'guide_id', 'library_id', 'culture_condition']].value_counts().reset_index()
    n_cells_target_contrast = guide_cell_counts.groupby(['perturbed_gene_id', 'culture_condition'])['count'].sum().reset_index()
    n_cells_target_contrast.index = n_cells_target_contrast['perturbed_gene_id'] + '_' + n_cells_target_contrast['culture_condition']
    combined_de_adata.obs['n_cells_target'] = n_cells_target_contrast['count'].loc[combined_de_adata.obs_names]

    # Add gene names
    combined_de_adata.var = var_df.loc[combined_de_adata.var_names]

    # Save as anndata object
    combined_de_adata.write_h5ad(datadir + f'/DE_results/{experiment_name}.merged_DE_results.h5ad')

if __name__ == "__main__":
    main()