import os,sys
import numpy as np
import anndata
import pandas as pd
import scanpy as sc
import glob
from tqdm import tqdm
import yaml

from copy import deepcopy
import argparse

sys.path.append(os.path.abspath('../'))
from utils import _convert_oak_path

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
    parser.add_argument('--config', type=str, required=True, help='Path to config YAML file')
    parser.add_argument('--force_combine', action='store_true', help='Force recombination of DE results even if merged file exists')
    args = parser.parse_args()

    # Load configuration from YAML file
    with open(args.config, 'r') as config_file:
        config = yaml.safe_load(config_file)

    # Extract parameters from config
    datadir = _convert_oak_path(config['datadir'])
    experiment_name = config['experiment_name']
    datadir = f'{datadir}/{experiment_name}/'
    run_name = config.get('run_name', 'default')
    force_combine = args.force_combine

    # Read cell-level metadata
    sgrna_library_metadata = pd.read_csv('../../metadata/sgRNA_library_curated.csv', index_col=0)
    gene_name_to_id = dict(zip(sgrna_library_metadata['perturbed_gene_id'], sgrna_library_metadata['perturbed_gene_name']))
    var_df = sc.read_h5ad(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad', backed=True).var.copy()

    de_results_dir = datadir + f'/DE_results_{run_name}/tmp/'

    file_exists = os.path.exists(datadir + f'/DE_results_{run_name}/{experiment_name}.merged_DE_results.h5ad')
    
    if file_exists and not force_combine:
        combined_de_adata = sc.read_h5ad(datadir + f'/DE_results_{run_name}/{experiment_name}.merged_DE_results.h5ad')
    else:
        # Read all csv.gz files from the DE results directory
        de_results_files = glob.glob(de_results_dir + 'DE_results.*.csv.gz')
        de_results_adatas = []

        for file in tqdm(de_results_files, desc="Processing DE result files"):
            try:
                df = pd.read_csv(file, compression='gzip', index_col=0)
                df = df.rename({'contrast': 'target_contrast'}, axis=1)
                df['target_contrast_gene_name'] = df['target_contrast'].map(lambda x: gene_name_to_id.get(x, x))
                de_results_adatas.append(parse_DE_results_2_adata(df))
            except EOFError:
                continue

        combined_de_adata = anndata.concat(de_results_adatas, label='chunk', join='outer')
        combined_de_adata.obs_names = combined_de_adata.obs_names.str.split('-').str[0]
        assert combined_de_adata.obs_names.is_unique

    # Annotate number of cells per target gene
    pbulk_adata_obs = anndata.experimental.read_lazy(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad').obs.to_dataframe()
    n_cells_target_contrast = pbulk_adata_obs[['perturbed_gene_id', 'culture_condition', 'n_cells']].copy()
    n_cells_target_contrast['test'] = n_cells_target_contrast['perturbed_gene_id'].astype(str) + '_' + n_cells_target_contrast['culture_condition'].astype(str)
    n_cells_target_contrast = n_cells_target_contrast.groupby('test')['n_cells'].sum().reset_index().set_index('test')
    combined_de_adata.obs['n_cells_target'] = n_cells_target_contrast['n_cells'].loc[combined_de_adata.obs_names]

    # Add gene names
    combined_de_adata.var = var_df.loc[combined_de_adata.var_names]

    # Add MASH results
    mash_results_files = glob.glob(de_results_dir + 'MASH_results*.csv.gz')
    if mash_results_files:  # Only process if MASH results exist
        mash_results = pd.DataFrame()
        for f in mash_results_files:
            mash_results_gr = pd.read_csv(f, compression='gzip', index_col=0)
            mash_results = pd.concat([mash_results, mash_results_gr])

        mash_results['obs_names'] = mash_results['target_contrast'] + '_' + mash_results['culture_condition']
        for stat in ['PosteriorMean', 'PosteriorSD', 'lfsr']:
            stat_layer = mash_results.pivot(index='obs_names', columns = 'var_names', values=stat)
            combined_de_adata.layers[f'MASH_{stat}'] = stat_layer.loc[combined_de_adata.obs_names, combined_de_adata.var_names].values

    # Save as anndata object
    combined_de_adata.write_h5ad(datadir + f'/DE_results_{run_name}/{experiment_name}.merged_DE_results.h5ad')

if __name__ == "__main__":
    main()