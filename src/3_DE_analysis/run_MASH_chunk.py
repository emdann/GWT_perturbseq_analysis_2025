import os,sys
import numpy as np
import anndata
import pandas as pd

from copy import deepcopy
import argparse
import yaml
import shutil
import pymash


def parse_chunk_4_mash(
    chunk_ix,
    de_results_tmp_dir,
    all_conditions=['Rest', 'Stim8hr'],
    state_col='culture_condition',
    target_col='contrast'
):
    """
    Parse a chunk of differential expression results for MASH analysis.
    
    Parameters:
    -----------
    chunk_ix : int
        Index of the chunk to parse
    de_results_tmp_dir : str
        Directory containing the DE results files
    all_conditions : list
        List of experimental conditions
    state_col : str
        Column name for the condition/state
    target_col : str
        Column name for the contrast/target
        
    Returns:
    --------
    tuple
        (Bhat, Shat) matrices for MASH analysis
    """
    results = pd.DataFrame()
    for cond in all_conditions:
        results_c = pd.read_csv(de_results_tmp_dir + f"DE_results.{cond}.chunk_{chunk_ix}.csv.gz", 
                               compression='gzip', index_col=0)
        results = pd.concat([results, results_c])

    results['test'] = results[target_col].astype(str) + '-' + results['variable'].astype(str)
    Bhat = results.pivot(columns=state_col, index='test', values='log_fc')
    Shat = results.pivot(columns=state_col, index='test', values='lfcSE')
    return(Bhat, Shat)

def save_mash_h5ad(mash, mash_file_h5ad):
    all_attrs = [ attr for attr in dir(mash) if not callable(getattr(mash, attr)) and not attr.startswith('__')]
    all_attrs.remove('adata')
    mash.adata.uns['mash_attributes'] = {a:getattr(mash, a) for a in all_attrs}
    mash.adata.write(mash_file_h5ad)

def _mash_group_recipe(Bhat, Shat, n_random_genes = 2000, random_genes_seed=42):
    np.random.seed(random_genes_seed)
    all_genes = Bhat.index.str.split('-').str[-1].unique()
    random_genes = np.random.choice(all_genes, size=n_random_genes, replace=False)

    mash = pymash.tl.Mash(Bhat.T, Shat.T)
    mash.impute_missing()
    strong_tests = mash.get_strong_subset(thresh=0.05)

    # Canonical covariances
    Ulist_canonical = pymash.covariance_utils.cov_canonical(mash)
    mash.extend_Ulist(Ulist_canonical, overwrite=False)

    # Subset to strong tests
    mash_strong = mash.copy()
    mash_strong.adata = mash_strong.adata[:, strong_tests].copy()

    # Get PC based covariances
    Ulist_pca = pymash.covariance_utils.cov_pca(mash_strong)
    Ulist_ed = pymash.covariance_utils.cov_ed(mash_strong, Ulist_pca)
    mash.extend_Ulist(Ulist_ed)

    random_tests = mash.adata.var_names[mash.adata.var_names.str.split('-').str[1].isin(random_genes)].tolist()
    mash.fit(subset_vars=random_tests, return_full=False)
    mash.predict()
    return(mash)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('--config', type=str, required=True, help='Path to config YAML file')
    parser.add_argument('--test_group', type=int, required=True, help='Chunk of targets to test')
    parser.add_argument('--n_chunks_per_group', type=int, default=5, help='Chunk of targets to test')
    args = parser.parse_args()

    # Load configuration from YAML file
    with open(args.config, 'r') as config_file:
        config = yaml.safe_load(config_file)

    # Extract parameters from config
    datadir = config['datadir']
    experiment_name = config['experiment_name']
    run_name = config.get('run_name', 'default')
    
    datadir = f'{datadir}/{experiment_name}'
    group_ix = args.test_group
    n_chunks_per_group = args.n_chunks_per_group

    # Create directories for results
    de_results_dir = f'{datadir}/DE_results_{run_name}/'
    de_results_tmp_dir = f'{de_results_dir}/tmp/'

    # Get targets to test
    target_chunk_matrix = pd.read_csv(f'{datadir}/DE_target2chunk.csv.gz', compression='gzip', index_col=0)
    n_chunks = target_chunk_matrix.shape[1]

    # Select chunks to process
    start_idx = group_ix * n_chunks_per_group
    end_idx = min(start_idx + n_chunks_per_group, n_chunks)
    chunk_ixs = range(start_idx, end_idx)

    # Initialize empty dataframes to store concatenated results
    Bhat_all = pd.DataFrame()
    Shat_all = pd.DataFrame()

    for i in chunk_ixs:
        Bhat_c, Shat_c = parse_chunk_4_mash(chunk_ix=i, de_results_tmp_dir=de_results_tmp_dir)
        Bhat_all = pd.concat([Bhat_all, Bhat_c], axis=0)
        Shat_all = pd.concat([Shat_all, Shat_c], axis=0)
    
    mash = _mash_group_recipe(Bhat_all, Shat_all)
    save_mash_h5ad(mash, f"{de_results_tmp_dir}/mash_results_group_{group_ix}.h5ad")

    res_df = mash.get_results()
    res_df[['target_contrast', 'var_names']] = res_df['test'].str.split('-', expand=True)
    res_df = res_df.rename({'index':'culture_condition'}, axis=1)
    res_df.to_csv(de_results_tmp_dir + f"MASH_results.group_{group_ix}.csv.gz", compression='gzip')


