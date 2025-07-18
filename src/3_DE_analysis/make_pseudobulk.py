import os,sys
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import yaml
import argparse
from scipy import sparse
from tqdm import tqdm

# Add the parent directory to the path to import from sibling directory
sys.path.append(os.path.abspath('../'))
from utils import _convert_oak_path

def make_pseudobulk(h5ad_file, sample_metadata=None, condition_col='culture_condition', sgrna_col='guide_id'):
    """
    Create pseudobulk data from single-cell RNA-seq data.
    
    Parameters:
    -----------
    h5ad_file : str
        single-cell H5AD file path
    experiment_name : str
        Name of the experiment
    condition_col : str, optional
        Column name for condition information, default is 'culture_condition'
    sgrna_col : str, optional
        Column name for sgRNA information, default is 'guide_id'
    """
    adata = sc.read_h5ad(h5ad_file)

    if sample_metadata is not None:
        # Merge sample metadata
        merged = pd.merge(adata.obs, sample_metadata, how='left')
        merged.index = adata.obs_names.values
        adata.obs = merged.copy()

    sample_cols = sample_metadata.columns.tolist() + ['lane_id', 'guide_id', 'sequence', 'perturbed_gene_name','perturbed_gene_id', 'guide_type']
    adata.obs["sample_id"] = adata.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)

    n_cells_obs = adata.obs.value_counts(['sample_id'] + sample_cols).reset_index()
    n_cells_obs = n_cells_obs.set_index('sample_id').rename({'count':'n_cells'}, axis=1)
    pbulk_adata = sc.get.aggregate(adata, by=['sample_id'], func=['sum'])
    pbulk_adata.obs = n_cells_obs.loc[pbulk_adata.obs_names].copy()
    pbulk_adata.layers['sum'] = sparse.csr_matrix(pbulk_adata.layers['sum'])
    pbulk_adata.write_h5ad(h5ad_file.replace('.h5ad', '.DE_pseudobulk.h5ad'))
    return pbulk_adata


def merge_pseudobulks(h5ad_files, sample_cols = ['cell_sample_id', 'guide_id']):
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

    # Read and sum pseudobulks
    merged_adata = sc.read_h5ad(h5ad_files[0])
    merged_adata.obs["sample_id"] = merged_adata.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)
    n_cells_merged = merged_adata.obs[['sample_id', 'n_cells']].set_index('sample_id')

    for f in tqdm(h5ad_files[1:], desc="Merging pseudobulk files"):
        a = sc.read_h5ad(f)
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

def main_aggregate(args):
    # Remove the parser creation - args are already parsed
    sample_metadata = pd.read_csv(args.sample_metadata_csv, index_col=0)
    make_pseudobulk(
        args.h5ad_file,
        sample_metadata,
        args.condition_col,
        args.sgrna_col
    )

def main_merge(args):
    # Remove the parser creation - args are already parsed
    with open(args.DE_config, 'r') as config_file:  # Fixed: was args.config
        config = yaml.safe_load(config_file)
    
    # Extract parameters from config
    datadir = _convert_oak_path(config['datadir'])
    experiment_name = config['experiment_name']
    datadir = f'{datadir}/{experiment_name}'

    # Find all pseudobulk files
    pseudobulk_files = glob.glob('/mnt/oak/users/emma/data/GWT/CD4iR1_Psomagen/tmp/*DE_pseudobulk.h5ad') + glob.glob('/mnt/oak/users/emma/data/GWT/CD4iR2_Psomagen/tmp/*DE_pseudobulk.h5ad')
    print(f"Found {len(pseudobulk_files)} pseudobulk files")
    sample2files = {}
    for f in pseudobulk_files:
        sample_id = os.path.basename(f).split('.')[0]
        if sample_id not in sample2files:
            sample2files[sample_id] = []
        sample2files[sample_id].append(f)
    
    s = args.sample_id
    pbulk_adata = merge_pseudobulks(sample2files[s])
    pbulk_adata.write_h5ad(f'{datadir}/{experiment_name}_{s}.DE_pseudobulk.h5ad')

def main():
    parser = argparse.ArgumentParser(description='Pseudobulk processing script')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Aggregate subcommand with its specific arguments
    aggregate_parser = subparsers.add_parser('aggregate', help='Create pseudobulk data from single-cell RNA-seq data')
    aggregate_parser.add_argument('h5ad_file', type=str, help='Path to single cell h5ad file to pseudobulk')
    aggregate_parser.add_argument('--sample_metadata_csv', type=str, default=None, help='Path to sample metadata CSV file')
    aggregate_parser.add_argument('--condition_col', type=str, default='culture_condition', help='Column name for condition information')
    aggregate_parser.add_argument('--sgrna_col', type=str, default='guide_id', help='Column name for sgRNA information')
    
    # Merge subcommand with its specific arguments
    merge_parser = subparsers.add_parser('merge', help='Merge pseudobulk data')
    merge_parser.add_argument('sample_id', type=str, help='Sample ID to merge pseudobulk')
    merge_parser.add_argument('--DE_config', type=str, default='DE_config_full.yaml', help='Path to DE config')
    
    args = parser.parse_args()
    
    if args.command == 'aggregate':
        main_aggregate(args)
    elif args.command == 'merge':
        main_merge(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()