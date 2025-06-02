import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from scipy import sparse

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


def main():
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('h5ad_file', type=str, help='Path to single cell h5ad file to pseudobulk')
    parser.add_argument('--sample_metadata_csv', type=str, default=None, help='Path to sample metadata CSV file')
    parser.add_argument('--condition_col', type=str, default='culture_condition', help='Column name for condition information')
    parser.add_argument('--sgrna_col', type=str, default='guide_id', help='Column name for sgRNA information')
    args = parser.parse_args()
    
    sample_metadata = pd.read_csv(args.sample_metadata_csv, index_col=0)
    make_pseudobulk(
        args.h5ad_file,
        sample_metadata,
        args.condition_col,
        args.sgrna_col
    )


if __name__ == "__main__":
    main()