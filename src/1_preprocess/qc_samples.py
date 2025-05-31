import yaml
import anndata as ad
import pandas as pd
import numpy as np

from preprocess import _convert_oak_path
from qc_plots import get_qc_summary
from sgrna_assignment import sgrna_assignments2adata
import argparse

import scipy
import dask.array as da

def add_target_expression(adata, perturbed_gene_col='perturbed_gene_id', 
                           var_names='gene_ids', ntc_col='perturbed_gene_name', 
                           ntc_label='NTC'):
    """
    Calculate target gene expression statistics per guide for a given sample.
    
    Parameters
    ----------
    adata : str
        Anndata object (expects log-normalized counts in X)
    perturbed_gene_col : str, optional (default: 'perturbed_gene_id')
        Column name containing perturbed gene IDs
    var_names : str, optional (default: 'gene_ids')
        Column name containing gene IDs in var
    ntc_col : str, optional (default: 'perturbed_gene_name')
        Column name containing perturbed gene names
    ntc_label : str, optional (default: 'NTC')
        Label for non-targeting control guides
        
    Returns
    -------
    anndata 
    """
    # Get all perturbed genes that are also measured
    all_perturbed_gene_ids = adata.obs[perturbed_gene_col].dropna().unique()
    measured_perturbed_gene_ids = np.intersect1d(all_perturbed_gene_ids, adata.var[var_names])
    adata.obs[perturbed_gene_col + '_test'] = adata.obs[perturbed_gene_col].copy()
    adata.obs.loc[~adata.obs[perturbed_gene_col + '_test'].isin(measured_perturbed_gene_ids), 
                 perturbed_gene_col + '_test'] = np.nan
    perturbed_gene_col = perturbed_gene_col + '_test'

    # Initialize binary matrix (cells × perturbed genes)
    cell_gene_matrix = scipy.sparse.lil_matrix((adata.n_obs, adata.n_vars), dtype=bool)
    # Create mapping from gene ID to column index
    gene_to_idx = {gene: i for i, gene in enumerate(adata.var[var_names])}

    # Fill the matrix: True if cell has that gene perturbed
    for cell_idx, gene in enumerate(adata.obs[perturbed_gene_col]):
        if pd.notna(gene):
            cell_gene_matrix[cell_idx, gene_to_idx[gene]] = True

    # Get NTC cells mask
    ntc_mask = adata.obs[ntc_col] == ntc_label

    # Sum of expression in NTCs
    adata_ntc = adata[ntc_mask].copy()
    ntc_mat = adata_ntc.X.map_blocks(lambda x: x.toarray(), 
                                            dtype=adata_ntc.X.dtype, 
                                            meta=np.array([]))
    # Calculate sum and sum of squares for NTC cells
    ntc_sums = ntc_mat.sum(0).compute()
    ntc_mat_squared = ntc_mat.map_blocks(lambda x: x ** 2, 
                                        dtype=ntc_mat.dtype,
                                        meta=np.array([]))
    ntc_sums_sq = ntc_mat_squared.sum(0).compute()
    n_cells_ntc = adata_ntc.n_obs

    # Get expression of perturbed gene for each cell
    cell_gene_matrix = cell_gene_matrix.tocsr()
    rows, cols = cell_gene_matrix.nonzero()
    target_mat = adata.X[rows, :][:, cols]
    target_mat_square = target_mat.map_blocks(lambda x: x.toarray(), 
                                            dtype=target_mat.dtype, 
                                            meta=np.array([]))
    target_gene_values = da.diag(target_mat_square).compute()

    # Save results in adata
    adata.var['sum_lognorm_expr_ntc'] = ntc_sums
    adata.var['sum_sq_lognorm_expr_ntc'] = ntc_sums_sq

    adata.uns['n_ntc'] = n_cells_ntc

    adata.obs['target_lognorm_expr'] = np.nan
    adata.obs.loc[adata.obs_names[rows], 'target_lognorm_expr'] = target_gene_values
    adata.obs['target_lognorm_expr'] = np.where(adata.obs['target_lognorm_expr'].isna() & (~ntc_mask), 0, adata.obs['target_lognorm_expr'])
    adata.obs.drop(perturbed_gene_col, axis=1, inplace=True)
    return(adata)

def parse_sample(f, datadir, sgrna_library_metadata, sample_id):
    adata = ad.experimental.read_lazy(f)
    if not isinstance(adata.obs, pd.DataFrame):
        adata.obs = adata.obs.to_dataframe()
    sgrna_assignments2adata(adata, datadir, sgrna_library_metadata=sgrna_library_metadata, sample_id=sample_id)
    qc_summary = get_qc_summary(adata)
    # Exclude low quality transcriptomes & perturbations
    exclude_filt = adata.obs['low_quality'] | (adata.obs['guide_id'] == 'multi_sgRNA') | (adata.obs['guide_id'].isna())
    adata = adata[~exclude_filt]
    # Count number of cells per perturbation/guide
    count_perturbs = pd.DataFrame({
        'guide_id': adata.obs['guide_id'].value_counts().index,
        'n_cells': adata.obs['guide_id'].value_counts().values,
        'sample_id': sample_id
    })
    return adata, qc_summary, count_perturbs
    
def main():
    parser = argparse.ArgumentParser(description='Process single cell RNA-seq data for a specific sample')
    parser.add_argument('--experiment_name', type=str, required=True, help='Name of the experiment')
    parser.add_argument('--sample_id', type=str, required=True, help='Sample ID to process')
    parser.add_argument('--config', type=str, default='../../metadata/experiments_config.yaml',
                       help='Path to the experiment configuration file')
    args = parser.parse_args()

    # Read config
    config_file = args.config
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    config = config[args.experiment_name]
    datadir = _convert_oak_path(config['datadir'])
    sample_metadata_csv = _convert_oak_path(config['sample_metadata'])

    sample_metadata = pd.read_csv(sample_metadata_csv, index_col=0)
    sgrna_library_metadata = pd.read_csv('../../metadata/sgRNA_library_curated.csv', index_col=0)

    # Process single sample
    f = f'{datadir}/tmp/{args.sample_id}.scRNA.h5ad'
    adata, qc_summary, count_perturbs = parse_sample(f, datadir, sgrna_library_metadata, args.sample_id)
    
    if qc_summary is not None:
        # Save QC summary
        qc_summary.to_csv(f'{datadir}/tmp/{args.sample_id}.qc_summary.csv')
        
        # Save perturbation counts
        count_perturbs.to_csv(f'{datadir}/tmp/{args.sample_id}.perturbation_counts.csv')
        
        # Save processed AnnData for merge
        if not isinstance(adata.var, pd.DataFrame):
            adata.var = adata.var.to_dataframe()
        if not isinstance(adata.obs, pd.DataFrame):
            adata.obs = adata.obs.to_dataframe()
        
        # Store log-normalized expression of target gene
        adata = add_target_expression(adata)

        adata.X = adata.layers['counts'].copy()
        del adata.layers['counts']
        output_file = f'{datadir}/tmp/{args.sample_id}.scRNA.postQC.h5ad'
        adata[:, adata.var['total_counts'] > 0].write_h5ad(output_file)

if __name__ == '__main__':
    main()