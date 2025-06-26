'''
Utility functions
'''

import os
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
from typing import Optional, Union, List
import genomic_features as gf
import glob
from tqdm.notebook import tqdm

SPARSE_CHUNK_SIZE = 1_000_000

# Try to import rapids_singlecell, fall back to scanpy if not available
try:
    import rapids_singlecell as rsc
    HAS_RAPIDS = True
except ImportError:
    HAS_RAPIDS = False


def _convert_oak_path(path):
        """Helper function to convert oak paths between different mount points"""
        if not os.path.exists(path):
            return path.replace('/oak/stanford/groups/pritch/', '/mnt/oak/')
        return path

# ---- Feature selection ---- # 
def get_ribo_mito_genes():
    '''Get ribosomal and mitochondrial genes to exclude.'''
    ensdb = gf.ensembl.annotation(species="Hsapiens", version="108")
    genes = ensdb.genes()

    ribo_gene_ids = genes[genes.description.str.startswith('ribosomal protein') & (genes.gene_biotype == 'protein_coding')].gene_id.values
    mito_gene_ids = genes[genes.seq_name == 'MT'].gene_id.values
    mito_ribo_gene_ids = genes[genes.description.str.startswith('mitochondrial ribosomal protein') & (genes.gene_biotype == 'protein_coding')].gene_id.values
    return ribo_gene_ids, mito_gene_ids, mito_ribo_gene_ids

def feature_selection(
    adata: anndata.AnnData, 
    n_hvgs: int = 5000, 
    filter_ribo_mito: bool = True, 
    subset_adata: bool = True, 
    subset_obs: str = 'integration_sample_id',
    highx_min_mean_counts = 1000,
    highx_min_pct_dropouts_by_counts = 0.5,
    lowx_min_counts = 1,
    lowx_max_pct_dropouts_by_counts = 99.9,
    return_all = False,
    use_rapids = True
    ):
    '''Save table of selected features to use for integration.'''
    filter_genes = []
    if filter_ribo_mito:
        ribo_gene_ids, mito_gene_ids, mito_ribo_gene_ids = get_ribo_mito_genes()
        filter_genes.extend(ribo_gene_ids.tolist())
        filter_genes.extend(mito_gene_ids.tolist())
        filter_genes.extend(mito_ribo_gene_ids.tolist())

    if subset_adata:
        adata = sc.pp.sample(adata, fraction=0.1, copy=True)
    
    # Use rapids if available, otherwise use scanpy
    if HAS_RAPIDS and use_rapids:
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.calculate_qc_metrics(adata)
    else:
        sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Filter out highly and lowly expressed outlier genes
    highly_expressed_outiers = adata.var_names[(adata.var['mean_counts'] > highx_min_mean_counts) & (adata.var['pct_dropout_by_counts'] < highx_min_pct_dropouts_by_counts)].tolist()
    lowly_expressed_outliers = adata.var_names[(adata.var['total_counts'] < lowx_min_counts) & (adata.var['pct_dropout_by_counts'] > lowx_max_pct_dropouts_by_counts)].tolist()
    filter_genes.extend(highly_expressed_outiers)
    filter_genes.extend(lowly_expressed_outliers)

    # Subset
    adata = adata[:, ~adata.var_names.isin(filter_genes)].copy()

    # Get HVGs
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    if HAS_RAPIDS and use_rapids:
        rsc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor='seurat_v3')
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor='seurat_v3')
    
    if return_all:
        hvg_df = adata.var.copy()
    else:
        hvg_df = adata.var[adata.var['highly_variable']][['gene_name', 'gene_ids']]
    
    return(hvg_df)

# -- Utilities to read dataset -- #

def _make_obs_mask(obs_data, obs_filt):
    mask = True
    for col, value in obs_filt.items():
        condition = obs_data[col].isin(value) if isinstance(value, list) else (obs_data[col] == value)
        mask = mask & condition
    return mask

def load_cells(h5ad_files: List[str], obs_filt: dict):
    """Load cells from multiple h5ad files based on obs filter

    Parameters
    ----------
    h5ad_files : List[str]
        List of paths to h5ad files to load
    obs_filt : dict
        Dictionary specifying filtering criteria for observations/cells
        Keys are column names in obs, values are either single values or lists
        of values to keep

    Returns
    -------
    anndata.AnnData
        Combined AnnData object containing cells matching filter criteria
        from all input files, with var indices containing gene_ids and gene_name
    """
    a_ls = []
    for f in tqdm(h5ad_files, 'Reading files'):
        a = anndata.read_h5ad(f, backed='r')
        mask = _make_obs_mask(a.obs, obs_filt)
        a = a[mask].to_memory()
        if a.n_obs > 0:
            a_ls.append(a)
    adata = anndata.concat(a_ls, join='outer')
    adata.var = pd.concat([a.var[['gene_ids', 'gene_name']] for a in a_ls]).drop_duplicates().loc[adata.var_names]
    return(adata)
