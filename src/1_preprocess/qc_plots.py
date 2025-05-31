from typing import Dict, List, Optional
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import scipy

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

TCELL_MARKERS_URL = 'https://docs.google.com/spreadsheets/d/1NSjrMofwd7etOJVvqHXFMNV0--gjzg8BzHMWGxK7SUI'

def load_Tcell_markers(adata: anndata.AnnData) -> Dict[str, List[str]]:
    Tcell_markers = pd.read_csv(f'{TCELL_MARKERS_URL}/export?format=csv')
    Tcell_markers = Tcell_markers[Tcell_markers['gene_name'].isin(adata.var['gene_name'])]
    Tcell_markers['state_program_name'] = Tcell_markers['state_program_name'].astype('category').cat.reorder_categories(Tcell_markers.state_program_name.unique())
    markers_dict = Tcell_markers.groupby('state_program_name')['gene_name'].apply(list).to_dict()
    return markers_dict

def plot_ncells_sample(adata: anndata.AnnData, color_by: str = 'culture_condition', sample_col: str = 'cell_sample_id') -> None:
    cols = [color_by, sample_col]
    pl_df = sc.get.obs_df(adata, cols)
    pl_df = pl_df.groupby(cols).size().reset_index(name='count')
    pl_df = pl_df[pl_df['count'] > 0].copy()

    sns.barplot(data=pl_df, x=sample_col, y='count', hue=color_by)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, title=color_by)
    plt.xticks(rotation=90)
    plt.ylabel('# cells (pre-QC)')

def plot_markers_sample_dotplot(adata: anndata.AnnData, group_by: str = 'culture_condition', sample_col: str = 'cell_sample_id', **kwargs) -> None:
    markers_dict = load_Tcell_markers(adata)
    order_samples = adata.obs.sort_values([group_by])[sample_col].astype(str).unique()
    adata.obs[sample_col] = adata.obs[sample_col].cat.reorder_categories(order_samples)
    sc.pl.dotplot(adata, markers_dict, groupby=sample_col, **kwargs)

def plot_markers_umap(adata: anndata.AnnData, markers_dict: Optional[Dict[str, List[str]]] = None, markers_keys: Optional[List[str]] = None, **kwargs) -> None:
    if markers_dict is None:
        markers_dict = load_Tcell_markers(adata)
    
    if markers_keys is not None:
        markers_dict = {k: markers_dict[k] for k in markers_keys if k in markers_dict}

    pl_genes = [g for k,v in markers_dict.items() for g in v]
    pl_titles = [f'{g} ({k})' for k,v in markers_dict.items() for g in v]

    sc.pl.umap(adata, color=pl_genes, title=pl_titles, **kwargs)


def calculate_perturbed_gene_expression(adata: anndata.AnnData, perturbed_gene_col: str = 'perturbed_gene_id', var_names: str = 'gene_ids',
    ntc_col: str = 'perturbed_gene_name', ntc_label: str = 'NTC') -> pd.DataFrame:
    """Calculate expression of perturbed genes in cells with perturbations.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with perturbation information in .obs
    perturbed_gene_col : str, optional (default='perturbed_gene_id')
        Column in adata.obs containing perturbed gene IDs
    var_names : str, optional (default='gene_ids') 
        Column in adata.var containing gene IDs that match perturbed_gene_col
        
    Returns
    -------
    Modifies AnnData in place, and returns pd.DataFrame
        DataFrame containing perturbed gene expression info:
        - perturbed_gene: gene ID that was perturbed
        - perturbed_gene_expr: expression in perturbed cells
        - perturbed_gene_mean_ntc: mean expression in NTC cells
        - perturbed_gene_std_ntc: std dev of expression in NTC cells
    """

    # Get all perturbed genes that are also measured
    all_perturbed_gene_ids = adata.obs[perturbed_gene_col].dropna().unique()
    measured_perturbed_gene_ids = np.intersect1d(all_perturbed_gene_ids, adata.var[var_names])
    adata.obs[perturbed_gene_col + '_test'] = adata.obs[perturbed_gene_col].copy()
    adata.obs.loc[~adata.obs[perturbed_gene_col + '_test'].isin(measured_perturbed_gene_ids), perturbed_gene_col + '_test'] = np.nan
    perturbed_gene_col = perturbed_gene_col + '_test'

    # Initialize binary matrix (cells × perturbed genes)
    n_cells = adata.n_obs
    cell_gene_matrix = scipy.sparse.lil_matrix((adata.n_obs, adata.n_vars), dtype=bool)
    # Create mapping from gene ID to column index
    gene_to_idx = {gene: i for i, gene in enumerate(adata.var[var_names])}

    # Fill the matrix: True if cell has that gene perturbed
    for cell_idx, gene in enumerate(adata.obs[perturbed_gene_col]):
        if pd.notna(gene):
            cell_gene_matrix[cell_idx, gene_to_idx[gene]] = True

    # Get NTC cells mask
    ntc_mask = adata.obs[ntc_col] == ntc_label

    # Convert to more efficient format for operations
    cell_gene_matrix = cell_gene_matrix.tocsr()
    rows, cols = cell_gene_matrix.nonzero()
    result = np.full(adata.X.shape[0], np.nan)
    values = np.array(adata.X[rows, cols]).flatten()
    perturbed_gene_expr = pd.Series(values, index=adata.obs_names[rows])
    
    # Get summary of expression in NTCs
    ntc_values = adata.X[np.where(ntc_mask)[0],:][:, cols]
    perturbed_gene_mean_ntc = ntc_values.mean(axis=0)
    perturbed_gene_mean_ntc = np.array(perturbed_gene_mean_ntc).flatten()
    perturbed_gene_mean_ntc = pd.Series(perturbed_gene_mean_ntc, index=adata.obs_names[rows])
    perturbed_gene_std_ntc = np.sqrt(ntc_values.power(2).mean(axis=0).A1 - np.power(ntc_values.mean(axis=0).A1, 2))
    perturbed_gene_std_ntc = pd.Series(perturbed_gene_std_ntc, index=adata.obs_names[rows])

    # Add to adata
    adata.obs['perturbed_gene_expr'] = perturbed_gene_expr
    adata.obs['perturbed_gene_mean_ntc'] = perturbed_gene_mean_ntc
    adata.obs['perturbed_gene_std_ntc'] = perturbed_gene_std_ntc

    perturbed_gene_expr_df = pd.DataFrame({
        'perturbed_gene':adata.obs[perturbed_gene_col],
        'perturbed_gene_expr':perturbed_gene_expr, 
        'perturbed_gene_mean_ntc':perturbed_gene_mean_ntc,
        'perturbed_gene_std_ntc':perturbed_gene_std_ntc,
        })
    perturbed_gene_expr_df['n_ntc_cells'] = ntc_values.shape[0]
    return(perturbed_gene_expr_df)


def test_knockdown_simple(perturbed_gene_expr_df, group_col='perturbed_gene'):
    """
    Test for significant knockdown of perturbed genes using Welch's t-test.
    
    Args:
        perturbed_gene_expr_df: DataFrame with perturbed gene expression data (output of `calculate_perturbed_gene_expression`)
        group_col: Column name to use for grouping (default: 'perturbed_gene')
        
    Returns:
        DataFrame with test statistics and knockdown significance results
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    mean_perturbed_gene_expr_df = perturbed_gene_expr_df.groupby(group_col).mean()
    mean_perturbed_gene_expr_df['perturbed_gene_expr_std'] = perturbed_gene_expr_df.groupby(group_col)['perturbed_gene_expr'].std()
    mean_perturbed_gene_expr_df['n_perturbed_cells'] = perturbed_gene_expr_df.groupby(group_col).size()
    mean_perturbed_gene_expr_df['n_ntc_cells'] = perturbed_gene_expr_df['n_ntc_cells'].mean()
    mean_perturbed_gene_expr_df = mean_perturbed_gene_expr_df.reset_index()

    # Calculate Welch's t-test for each row
    t_stats = []
    p_values = []
    for _, row in mean_perturbed_gene_expr_df.iterrows():
        t_stat, p_value = stats.ttest_ind_from_stats(
            row['perturbed_gene_expr'], row['perturbed_gene_expr_std'],row['n_perturbed_cells'],
            row['perturbed_gene_mean_ntc'], row['perturbed_gene_std_ntc'], row['n_ntc_cells'],
            equal_var=False)
        t_stats.append(t_stat)
        p_values.append(p_value)

    mean_perturbed_gene_expr_df['t_statistic'] = t_stats
    mean_perturbed_gene_expr_df['p_value'] = p_values

    # Compute Benjamini-Hochberg adjusted p-values
    valid_pvals = mean_perturbed_gene_expr_df['p_value'].dropna()
    _, adj_pvals, _, _ = multipletests(valid_pvals, method='fdr_bh')
    mean_perturbed_gene_expr_df.loc[valid_pvals.index, 'adj_p_value'] = adj_pvals
    mean_perturbed_gene_expr_df['signif_knockdown'] = (mean_perturbed_gene_expr_df['adj_p_value'] < 0.1) & (mean_perturbed_gene_expr_df['t_statistic'] < 0)
    
    return mean_perturbed_gene_expr_df

def get_qc_summary(adata):
    '''
    Get QC metrics per adata object with sgRNA assignment 
    '''
    mean_qc_stats = adata.obs.groupby(['library_id', 'lane_id'])[['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'top_guide_UMI_counts']].mean().reset_index()
    mean_qc_stats = mean_qc_stats.rename(columns={
        'total_counts': 'mean_total_counts',
        'n_genes_by_counts': 'mean_n_genes', 
        'pct_counts_mt': 'mean_pct_counts_mt',
        'top_guide_UMI_counts': 'mean_top_guide_UMI_counts'
    })


    # Low quality mask
    adata.obs['low_quality'] = (adata.obs['pct_counts_mt'] > 5) | (adata.obs['n_genes'] < 500)

    mean_qc_stats['n_cells'] = adata.obs.groupby(['library_id', 'lane_id']).size().reset_index()[0]
    mean_qc_stats['n_low_quality_cells'] = adata.obs[adata.obs['low_quality']].groupby(['library_id', 'lane_id']).size().reset_index()[0]

    # Count groups of cells by guide assignment
    pl_df = adata.obs[['guide_id','perturbed_gene_id', 'low_quality', 'library_id', 'lane_id']]
    pl_df['group'] = ''
    pl_df.loc[adata.obs['guide_id'].isna(), 'group'] = 'no sgRNA (>= 3 UMIs)'
    pl_df.loc[adata.obs['guide_id'].str.startswith('NTC', na=False), 'group'] = 'NTC single sgRNA'
    pl_df.loc[(~adata.obs['guide_id'].isna()) & 
            (~adata.obs['guide_id'].str.startswith('NTC', na=False)) & 
            (adata.obs['guide_id'] != 'multi_sgRNA'), 'group'] = 'targeting single sgRNA'
    pl_df.loc[adata.obs['guide_id'] == 'multi_sgRNA', 'group'] = 'multi sgRNA'
    pl_df = pl_df[~pl_df['low_quality']].copy()
    group_counts = pl_df.groupby(['library_id', 'lane_id', 'group']).size().unstack().reset_index()
    mean_qc_stats = pd.merge(mean_qc_stats, group_counts)


    mean_qc_stats['n_unique_guides'] = pl_df[pl_df['group'].isin(['targeting single sgRNA', 'NTC single sgRNA'])]['guide_id'].nunique()
    mean_qc_stats['n_unique_perturbed_genes'] = pl_df[pl_df['group'].isin(['targeting single sgRNA', 'NTC single sgRNA'])]['perturbed_gene_id'].nunique()

    mean_qc_stats['mean_cells_x_guide'] = pl_df[pl_df['group'].isin(['targeting single sgRNA', 'NTC single sgRNA'])].groupby(['guide_id']).size().mean()
    mean_qc_stats['mean_cells_x_perturbed_gene'] = pl_df[pl_df['group'].isin(['targeting single sgRNA', 'NTC single sgRNA'])].groupby(['perturbed_gene_id']).size().mean()
    return mean_qc_stats