import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os

from preprocess import _convert_oak_path

def assign_sgrna(crispr_adata, min_sgrna_counts = 3, min_sgrna_counts_double = 10):
    """Assign sgRNAs to cells based on UMI count thresholds.
    
    This function assigns sgRNAs to cells using a two-step process:
    1. Identifies cells with a single sgRNA above min_sgrna_counts
    2. For cells with multiple sgRNAs, assigns the dominant sgRNA if:
       - The second highest sgRNA has < min_sgrna_counts_double UMIs
       - The highest sgRNA has > median UMIs of single sgRNA cells
    
    Params:
        crispr_adata: AnnData object containing sgRNA UMI counts
        min_sgrna_counts: Minimum UMI threshold for initial sgRNA detection (default: 3)
        min_sgrna_counts_double: UMI threshold for second highest sgRNA in multi-sgRNA cells (default: 10)
        
    Returns:
        None - modifies crispr_adata.obs in place by adding:
            - guide_id: Assigned sgRNA ID for each cell (NaN if unassigned)
            - top_guide_umi_counts: UMI count of highest sgRNA for each cell
    """
    # Exclude blacklisted sgRNAs
    crispr_adata.var['exclude_sgrna'] = crispr_adata.var['inefficient'] | crispr_adata.var['nonspecific'] | (crispr_adata.var['sgrna_type'] == 'ProbeNTC')
    # crispr_adata = crispr_adata[:, ~crispr_adata.var['exclude_sgrna']]

    # Count sgRNAs at UMI threshold t
    sgrna_assignment_mat = crispr_adata.X.copy()
    sgrna_assignment_mat[:, crispr_adata.var['exclude_sgrna']] = 0
    sgrna_assignment_mat[sgrna_assignment_mat < min_sgrna_counts] = 0

    sgrna_assignment_bin = sgrna_assignment_mat.copy()
    sgrna_assignment_bin[sgrna_assignment_bin > 0] = 1
    
    # Convert sparse matrices to dataframes
    sgrna_assignment_bin = pd.DataFrame(sgrna_assignment_bin.toarray(), 
                                      index=crispr_adata.obs_names,
                                      columns=crispr_adata.var_names)
    sgrna_assignment_mat = pd.DataFrame(sgrna_assignment_mat.toarray(),
                                      index=crispr_adata.obs_names, 
                                      columns=crispr_adata.var_names)
    
    # Store as layers
    crispr_adata.layers['binary_assignment'] = sgrna_assignment_bin.values
    crispr_adata.layers['umi_assignment'] = sgrna_assignment_mat.values
    
    # Define cells with single sgRNA
    single_sgrna_cells = sgrna_assignment_bin.index[sgrna_assignment_bin.sum(1) == 1].tolist()
    multi_sgrna_cells = sgrna_assignment_bin.index[sgrna_assignment_bin.sum(1) >= 2].tolist()
    sgrna_UMI_median = sgrna_assignment_mat.loc[single_sgrna_cells].max().median()
    top2_sgrnas = sgrna_assignment_mat.loc[multi_sgrna_cells].T.apply(lambda x: pd.Series(sorted(x, reverse=True)[:2])).T
    single_sgrna_cells.extend(top2_sgrnas[ (top2_sgrnas.min(axis=1) < min_sgrna_counts_double) & (top2_sgrnas.max(axis=1) > sgrna_UMI_median) ].index.tolist())
    
    # Assing top sgRNA to cell with unique target
    assigned_sgrna = sgrna_assignment_mat.loc[single_sgrna_cells].idxmax(axis=1)
    max_sgrna_umi = sgrna_assignment_mat.max(axis=1)

    crispr_adata.obs['guide_id'] = np.nan
    crispr_adata.obs.loc[assigned_sgrna.index, 'guide_id'] = assigned_sgrna
    crispr_adata.obs.loc[multi_sgrna_cells, 'guide_id'] = np.where(crispr_adata.obs.loc[multi_sgrna_cells, 'guide_id'].isna(), 'multi_sgRNA', crispr_adata.obs.loc[multi_sgrna_cells, 'guide_id']) 
    crispr_adata.obs['top_guide_umi_counts'] = np.nan
    crispr_adata.obs.loc[max_sgrna_umi.index, 'top_guide_umi_counts'] = max_sgrna_umi

def plot_sgrna_assignment(crispr_adata, min_sgrna_counts = 3, figsize=(15,5)):
    """Plot diagnostic figures for sgRNA assignment.
    
    Args:
        crispr_adata: AnnData object containing sgRNA data
        min_sgrna_counts: UMI threshold used for initial assignment
        
    Returns:
        matplotlib figure with 3 subplots
    """
    binary_assignment = crispr_adata.layers['binary_assignment']
    umi_assignment = crispr_adata.layers['umi_assignment']
    
    fig, axs = plt.subplots(1, 3, figsize=figsize)
    
    # Plot number of cells with # sgRNAs
    n_sgrnas_cells = pd.Series(binary_assignment.sum(axis=1)).value_counts()
    n_sgrnas_cells.index = n_sgrnas_cells.index.astype(int)
    sns.barplot(n_sgrnas_cells, ax=axs[0])
    axs[0].set_xlabel(f'# sgRNAs over threshold (>= {min_sgrna_counts} UMIs)')
    axs[0].set_ylabel('# cells')

    # Plot UMI counts for cells with one sgRNA
    single_mask = binary_assignment.sum(axis=1) == 1
    single_cell_umis = umi_assignment[single_mask].max(axis=1)
    values = np.sort(single_cell_umis)
    sgrna_UMI_median = np.median(values)
    axs[1].plot(values, '.')
    axs[1].set_yscale('log')
    axs[1].axhline(y=sgrna_UMI_median, color='r', linestyle='--', label='median')
    axs[1].set_xlabel('single sgRNA cell rank')
    axs[1].set_ylabel('sgRNA UMI counts')
    axs[1].set_title('UMI counts for cells with 1 sgRNA over threshold')
    
    # Plot UMI counts for cells with multiple sgRNAs
    multi_mask = binary_assignment.sum(axis=1) >= 2
    multi_cell_umis = umi_assignment[multi_mask]
    multi_cell_ranks = np.argsort(multi_cell_umis.max(axis=1))
    
    # Gather UMI counts and ranks for plotting
    plot_umis = []
    plot_ranks = []
    for i, cell_idx in enumerate(multi_cell_ranks):
        cell_umis = multi_cell_umis[cell_idx]
        nonzero_umis = cell_umis[cell_umis > 0]
        plot_umis.extend(nonzero_umis)
        plot_ranks.extend([i] * len(nonzero_umis))

    axs[2].hist2d(
        plot_ranks,
        np.log10(plot_umis),
        bins=100,
        norm=matplotlib.colors.LogNorm()
    )
    axs[2].axhline(y=np.log10(sgrna_UMI_median), color='r', linestyle='--', label='median')
    axs[2].set_xlabel('multi sgRNA cell rank')
    axs[2].set_ylabel('sgRNA UMI counts')
    axs[2].set_title('UMI counts for cells with >1 sgRNA over threshold')
    
    return fig

if __name__ == "__main__":
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser(description='Assign sgRNAs for GWT experiment')
    parser.add_argument('experiment', type=str,
                       help='Experiment ID to process. If not specified, processes all experiments')
    parser.add_argument('--config', type=str,
                       default='/oak/stanford/groups/pritch/users/emma/bin/GWT_perturbseq_analysis/metadata/experiment_config.yaml',
                       help='Path to experiment config YAML file')
    parser.add_argument('--plot_dir', type=str, default='../../results/',
                       help='Path to store plots')
    args = parser.parse_args()

    # Load config file
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    experiment = args.experiment
    config = config[experiment]    

    datadir = _convert_oak_path(config['datadir'])
    sample_metadata_csv = _convert_oak_path(config['sample_metadata'])
    sample_metadata = pd.read_csv(sample_metadata_csv, index_col=0)

    results_dir = f'{args.plot_dir}/{experiment}/'
    os.makedirs(results_dir, exist_ok=True)

    assignment_merged = pd.DataFrame()
    for sample_id in sample_metadata['sample_id']:
        sgrna_h5ad = f"{datadir}{sample_id}.sgRNA.h5ad"
        crispr_a = sc.read_h5ad(sgrna_h5ad)
        assign_sgrna(crispr_a)
        fig = plot_sgrna_assignment(crispr_a)
        fig.savefig(f'{results_dir}/{experiment}_{sample_id}.sgRNA_assignment.pdf')
        fig.savefig(f'{results_dir}/{experiment}_{sample_id}.sgRNA_assignment.png')
        assignment_s = crispr_a.obs[['guide_id', 'top_guide_umi_counts']]
        assignment_merged = pd.concat([assignment_merged, assignment_s])

    assignment_merged.to_csv(f'{datadir}/sgrna_assignment.csv')