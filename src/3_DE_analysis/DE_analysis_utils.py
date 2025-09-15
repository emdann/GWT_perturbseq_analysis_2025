'''
Utility functions to explore differential expression analysis results (stored in AnnData format after merge_DE_results.py)
'''
from typing import List, Union, Optional
import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc
import glob
from tqdm import tqdm
import scipy

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

from copy import deepcopy


## --- GETTERS --- ## 

def get_DE_results_long(
    adata_de: anndata.AnnData, 
    targets: Union[List[str], None] = None, 
    genes: Union[List[str], None] = None, 
    effect_estimates: Union[str, List[str]] = ['log_fc', 'zscore'],
    signif_estimate: Union[None, str] = 'adj_p_value',
    signif_alpha = 0.1,
    gene_id_col: str = 'gene_name',
    target_id_col: str = 'target_contrast_gene_name',
    target_metadata_cols: List[str] = ['culture_condition', 'target_contrast']
    ):
    """
    Extract differential expression results from AnnData object in long table format.
    
    Parameters
    ----------
    adata_de : AnnData
        AnnData object containing DE results
    targets : List[str] or None
        List of target contrasts to include, or None for all
    genes : List[str] or None
        List of genes to include, or None for all
    effect_estimates : str or List[str]
        Names of layers containing effect estimates (e.g., log fold changes, z-scores)
    signif_estimate : str or None
        Name of layer containing significance values, or None if not needed
    signif_alpha : float
        Significance threshold for filtering
    gene_id_col : str
        Column in var containing gene identifiers
    target_id_col : str
        Column in obs containing target identifiers
    target_metadata_cols : List[str]
        Columns in obs to include as target metadata in the result
        
    Returns
    -------
    pd.DataFrame
        Long-format dataframe with DE results
    """
    # Convert effect_estimates to list if it's a string
    if isinstance(effect_estimates, str):
        effect_estimates = [effect_estimates]
    
    # Filter targets if specified
    if targets is not None:
        target_mask = adata_de.obs[target_id_col].isin(targets)
        adata_subset = adata_de[target_mask]
    else:
        adata_subset = adata_de
    
    # Filter genes if specified
    if genes is not None:
        gene_mask = adata_subset.var[gene_id_col].isin(genes)
        adata_subset = adata_subset[:, gene_mask]
    
    # Initialize results dataframe with proper structure
    results = None
    
    # Extract effect estimates
    for effect in effect_estimates:
        if effect not in adata_subset.layers:
            continue
            
        df = sc.get.obs_df(adata_subset, 
                          keys=adata_subset.var_names.tolist(),
                          layer=effect)
        
        # Melt to long format
        df_long = df.melt(ignore_index=False, 
                          var_name='gene', 
                          value_name=effect)
        
        # Add to results
        if results is None:
            results = df_long.reset_index()
        else:
            # Merge with existing results on index and gene
            results = results.merge(df_long.reset_index(), 
                                   on=['index', 'gene'], 
                                   how='outer')
    
    # If no valid effect estimates were found, create empty dataframe with required columns
    if results is None:
        results = pd.DataFrame(columns=['index', 'gene'] + effect_estimates)
    
    # Add significance values if requested
    if signif_estimate is not None and signif_estimate in adata_subset.layers:
        signif_df = sc.get.obs_df(adata_subset, 
                                 keys=adata_subset.var_names.tolist(),
                                 layer=signif_estimate)
        
        # Melt to long format
        signif_long = signif_df.melt(ignore_index=False, 
                                    var_name='gene', 
                                    value_name=signif_estimate)
        
        # Add to results
        results = results.merge(signif_long.reset_index()[['index', 'gene', signif_estimate]], 
                               on=['index', 'gene'], 
                               how='left')
        
        # Add significance indicator
        results['significant'] = results[signif_estimate] < signif_alpha
    
    # Add gene metadata
    gene_meta = adata_subset.var[[gene_id_col]].reset_index()
    results = results.merge(gene_meta, left_on='gene', right_on='index', suffixes=('', '_gene'))
    results = results.drop('index_gene', axis=1)
    
    # Add target metadata
    if target_id_col not in target_metadata_cols:
        target_metadata_cols.append(target_id_col)
    target_meta = adata_subset.obs[target_metadata_cols].reset_index()
    results = results.merge(target_meta, left_on='index', right_on='index')
    
    return results

def get_de_stats(adata_de, alpha=0.1, signif_col='MASH_lfsr', effect_col='MASH_PosteriorMean', axis='targets'):
    """
    Calculate differential expression stats by target or gene in the dataset.
    
    Parameters:
    -----------
    adata_de : AnnData
        AnnData object containing differential expression results
    alpha : float, default=0.1
        Significance threshold for the specified significance column
    signif_col : str, default='MASH_lfsr'
        Column name in adata_de.layers for significance values
    effect_col : str, default='MASH_PosteriorMean'
        Column name in adata_de.layers for effect size values
    axis:
        either targets or genes
    ontarget_df : DataFrame
        (optional) DataFrame containing on-target effect information
    
    
    Returns:
    --------
    DataFrame
        DataFrame with DE gene counts by target or gene
    """
    # Create significance matrix
    signif_mat = (adata_de.layers[signif_col] < alpha).astype(int)
    effect_sign = np.sign(adata_de.layers[effect_col])
    
    if axis == 'targets':
        # Count significant upregulated genes (significant and positive log_fc)
        n_up_genes = ((signif_mat == 1) & (effect_sign > 0)).sum(axis=1)
        
        # Count significant downregulated genes (significant and negative log_fc)
        n_down_genes = ((signif_mat == 1) & (effect_sign < 0)).sum(axis=1)
        
        # Total number of DE genes per sample
        n_de_genes = signif_mat.sum(axis=1)
        
        # Create a DataFrame to display the counts
        de_counts = pd.DataFrame({
            'target_contrast': adata_de.obs['target_contrast'],
            'target_name': adata_de.obs['target_contrast_gene_name'],
            'condition': adata_de.obs['culture_condition'],
            'n_cells_target': adata_de.obs['n_cells_target'],
            'n_up_genes': n_up_genes,
            'n_down_genes': n_down_genes,
            'n_total_de_genes': n_de_genes
        })
    elif axis == 'genes':
        # Count significant upregulated genes (significant and positive log_fc)
        n_up_targets = ((signif_mat == 1) & (effect_sign > 0)).sum(axis=0)
        
        # Count significant downregulated genes (significant and negative log_fc)
        n_down_targets = ((signif_mat == 1) & (effect_sign < 0)).sum(axis=0)
        
        # Total number of DE genes per sample
        n_de_targets = signif_mat.sum(axis=0)
        
        # Create a DataFrame to display the counts
        de_counts = pd.DataFrame({
            'gene_id': adata_de.var['gene_ids'],
            'gene_name': adata_de.var['gene_name'],
            'n_up_targets': n_up_targets,
            'n_down_targets': n_down_targets,
            'n_total_de_targets': n_de_targets
        })
    else:
        raise ValueError('axis is either genes or targets')
    return de_counts

def get_ontarget_effect(adata_de, signif_estimate='MASH_lfsr', signif_alpha=0.1):
    """
    Extract on-target effects from differential expression results.
    
    Parameters:
    -----------
    adata_de : AnnData
        AnnData object containing differential expression results
    signif_estimate : str, default='MASH_lfsr'
        Column name for significance estimates
    signif_alpha : float, default=0.1
        Significance threshold
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing on-target effects
    """
    all_targets = adata_de.obs['target_contrast'].astype(str).unique()
    measured_targets = adata_de.var_names[adata_de.var_names.isin(all_targets)].tolist()
    measured_targets_names = adata_de.var['gene_name'][adata_de.var_names.isin(all_targets)].tolist()

    ontarget_df = get_DE_results_long(
        adata_de, 
        targets=measured_targets, genes=measured_targets_names, 
        effect_estimates=['log_fc', 'zscore', 'baseMean'],
        signif_estimate=signif_estimate,
        signif_alpha=signif_alpha,
        target_id_col='target_contrast', 
        target_metadata_cols=['culture_condition']
    )

    ontarget_df = ontarget_df[ontarget_df['gene'] == ontarget_df['target_contrast']].copy()
    return(ontarget_df)


## --- PLOTTING UTILS --- ##

def plot_gene_expression_by_target(pbulk_adata, target_id, gene_id, condition='Stim8hr', 
                                  target_name_col='target', gene_name_col='gene_name', hue=None, ax=None):
    """
    Plot target vs NTC expression of a single gene with means highlighted. 
    
    Pass log-normalized counts in pbulk_adata.X to visualize normalized values.
    
    Parameters:
    -----------
    pbulk_adata : AnnData
        The AnnData object containing the pseudobulked gene expression data
    target_id : str
        The target ID to compare with NTC
    gene_id : str
        The gene ID to visualize expression for
    condition : str, default='Stim8hr'
        The condition to filter on
    target_name_col : str, default='target'
        Column in obs containing target identifiers
    gene_name_col : str, default='gene_name'
        Column in var containing gene identifiers
    hue : str, default=None
        Column in obs to use for grouping/coloring points. If None, no grouping is applied.
    ax : matplotlib.axes.Axes, default=None
        The axes on which to draw the plot. If None, a new figure and axes will be created.
    
    Returns:
    --------
    matplotlib.axes.Axes
        The plot axes
    """
    # Get gene name if gene_id is not directly in var_names
    if gene_id not in pbulk_adata.var_names:
        gene_var_names = pbulk_adata.var_names[pbulk_adata.var[gene_name_col] == gene_id].tolist()
        if not gene_var_names:
            raise ValueError(f"Gene {gene_id} not found in {gene_name_col}")
    else:
        gene_var_names = [gene_id]
    
    # Get data for plotting
    columns_to_get = ['culture_condition', target_name_col, gene_var_names[0]]
    if hue is not None:
        columns_to_get.append(hue)
    
    df = sc.get.obs_df(pbulk_adata[pbulk_adata.obs['keep_for_DE']], columns_to_get)
    
    # Filter for target of interest and NTC
    pl_df = df[df[target_name_col].isin(['NTC', target_id])]
    pl_df[target_name_col] = pl_df[target_name_col].astype(str)
    
    # Create figure if ax is not provided
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Increase font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12
    })
    
    # Create stripplot
    if hue is not None:
        sns.stripplot(data=pl_df[pl_df['culture_condition'] == condition], 
                     x=target_name_col, y=gene_var_names[0], hue=hue, dodge=True, ax=ax)
        
        # Calculate means for each target-hue pair
        means = pl_df[pl_df['culture_condition'] == condition].groupby(
            [target_name_col, hue])[gene_var_names[0]].mean().reset_index()
    else:
        sns.stripplot(data=pl_df[pl_df['culture_condition'] == condition], 
                     x=target_name_col, y=gene_var_names[0], ax=ax)
        
        # Calculate means for each target
        means = pl_df[pl_df['culture_condition'] == condition].groupby(
            [target_name_col])[gene_var_names[0]].mean().reset_index()

    means[gene_var_names[0]] = means[gene_var_names[0]].fillna(0)
    
    # Get the order of targets as they appear in the stripplot
    target_order = ax.get_xticklabels()
    target_order = [t.get_text() for t in target_order]
    
    # Plot mean for each group as a large dot
    if hue is not None:
        for _, row in means.iterrows():
            target = row[target_name_col]
            hue_val = row[hue]
            mean_value = row[gene_var_names[0]]
            target_idx = target_order.index(target)
            
            # Determine the x-position with dodge based on positions in sns.stripplot
            hue_levels = list(means[hue].unique())
            n_hue = len(hue_levels)
            hue_idx = hue_levels.index(hue_val)
            # Get dodge width from the stripplot (default is 0.8 in seaborn)
            dodge_width = 0.8 / n_hue if n_hue > 1 else 0
            # Center the means within each target group
            x_pos = target_idx - 0.4 + (hue_idx + 0.5) * dodge_width
            
            ax.scatter(x_pos, mean_value, 
                      s=150,  # larger size for mean dots
                      color='red',
                      edgecolor='black', 
                      alpha=0.5,
                      linewidth=1.5,
                      zorder=10)  # ensure means are plotted on top
    else:
        for _, row in means.iterrows():
            target = row[target_name_col]
            mean_value = row[gene_var_names[0]]
            target_idx = target_order.index(target)
            
            ax.scatter(target_idx, mean_value, 
                      s=150,  # larger size for mean dots
                      color='red',
                      edgecolor='black', 
                      alpha=0.5,
                      linewidth=1.5,
                      zorder=10)  # ensure means are plotted on top
    
    # Set title and y-axis label with larger font sizes
    ax.set_title(f"{gene_id} expression by {target_id} vs NTC in {condition}", fontsize=16)
    ax.set_ylabel(f"{gene_id} expression", fontsize=14)
    ax.set_xlabel(ax.get_xlabel(), fontsize=14)
    
    # If there's a legend, increase its font size
    if hue is not None and ax.get_legend() is not None:
        ax.legend(fontsize=12)
    
    return ax

def plot_effect_comparison(
    adata_de,
    comparison_params,
    n_top_genes=15,
    plot_correlation=False,
    axis_label='DE effect',
    corr_coords_xy=(0.01, 0.95),
    figsize=(8, 8),
    ax=None,
    annotate_genes=None,
    annotate_significant=False
):
    """
    Scatter plot comparing DE effects on all tested genes for a pair of conditions, targets or stats.

    Parameters:
    -----------
    adata_de : AnnData
        AnnData object containing DE analysis results
    comparison_params : dict
        Dictionary with parameters for comparison (put 2 values for parameter you want to compare on x or y axis)
        - 'target_contrast_gene_name': list of target genes to analyze
        - 'stat': list of statistics to use (e.g., 'MASH_PosteriorMean')
        - 'culture_condition': list of conditions to compare
    n_top_genes : int, default=15
        Number of top/bottom genes to annotate on each axis (ignored if annotate_genes is provided)
    plot_correlation : bool, default=False
        Whether to plot correlation statistics
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. If None, a new figure and axes will be created.
    annotate_genes : list or None, optional
        List of gene names to annotate on the plot. If None, will annotate top/bottom genes.
    annotate_significant : bool, default=False
        If True, color points based on significance in one or both conditions.

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure (None if ax was provided)
    ax : matplotlib.axes.Axes
        The axes containing the plot
    pl_df : pandas.DataFrame
        The dataframe used for plotting
    """
    # Validate parameters
    compare = [k for k, v in comparison_params.items() if len(v) == 2]
    if len(compare) > 1:
        raise ValueError(f"More than one parameter has length 2: {', '.join(compare)}")
    else:
        compare = compare[0]

    # Get DE results
    res_df = get_DE_results_long(
        adata_de,
        targets=comparison_params['target_contrast_gene_name'],
        effect_estimates=comparison_params['stat'],
        gene_id_col='gene_name',
        target_id_col='target_contrast_gene_name',
        target_metadata_cols=['culture_condition']
    )
    res_df = res_df[res_df['culture_condition'].isin(comparison_params['culture_condition'])]

    # Create pivot table
    if compare != 'stat':
        pl_df = res_df.pivot(columns=compare, index='gene', values=comparison_params['stat'][0])
        if annotate_significant:
            # Get significant genes for each condition
            sig_df = res_df.pivot(columns=compare, index='gene', values='significant')
    else:
        pl_df = res_df[comparison_params['stat'] + ['gene']].set_index('gene')

    # Create figure if ax is not provided
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Get column names for x and y axis labels
    x_col = pl_df.columns[0]  # First column (likely 'Rest')
    y_col = pl_df.columns[1]  # Second column (likely 'Stim8hr')

    if annotate_significant:
        # Create color mapping for significance
        sig_both = sig_df[sig_df[x_col] & sig_df[y_col]].index
        sig_x_only = sig_df[sig_df[x_col] & ~sig_df[y_col]].index
        sig_y_only = sig_df[~sig_df[x_col] & sig_df[y_col]].index
        
        # Plot points with different colors based on significance
        ax.scatter(pl_df.loc[sig_both, x_col], pl_df.loc[sig_both, y_col], 
                  s=4, color='purple', alpha=0.6, label='Significant in both')
        ax.scatter(pl_df.loc[sig_x_only, x_col], pl_df.loc[sig_x_only, y_col], 
                  s=4, color='darkred', alpha=0.6, label=f'Significant in {x_col}')
        ax.scatter(pl_df.loc[sig_y_only, x_col], pl_df.loc[sig_y_only, y_col], 
                  s=4, color='darkblue', alpha=0.6, label=f'Significant in {y_col}')
        
        # Plot non-significant points in gray
        non_sig = pl_df.index.difference(sig_both.union(sig_x_only).union(sig_y_only))
        ax.scatter(pl_df.loc[non_sig, x_col], pl_df.loc[non_sig, y_col], 
                  s=4, color='gray', alpha=0.3, label='Not significant')
    else:
        # Create scatter plot using column names from pl_df
        sns.scatterplot(data=pl_df, x=x_col, y=y_col, s=4, edgecolor='none', color='black', ax=ax)

    # Add dotted lines for x and y axes
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=1)
    ax.axvline(x=0, color='gray', linestyle=':', linewidth=1)

    if plot_correlation:
        # Calculate correlation between the two conditions, ignoring NAs
        corr, pval = scipy.stats.pearsonr(pl_df.dropna()[x_col], pl_df.dropna()[y_col])

        # Add correlation information as text
        ax.annotate(
            f'Correlation: {corr:.3f}\n(p{" < 10e-16" if pval < 1e-16 else f" = {pval:.3e}"})',
            xy=corr_coords_xy, xycoords='axes fraction', fontsize=12
        )

    # Add axis labels and title
    ax.set_xlabel(f'{axis_label} ({x_col})', fontsize=12)
    ax.set_ylabel(f'{axis_label} ({y_col})', fontsize=12)
    ax.set_title(f'{comparison_params["target_contrast_gene_name"][0]} knock-out effect\n{x_col} vs {y_col} comparison', fontsize=14)

    # Annotate genes
    texts = []
    if annotate_genes is not None:
        # Only annotate the specified genes if they are present in pl_df
        for gene in annotate_genes:
            if gene in pl_df.index:
                row = pl_df.loc[gene]
                texts.append(ax.text(row[x_col], row[y_col], gene, fontsize=8, color='black'))
    else:
        # Annotate genes with extreme values (top/bottom)
        top_genes_y = pl_df.nlargest(n_top_genes, y_col)
        bottom_genes_y = pl_df.nsmallest(n_top_genes, y_col)
        top_genes_x = pl_df.nlargest(n_top_genes, x_col)
        bottom_genes_x = pl_df.nsmallest(n_top_genes, x_col)
        extreme_genes = pd.concat([top_genes_y, bottom_genes_y, top_genes_x, bottom_genes_x]).drop_duplicates()
        for idx, row in extreme_genes.iterrows():
            texts.append(ax.text(row[x_col], row[y_col], idx, fontsize=8, color='black'))

    # Add white background to text labels
    for text in texts:
        text.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0))

    if annotate_significant:
        ax.legend(fontsize=10)

    if fig is not None:
        return fig, ax, pl_df
    else:
        return ax, pl_df

def plot_volcano(adata_de, target, condition, on='perturbation', significance_threshold=0.1, 
                 log_fc_threshold=1, top_n_genes=10, figsize=(10, 8),
                 effect_estimates=['log_fc'], target_metadata_cols=['culture_condition']):
    """
    Create a volcano plot from differential expression results for a specific target and condition.
    
    Parameters:
    -----------
    adata_de : AnnData
        AnnData object containing differential expression results
    target : str
        Target gene ID or name to plot results for
    condition : str
        Condition to filter results (e.g., 'Rest', 'Stim8hr')
    on: str
        'perturbation' if plotting DE effects of certain target contrast perturbation
        'downstream' if plotting regulator effects of certain downstream gene
    significance_threshold : float, optional
        P-value threshold for significance
    log_fc_threshold : float, optional
        Log fold change threshold to highlight
    top_n_genes : int, optional
        Number of top genes to label
    figsize : tuple, optional
        Figure size as (width, height)
    effect_estimates : list, optional
        List of effect estimate columns to include
    target_metadata_cols : list, optional
        List of target metadata columns to include
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure object
    pl_df : pandas.DataFrame
        The data used for plotting
    """
    # Get differential expression results for the target and condition
    if on=='perturbation':
        pl_df = get_DE_results_long(adata_de, targets=[target], 
                                   effect_estimates=effect_estimates, 
                                   target_metadata_cols=target_metadata_cols)
    elif on=='downstream':
        pl_df = get_DE_results_long(adata_de, genes=[target], 
                                   effect_estimates=effect_estimates, 
                                   target_metadata_cols=target_metadata_cols)        
    
    pl_df = pl_df[pl_df['culture_condition'] == condition]
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    
    # Calculate -log10(adj_p_value) for y-axis
    pl_df['-log10_adj_p_value'] = -np.log10(pl_df['adj_p_value'])
    
    # Set a threshold for significance
    pl_df['significant'] = pl_df['adj_p_value'] < significance_threshold
    
    # Create the scatter plot
    scatter = plt.scatter(
        pl_df['log_fc'], 
        pl_df['-log10_adj_p_value'],
        c=pl_df['significant'].map({True: 'red', False: 'black'}),
        alpha=0.7,
        s=15
    )
    
    # Add labels for the most significant genes
    top_genes = pl_df.nlargest(top_n_genes, '-log10_adj_p_value')
    texts = []
    for _, gene in top_genes.iterrows():
        if on=='perturbation':
            txt = plt.text(
                gene['log_fc'], 
                gene['-log10_adj_p_value'],
                gene['gene_name'],
                fontsize=8,
                ha='center',
                va='bottom',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
            )
        elif on=='downstream':
            txt = plt.text(
                gene['log_fc'], 
                gene['-log10_adj_p_value'],
                gene['target_contrast_gene_name'],
                fontsize=8,
                ha='center',
                va='bottom',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1)
            )
        texts.append(txt)

    adjustText.adjust_text(texts, arrowprops={"arrowstyle": "-", "color": "k", "zorder": 5})

    # Add horizontal line for significance threshold
    plt.axhline(y=-np.log10(significance_threshold), color='gray', linestyle='--', alpha=0.7)
    # Add vertical lines for log fold change thresholds
    plt.axvline(x=-log_fc_threshold, color='gray', linestyle='--', alpha=0.7)
    plt.axvline(x=log_fc_threshold, color='gray', linestyle='--', alpha=0.7)
    
    # Add labels and title
    plt.xlabel('Log2 Fold Change', fontsize=12)
    plt.ylabel('-log10(Adjusted p-value)', fontsize=12)
    plt.title(f'Volcano Plot of Differential Expression for {target} in {condition}', fontsize=14)
    
    # Add legend outside the plot
    plt.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f'Significant ({significance_threshold*100}% FDR)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label='Not Significant')
    ], loc='upper left', bbox_to_anchor=(1.05, 1))
    
    plt.tight_layout()
    return fig, pl_df