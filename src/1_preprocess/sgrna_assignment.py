import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import time
import sys
import shutil
import multiprocessing as mp
import re

from preprocess import _convert_oak_path
from crispat.poisson_gauss import fit_PGMM

def process_batch(args):
    """
    Process a batch of gRNAs - redesigned to take a single argument for better pickle compatibility
    """
    gRNA_list, adata_crispr, output_dir, n_iter, start_idx, step, batch_id = args
    
    # Add debugging information
    print(f"[Worker {batch_id}] Starting batch with {min(step, len(gRNA_list) - start_idx)} gRNAs")
    sys.stdout.flush()  # Force output to be displayed immediately
    
    batch_perturbations = pd.DataFrame()
    batch_thresholds = pd.DataFrame()
    
    end_idx = min(start_idx + step, len(gRNA_list))
    for i in range(start_idx, end_idx):
        gRNA = gRNA_list[i]
        # Removed tqdm from inside the worker function
        try:
            if i % 5 == 0:  # Print progress every few gRNAs
                print(f"[Worker {batch_id}] Processing gRNA {i-start_idx+1}/{end_idx-start_idx}: {gRNA}")
                sys.stdout.flush()
                
            perturbed_cells, threshold, loss, map_estimates = fit_PGMM(
                gRNA, adata_crispr, output_dir, 2024, n_iter
            )
            if len(perturbed_cells) != 0:
                # get UMI_counts of assigned cells
                UMI_counts = adata_crispr[perturbed_cells, [gRNA]].X.toarray().reshape(-1)
                df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA, 'UMI_counts': UMI_counts})
                batch_perturbations = pd.concat([batch_perturbations, df], ignore_index=True)
                batch_thresholds = pd.concat([batch_thresholds, pd.DataFrame({'gRNA': [gRNA], 'threshold': [threshold]})])
        except Exception as e:
            print(f"[Worker {batch_id}] Error processing gRNA {gRNA}: {str(e)}")
            sys.stdout.flush()
    
    print(f"[Worker {batch_id}] Finished batch with {batch_perturbations.shape[0]} perturbations")
    sys.stdout.flush()
    
    return batch_perturbations, batch_thresholds

def assign_sgrna_crispat(adata_crispr, output_dir, start_idx=0, end_idx=None, UMI_threshold=3, n_iter=500, n_guides_parallel=4, num_cores=5):
    """
    Assign sgRNAs to cells using the CRISPAT Poisson-Gaussian Mixture Model with parallel processing.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Start timer to measure performance
    start_time = time.time()
    
    gRNA_list = adata_crispr.var_names.tolist()
    
    # Set up parallel processing
    if num_cores is None:
        num_cores = min(mp.cpu_count(), n_guides_parallel)
    else:
        num_cores = min(num_cores, mp.cpu_count(), len(gRNA_list))

    # Determine chunk boundaries
    if end_idx is None:
        end_idx = len(gRNA_list)
    else:
        end_idx = min(end_idx, len(gRNA_list))
    
    # Extract the subset of guides we'll process
    chunk_gRNAs = gRNA_list[start_idx:end_idx]
    chunk_adata = adata_crispr[:, chunk_gRNAs].copy()

    # Create batches for parallel processing
    batch_size = max(1, len(chunk_gRNAs) // num_cores)
    batch_indices = list(range(0, len(chunk_gRNAs), batch_size))
    
    # Prepare arguments for multiprocessing
    process_args = [(chunk_gRNAs, chunk_adata, output_dir, n_iter, start_batch, batch_size, idx) 
                   for idx, start_batch in enumerate(batch_indices)]
    
    # Process batches in parallel with better progress reporting
    print(f'Fitting Poisson-Gaussian Mixture Model for {len(gRNA_list)} gRNAs using {num_cores} cores')
    print(f'Each core will process approximately {batch_size} gRNAs')
    
    # Use 'fork' method instead of 'spawn' for Linux/macOS to avoid module import overhead
    ctx = mp.get_context('fork' if sys.platform != 'win32' else 'spawn')
    with ctx.Pool(processes=num_cores) as pool:
        results = []
        for i, res in enumerate(pool.imap(process_batch, process_args)):
            print(f"Completed batch {i+1}/{len(process_args)}")
            results.append(res)
    
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    print(f"All batches completed in {elapsed_time:.2f} seconds")
    
    # Combine results
    perturbations = pd.DataFrame()
    thresholds = pd.DataFrame()
    for batch_perturbations, batch_thresholds in results:
        perturbations = pd.concat([perturbations, batch_perturbations], ignore_index=True)
        thresholds = pd.concat([thresholds, batch_thresholds], ignore_index=True)

    # Filter by UMI threshold
    perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]

    # Make unique cell assignment - handle empty DataFrame case
    if len(perturbations) == 0:
        print("Warning: No perturbations passed the UMI threshold filter")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    assignment_size = perturbations.groupby('cell').size() 
    # Use drop=True to avoid duplicate 'cell' column
    assignment_crispat = perturbations.groupby('cell').apply(lambda x: x.loc[x['UMI_counts'].idxmax()]).reset_index(drop=True)
    # Add cell column back
    assignment_crispat['cell'] = perturbations.groupby('cell').apply(lambda x: x.name).values
    assignment_crispat['guide_id'] = np.where(assignment_size[assignment_crispat['cell']].values > 1, 
                                             'multi_sgRNA', 
                                             assignment_crispat['gRNA'])
    
    return assignment_crispat, perturbations, thresholds

def assign_sgrna_naive(crispr_adata, min_sgrna_counts = 3, min_sgrna_counts_double = 10):
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
    crispr_adata.var['exclude_sgrna'] = crispr_adata.var['nonspecific'] | (crispr_adata.var['sgrna_type'] == 'ProbeNTC')
    # crispr_adata = crispr_adata[:, ~crispr_adata.var['exclude_sgrna']]

    # Count sgRNAs at UMI threshold t
    sgrna_assignment_mat = crispr_adata.X.copy()
    sgrna_assignment_mat[:, crispr_adata.var['exclude_sgrna']] = 0
    if scipy.sparse.issparse(sgrna_assignment_mat):
        mask = sgrna_assignment_mat >= min_sgrna_counts
        sgrna_assignment_mat = sgrna_assignment_mat.multiply(mask)
    else:
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

def sgrna_assignments2adata(
    adata, 
    datadir, 
    sgrna_library_metadata=None,
    sample_id=None):
    '''Utility function to merge finalized assignments with anndata for scRNA-seq data.
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing scRNA-seq data
    datadir : str
        Directory containing sgRNA assignment files
    sgrna_library_metadata : pandas.DataFrame, optional
        DataFrame containing sgRNA library metadata. If None, loads from default location
    sample_id : str, optional
        If provided, only process sgRNA assignments for this sample
    '''

    if sgrna_library_metadata is None:
        sgrna_library_metadata = pd.read_csv('../../metadata/sgRNA_library_curated.csv', index_col=0)
    
    # Get list of assignment files to process
    if sample_id is not None:
        assignment_files = [f for f in os.listdir(datadir) if f.endswith('.sgrna_assignment.csv') and sample_id in f]
        if not assignment_files:
            raise ValueError(f"No sgRNA assignment files found for sample {sample_id}")
    else:
        assignment_files = [f for f in os.listdir(datadir) if f.endswith('.sgrna_assignment.csv')]
    
    # Concatenate sgrna_assignment files
    sgrna_assignment = pd.concat([pd.read_csv(os.path.join(datadir, f)) for f in assignment_files])
    assert sgrna_assignment.cell.is_unique

    # Get guide metadata
    sgrna_assignment = sgrna_assignment[~sgrna_assignment['guide_id'].str.startswith('ProbeNTC-')].copy() # Exclude probe NTCs
    sgrna_assignment = pd.merge(sgrna_assignment, sgrna_library_metadata.rename({"sgrna_id":'guide_id'}, axis=1), how='left')
    sgrna_assignment['perturbed_gene_id'] = np.where(sgrna_assignment['perturbed_gene_id'].isna(), sgrna_assignment['perturbed_gene_name'], sgrna_assignment['perturbed_gene_id'])
    sgrna_assignment = sgrna_assignment.drop('gRNA', axis=1).rename({"UMI_counts":'top_guide_UMI_counts'}, axis=1)
    sgrna_assignment['guide_type'] = np.where(sgrna_assignment['guide_id'].str.startswith("NTC-"), 'non-targeting', 'targeting')
    sgrna_assignment = sgrna_assignment.set_index('cell')
    sgrna_assignment = sgrna_assignment[sgrna_assignment.index.isin(adata.obs_names)].copy()

    # Merge with scRNA-seq anndata
    existing_cols = [col for col in sgrna_assignment.columns if col in adata.obs.columns]
    adata.obs.drop(columns=existing_cols, inplace=True)
    adata.obs = pd.concat([adata.obs, sgrna_assignment], axis=1)

    return sgrna_assignment

def get_long_sgrna_umi_counts(grna_ad):
    """
    Convert sgRNA AnnData object to long format DataFrame with UMI counts.
    
    Parameters
    ----------
    grna_ad : anndata.AnnData
        AnnData object containing sgRNA UMI counts
        
    Returns
    -------
    pandas.DataFrame
        Long format DataFrame with columns:
        - cell: cell barcode
        - guide: guide ID
        - UMI_counts: number of UMIs
    """
    # Get the UMI count matrix
    X = grna_ad.X

    # Convert sparse matrix to coordinate format if needed
    if hasattr(X, 'tocoo'):
        X_coo = X.tocoo()
    else:
        # If X is dense, convert to sparse first
        X_coo = csr_matrix(X).tocoo()

    long_df = pd.DataFrame({
        'cell': grna_ad.obs_names[X_coo.row],
        'gRNA': grna_ad.var_names[X_coo.col], 
        'UMI_counts': X_coo.data
    })
    
    # Filter out zero counts and reset index
    long_df = long_df[long_df['UMI_counts'] > 0].reset_index(drop=True)
    
    return long_df

def get_background_vs_signal_guide_counts(datadir, return_summary=False):
    '''
    Analyze sgRNA counts to distinguish between background and signal for all samples in a directory.
    
    Args:
        datadir: Directory containing sgRNA assignment files
        return_summary: If True, also return a summary dataframe with mean UMI counts
        
    Returns:
        If return_summary is False: DataFrame with cell, gRNA, UMI counts and signal/background classification
        If return_summary is True: Tuple of (all_guides_df, summary_df)
    '''
    # Find all assignment files
    all_assignment_files = [f for f in os.listdir(datadir) if f.endswith('.sgrna_assignment_all.csv')]
    samples = set([f.split('.sgrna_assignment_all.csv')[0] for f in all_assignment_files])
    
    all_guides_combined = pd.DataFrame()
    all_multi_sgrna_cells = []
    for sample in samples:
        # Get all UMI counts (not just passing assignment thresholds)
        grna_ad = sc.read_h5ad(f'{datadir}/{sample}.sgRNA.h5ad')
        all_umi_counts = get_long_sgrna_umi_counts(grna_ad)
        
        clean_assignment = f'{datadir}/{sample}.sgrna_assignment_all.csv'
        assigned_guides = pd.read_csv(clean_assignment)

        all_guides = pd.merge(all_umi_counts, assigned_guides, on=['cell','gRNA'], how='left', suffixes=['_all', '_assigned'])
        all_guides['signal_vs_bg'] = np.where(all_guides['UMI_counts_assigned'].isna(), 'background', 'signal')
        all_guides['sample'] = sample
        
        all_guides_combined = pd.concat([all_guides_combined, all_guides], ignore_index=True)
    
    if return_summary:
        # Exclude multi guide, because we can't distinguish background here
        all_guides_combined_sum = all_guides_combined[~all_guides_combined['cell'].isin(all_multi_sgrna_cells)]
        all_guides_combined_sum = all_guides_combined_sum[['cell', 'gRNA', 'UMI_counts_all', 'signal_vs_bg']].copy()
        summary_df = all_guides_combined_sum.groupby(['signal_vs_bg', 'gRNA'])['UMI_counts_all'].median().reset_index()
        summary_df = summary_df.pivot_table(columns='signal_vs_bg', index=['gRNA'], values='UMI_counts_all')
        summary_df = summary_df.fillna(0)
        return(all_guides_combined, summary_df)
    else:
        return(all_guides_combined)


if __name__ == "__main__":
    import argparse
    import yaml
    import glob
    
    parser = argparse.ArgumentParser(description='Assign sgRNAs for GWT experiment')
    parser.add_argument('experiment', type=str,
                       help='Experiment ID to process. If not specified, processes all experiments')
    parser.add_argument('--config', type=str,
                       default='../../metadata/experiments_config.yaml',
                       help='Path to experiment config YAML file')
    parser.add_argument('--plot_dir', type=str, default='/scratch/groups/pritch/emma/sgrna_assignment_outs/',
                       help='Path to store plots')
    parser.add_argument('--n_cores', type=int, default=5,
                       help='Number of cores')
    parser.add_argument('--n_guides_parallel', type=int, default=10,
                       help='Number of gRNAs to process in parallel')
    parser.add_argument('--crispr_h5ad', type=str, default=None,
                       help='Process only a specific CRISPR h5ad file')
                       
    # arguments for chunking
    parser.add_argument('--start_idx', type=int, default=0,
                       help='Start index for variable chunk processing')
    parser.add_argument('--end_idx', type=int, default=None,
                       help='End index for variable chunk processing')
    parser.add_argument('--chunk_id', type=str, default=None,
                       help='Identifier for this chunk (used in output filenames)')
    parser.add_argument('--merge', action='store_true',
                       help='Merge chunk results into final assignment file.')
    parser.add_argument('--expected_n_chunks',type=int, default=18,
                       help='Expected chunk number for merge.')
    args = parser.parse_args()

    # Load config file
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    experiment = args.experiment
    config = config[experiment]    

    datadir = _convert_oak_path(config['datadir'])
 
    # Find all CRISPR h5ad files or use the specified one
    if args.crispr_h5ad:
        crispr_files = [args.crispr_h5ad]
    else:
        crispr_files = glob.glob(f"{datadir}/**/*.sgRNA.h5ad", recursive=True)
        if not crispr_files:
            raise ValueError(f"No .sgRNA.h5ad files found in {datadir}")
        print(f"Found {len(crispr_files)} CRISPR h5ad files")


    for sgrna_h5ad in crispr_files:
        # Extract sample identifier from file path
        sample_lane_id = os.path.basename(sgrna_h5ad).replace('.sgRNA.h5ad', '')
        print(f"\nProcessing {sample_lane_id} from {sgrna_h5ad}")
        
        if args.merge:
            # Merge all chunk results
            print(f"Merging chunk results for sample {sample_lane_id}...")
            chunk_files = glob.glob(f'{datadir}/tmp/{sample_lane_id}.sgrna_assignment_chunk_*.csv')
            if not chunk_files:
                print(f"No chunk files found for sample {sample_lane_id}")
                continue
            if len(chunk_files) < args.expected_n_chunks:
                print(f"Missing chunk files found for sample {sample_lane_id}")
                continue
                
            try:
                # Check chunk numbers and expected chunk size from filenames
                chunk_info = [re.search(r'chunk_(\d+)_(\d+)\.csv$', f) for f in chunk_files]
                chunk_numbers = [int(m.group(1)) for m in chunk_info if m]
                chunk_numbers.sort()
                # Compute cumulative difference between chunk numbers to check for gaps
                chunk_diffs = np.diff(chunk_numbers)
                if len(np.unique(chunk_diffs)) > 1:
                    raise FileNotFoundError(f'missing chunk files for sample {sample_lane_id}')
                
                all_perturbations = []
                for chunk_file in sorted(chunk_files):
                    chunk_data = pd.read_csv(chunk_file)
                    if len(chunk_data) > 0:
                        all_perturbations.append(chunk_data)
                    else:
                        raise ValueError(f"Warning: Empty chunk file found: {chunk_file}")
                
                if not all_perturbations:
                    print(f"No valid perturbation data found for sample {sample_lane_id}")
                    continue
                
                perturbations = pd.concat(all_perturbations, ignore_index=True)
                
                # Make unique cell assignment
                assignment_size = perturbations.groupby('cell').size()
                assignment_crispat = perturbations.groupby('cell').apply(
                    lambda x: x.loc[x['UMI_counts'].idxmax()]
                ).reset_index(drop=True)
                assignment_crispat['cell'] = perturbations.groupby('cell').apply(lambda x: x.name).values
                assignment_crispat['guide_id'] = np.where(
                    assignment_size[assignment_crispat['cell']].values > 1,
                    'multi_sgRNA',
                    assignment_crispat['gRNA']
                )
                
                # Save final merged results
                output_dir = os.path.dirname(sgrna_h5ad)
                perturbations.to_csv(f'{output_dir}/{sample_lane_id}.sgrna_assignment_all.csv', index=False)
                assignment_crispat.to_csv(f'{output_dir}/{sample_lane_id}.sgrna_assignment.csv', index=False)
                print(f"Successfully saved merged assignments to {output_dir}/{sample_lane_id}.sgrna_assignment.csv")
                print(f"Total cells assigned: {len(assignment_crispat)}")
                print(f"Unique guides assigned: {assignment_crispat['guide_id'].nunique()}")

                # # Remove chunk files after successful merge
                # for chunk_file in chunk_files:
                #     os.remove(chunk_file)
                    
            except Exception as e:
                print(f"Error merging results for sample {sample_lane_id}: {str(e)}")
                continue
        
        else:
            # Create sample-specific output directory
            results_dir = f'{args.plot_dir}/{experiment}/'
            os.makedirs(results_dir + 'sgrna_assignment_crispat/', exist_ok=True)
    
            sample_output_dir = results_dir + f'sgrna_assignment_crispat/{sample_lane_id}/'
            os.makedirs(sample_output_dir, exist_ok=True)
            os.makedirs(sample_output_dir + "loss_plots/", exist_ok=True)
            os.makedirs(sample_output_dir + "fitted_model_plots/", exist_ok=True)

            crispr_a = sc.read_h5ad(sgrna_h5ad)
            crispr_a = crispr_a[:, crispr_a.var['n_cells'] > 3].copy()
            n_guides_parallel=args.n_guides_parallel
            num_cores=args.n_cores

            assignment_s, perturbations, thresholds = assign_sgrna_crispat(
                crispr_a,  
                output_dir=sample_output_dir, 
                start_idx=args.start_idx,
                end_idx=args.end_idx,
                n_guides_parallel=n_guides_parallel, num_cores=num_cores)

            # Save chunk-specific results
            chunk_suffix = f"chunk_{args.start_idx}_{args.end_idx}" if args.chunk_id is not None else ""
            
            if len(perturbations) > 0:
                if args.end_idx is not None:
                    perturbations.to_csv(f'{datadir}/tmp/{sample_lane_id}.sgrna_assignment_{chunk_suffix}.csv', index=False)
                else:
                    assignment_s.to_csv(f'{datadir}/tmp/{sample_lane_id}.sgrna_assignment.csv', index=False)