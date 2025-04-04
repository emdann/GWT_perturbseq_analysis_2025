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
import multiprocessing as mp

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
    parser.add_argument('--plot_dir', type=str, default='../../results/',
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
    args = parser.parse_args()

    # Load config file
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    experiment = args.experiment
    config = config[experiment]    

    datadir = _convert_oak_path(config['datadir'])
    results_dir = f'{args.plot_dir}/{experiment}/'
    os.makedirs(results_dir + 'sgrna_assignment_crispat/', exist_ok=True)
    # os.makedirs(results_dir + 'sgrna_assignment_crispat/' + "loss_plots/", exist_ok=True)
    # os.makedirs(results_dir + 'sgrna_assignment_crispat/' + "fitted_model_plots/", exist_ok=True)

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
        
        # Create sample-specific output directory
        sample_output_dir = results_dir + f'sgrna_assignment_crispat/{sample_lane_id}/'
        os.makedirs(sample_output_dir, exist_ok=True)
        os.makedirs(sample_output_dir + "loss_plots/", exist_ok=True)
        os.makedirs(sample_output_dir + "fitted_model_plots/", exist_ok=True)
        
        if args.merge:
            # Merge all chunk results
            print(f"Merging chunk results for sample {sample_lane_id}...")
            chunk_files = glob.glob(f'{datadir}/tmp/{sample_lane_id}.sgrna_assignment_chunk_*.csv')
            if not chunk_files:
                print(f"No chunk files found for sample {sample_lane_id}")
                continue
                
            try:
                all_perturbations = []
                for chunk_file in sorted(chunk_files):
                    chunk_data = pd.read_csv(chunk_file)
                    if len(chunk_data) > 0:
                        all_perturbations.append(chunk_data)
                    else:
                        print(f"Warning: Empty chunk file found: {chunk_file}")
                
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
                
            except Exception as e:
                print(f"Error merging results for sample {sample_lane_id}: {str(e)}")
                continue
        
        else:
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