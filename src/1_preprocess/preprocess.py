'''
Preprocess and merge CITE-seq data.

conda activate gwc-env
cd /oak/stanford/groups/pritch/users/emma/bin/tcell_perturbseq_analysis/_TcellsGW_pilot_analysis/
EXP_DIR=/oak/stanford/groups/pritch/users/emma/data/TcellsGW_PilotD2Redo_newRHSprobe_sample1/
sbatch --partition=pritch \
       --job-name=preprocess-TcellsGW \
       --cpus-per-task=3 \
       --mem=75000 \
       --time=2:00:00 \
       --output=$GROUP_SCRATCH/emma/slurm-pp-%j.out \
       --error=$GROUP_SCRATCH/emma/slurm-pp-%j.err \
       --wrap="python preprocess.py --datadir $EXP_DIR"
'''
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Define valid values
TIMEPOINTS = ["Day7", "Rest", "Restim_6hr", "Restim_24hr"]
HYBRID_PARAMS = ["half_million", "one_million", "one_and_half_million", "three_quarter_million"]

def ga_umi(adata_crispr, threshold = 3):   
    '''
    Guide assignment with fixed UMI thresholds 
    
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        thresholds (list): list of integers to use as thresholds (create assignment output file for each t in the list)
        output_dir (str): directory in which to store the resulting assignment
        
    Returns:
        None
    '''
    gRNA_list = adata_crispr.var_names.tolist()

    # Get perturbed cells for each gRNA based on a fixed UMI threshold
    perturbations = pd.DataFrame({'cell': [], 'gRNA': []})
    print('Get perturbed cells for each gRNA with UMI threshold = ' + str(threshold))
    for gRNA in gRNA_list:
        # Get cells with UMI counts higher than the threshold for specified gRNA
        selected_guide = adata_crispr[:,[gRNA]].X
        perturbed_cells = adata_crispr.obs_names[selected_guide.toarray().reshape(-1) >= threshold].tolist()
        UMI_counts = adata_crispr[selected_guide.toarray().reshape(-1) >= threshold, [gRNA]].X.toarray().reshape(-1)
        
        if len(perturbed_cells) != 0:
            df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA, 'UMI_counts': UMI_counts})
            perturbations = pd.concat([perturbations, df], ignore_index = True)
    return(perturbations)

def extract_metadata(sample_name):
    """
    Extract timepoint and hybridization ID from sample names using predefined valid values.
    
    Args:
        sample_name (str): Full sample name (e.g., "Day7_half_million")
        
    Returns:
        tuple: (timepoint, hybridization_id)
        
    Raises:
        ValueError: If timepoint or hybridization parameter is not recognized
    """
    # First try to find which timepoint is in the sample name
    found_timepoint = None
    for timepoint in TIMEPOINTS:
        if sample_name.startswith(timepoint+"_"):
            found_timepoint = timepoint
            remaining = sample_name[len(timepoint)+1:]  # +1 for the underscore
            break
    
    if not found_timepoint:
        raise ValueError(f"Unknown timepoint in sample name: {sample_name}\nValid timepoints are: {TIMEPOINTS}")
    
    # Find hybridization parameter in remaining string
    found_hybid = None
    for hybid in HYBRID_PARAMS:
        if remaining == hybid:
            found_hybid = hybid
            break
    
    if not found_hybid:
        raise ValueError(f"Unknown hybridization parameter in sample name: {remaining}\nValid parameters are: {HYBRID_PARAMS}")
    
    return found_timepoint, found_hybid

def _process_cellranger(f):
    try:
        f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
    except IndexError:
        raise ValueError(f"Filename {f} doesn't match expected format")
    
    a = sc.read_10x_h5(f, gex_only=False)
    
    # f_timepoint,f_hybrid = extract_metadata(f_sample_name)
    # a.obs['hybridization_condition'] = f_hybrid
    # a.obs['timepoint'] = f_timepoint
    # a.obs['sample_id'] = f'{f_timepoint}_{f_hybrid}'
    a.obs['sample_id'] = f'{f_sample_name}'
    a.obs_names = a.obs_names + "_" + a.obs['sample_id']

    # split by modality (for Perturb-seq)
    if not all(a.var['feature_types'] == 'Gene Expression'):
        gex_a = a[:, a.var['feature_types'] == 'Gene Expression'].copy()
        gex_a.var.drop(['pattern', 'read', 'sequence'], axis=1, inplace=True, errors='ignore')
        gex_a.var['gene_name'] = gex_a.var_names.values
        gex_a.var_names = gex_a.var['gene_ids'].values
        crispr_a = a[:, a.var['feature_types'] != 'Gene Expression'].copy()
        return(gex_a, crispr_a)
    else:
        a.var.drop(['pattern', 'read', 'sequence'], axis=1, inplace=True, errors='ignore')
        a.var['gene_name'] = a.var_names.values
        a.var_names = a.var['gene_ids'].values
        return(a)

def _compute_nonzero_means_v1(X_mat):
    if X_mat.format == 'csc':
        nnz_per_col = np.diff(X_mat.indptr)
    else:
        # Convert to CSC temporarily for column operations
        X_csc = X_mat.tocsc()
        nnz_per_col = np.diff(X_csc.indptr)
    
    col_sums = np.asarray(X_mat.sum(axis=0)).ravel()
    
    return np.divide(col_sums, nnz_per_col, 
                    out=np.zeros_like(col_sums), 
                    where=nnz_per_col != 0)


def get_sgrna_qc_metrics(crispr_a):
    var_cols = ['sgrna_id','perturbed_gene_name', 'perturbation_type','sgrna_type', 'feature_types', 'genome', 'pattern', 'read', 'sequence',
       'n_cells', 'mean_counts', 'total_counts', 'nonz_means']
    # Sanitize excel problems
    crispr_a.var_names = np.where(crispr_a.var_names == '1-Jun', 'JUN-1', crispr_a.var_names)
    crispr_a.var_names = np.where(crispr_a.var_names == '2-Jun', 'JUN-2', crispr_a.var_names)
    # Compute mean of non-zero UMIs
    crispr_a.var['nonz_means'] = _compute_nonzero_means_v1(crispr_a.X)
    sc.pp.calculate_qc_metrics(crispr_a, inplace=True)
    
    perturb_metadata = crispr_a.var.copy()
    # Annotate perturbed gene
    perturb_metadata['perturbed_gene_name'] = perturb_metadata.index.str.split('-').str[0:-1].str.join('-')
    perturb_metadata['perturbation_type'] = 'CRISPRi'
    perturb_metadata['sgrna_type'] = 'targeting'
    perturb_metadata['sgrna_type'] = np.where(perturb_metadata['perturbed_gene_name'] == 'NTC', 'NTC', perturb_metadata['sgrna_type'])
    perturb_metadata['sgrna_type'] = np.where(perturb_metadata['perturbed_gene_name'] == 'ProbeNTC', 'ProbeNTC', perturb_metadata['sgrna_type'])
    perturb_metadata['sgrna_id'] = perturb_metadata.index
    perturb_metadata = perturb_metadata.rename(
        {'n_cells_by_counts':'n_cells'}, 
        axis=1) 
    crispr_a.var = perturb_metadata[var_cols].copy()

def _basic_qc(adata, filter_cells=True):
    """
    Perform basic QC on gene expression data
    """
    ## Basic QC metrics
    adata.var["mt"] = adata.var['gene_name'].str.startswith("MT-")  # "MT-" for human, "Mt-" for mouse
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

    if filter_cells:
        print(f"Cells before filtering: {adata.n_obs}")
        adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
        sc.pp.filter_cells(adata, min_genes=200)
        print(f"Cells after filtering: {adata.n_obs}")
    
    return adata

def main(datadir):
    tmpdir = f'{datadir}/tmp/'
    os.makedirs(tmpdir, exist_ok=True)
    experiment_id = os.path.basename(os.path.normpath(datadir))

    # Chekc if f"{tmpdir}TcellsGW_pilot_merged.gex.h5ad" exists
    # if os.path.exists(f"{tmpdir}/{experiment_id}_merged.gex.h5ad"):
    #     print(f"{tmpdir}/{experiment_id}_merged.gex.h5ad already exists. Skipping...")
    #     adata = sc.read_h5ad(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")
    # else:        
    # Process files
    h5_files = [f'{datadir}/cellranger_outs/{f}' for f in os.listdir(f'{datadir}/cellranger_outs/') 
                if f.endswith('_sample_filtered_feature_bc_matrix.h5')]
    
    if not h5_files:
        raise ValueError(f"No .h5 files found in {datadir}")
    
    print(f"Processing {len(h5_files)} files...")
    adata = None
    
    for f in h5_files:
        print(f"Processing {f}")
        f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
        gex_a, crispr_a = _process_cellranger(f)
        gex_a = _basic_qc(gex_a)
        get_sgrna_qc_metrics(crispr_a)
        crispr_a.write_h5ad(f'{datadir}/{f_sample_name}.sgRNA.h5ad')
        
        if adata is None:
            adata = gex_a
        else:
            adata = adata.concatenate(gex_a, index_unique=None)
    
    # Save merged objects
    print("Saving merged objects...")
    adata.write(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")

    # Basic dim reduction analysis
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    sc.pp.pca(adata, n_comps=50)
    adata.write(f"{datadir}/{experiment_id}_merged.gex.lognorm.h5ad")
    
    return adata
if __name__ == "__main__":
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser(description='Process GWT experiment')
    parser.add_argument('--config', type=str,
                       default='/oak/stanford/groups/pritch/users/emma/bin/GWT_perturbseq_analysis/metadata/experiment_config.yaml',
                       help='Path to experiment config YAML file')
    parser.add_argument('--experiment', type=str,
                       help='Experiment ID to process. If not specified, processes all experiments')
    args = parser.parse_args()

    # Load config file
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    if args.experiment:
        if args.experiment not in config:
            raise ValueError(f"Experiment {args.experiment} not found in config file")
        experiments = [args.experiment]
    else:
        experiments = list(config.keys())
    
    for exp_id in experiments:
        print(f"\nProcessing experiment: {exp_id}")
        datadir = config[exp_id]['datadir']
        if not os.path.exists(datadir):
            datadir = datadir.replace('/oak/stanford/groups/pritch/', '/mnt/oak/')
        adata = main(datadir=datadir)
        