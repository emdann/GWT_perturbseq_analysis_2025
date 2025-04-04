'''
Preprocess and merge Perturb-seq data.
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

def _process_cellranger_h5(f):
    a = sc.read_10x_h5(f, gex_only=False)
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


def get_sgrna_qc_metrics(crispr_a, min_sgrna_counts=3, q=0.05):
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

def _convert_oak_path(path):
        """Helper function to convert oak paths between different mount points"""
        if not os.path.exists(path):
            return path.replace('/oak/stanford/groups/pritch/', '/mnt/oak/')
        return path

def process_experiment(exp_config):

    # Make directories
    datadir = _convert_oak_path(exp_config['datadir'])
    tmpdir = f'{datadir}/tmp/'
    os.makedirs(tmpdir, exist_ok=True)
    experiment_id = os.path.basename(os.path.normpath(datadir))
    
    # Read sample metadata
    sample_metadata_path = _convert_oak_path(exp_config['sample_metadata'])
    sample_metadata = pd.read_csv(sample_metadata_path)

    if exp_config['lane_ids'] is not None:
        all_lanes = exp_config['lane_ids']

    h5_files = []
    for lane in all_lanes:
        h5_files_lane = [f'{datadir}/cellranger_outs/{lane}/{f}' for f in os.listdir(f'{datadir}/cellranger_outs/{lane}/') 
                    if f.endswith('_sample_filtered_feature_bc_matrix.h5')]
        if not h5_files_lane:
            raise ValueError(f"No .h5 files found in {datadir}")
        
        # Check that h5_sample_names match values in sample_metadata['library_id']
        h5_sample_names = [f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0] for f in h5_files_lane]
        missing_samples = set(h5_sample_names) - set(sample_metadata['library_id'])
        if missing_samples:
            print(f"Warning: Found samples in data that are missing from metadata library_id: {missing_samples}")
            
            # Try to map sample names using sample_id_mapping if available
            if 'sample_id_mapping' in exp_config and exp_config['sample_id_mapping']:
                sample_id_mapping = exp_config['sample_id_mapping']
                
                # Check if all missing samples are in the mapping
                unmapped_samples = missing_samples - set(sample_id_mapping.keys())
                if unmapped_samples:
                    raise ValueError(f"Samples {unmapped_samples} not found in metadata or sample_id_mapping")
            else:
                raise ValueError(f"Samples {missing_samples} not found in metadata and no sample_id_mapping provided")
        
        h5_files.extend(h5_files_lane)

    print(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")
    try:
        adata = sc.read_h5ad(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")
    except FileNotFoundError:
        print(f"Processing {len(h5_files)} files...")
        adata = None
        
        for f in h5_files:
            print(f"Processing {f}")
            f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
            lane_id = f.split('/')[-2]
            f_sample_name = sample_id_mapping[f_sample_name]
            gex_a, crispr_a = _process_cellranger_h5(f)
            gex_a.obs['library_id'] = f_sample_name
            gex_a.obs['lane_id'] = lane_id
            crispr_a.obs['library_id'] = f_sample_name
            crispr_a.obs['lane_id'] = lane_id
            gex_a.obs_names = gex_a.obs_names + "_" + gex_a.obs['lane_id'] + "_" + gex_a.obs['library_id'] 
            crispr_a.obs_names = crispr_a.obs_names + "_" + crispr_a.obs['lane_id'] + "_" + crispr_a.obs['library_id']
            gex_a = _basic_qc(gex_a)
            
            # Process sgRNA adata
            get_sgrna_qc_metrics(crispr_a, min_sgrna_counts=3, q=0.05)
            crispr_a.write_h5ad(f'{datadir}/{f_sample_name}.{lane_id}.sgRNA.h5ad')
            
            if adata is None:
                adata = gex_a
            else:
                adata = adata.concatenate(gex_a, index_unique=None)
        
        # Add sample-level metadata
        missing_samples = set(adata.obs['library_id']) - set(sample_metadata['library_id'])
        if missing_samples:
            raise ValueError(f"Found samples in data that are missing from metadata: {missing_samples}")
        
        obs_names = adata.obs_names.copy()
        adata.obs = adata.obs.merge(sample_metadata, on='library_id', how='left')
        adata.obs.index = obs_names

        # Save merged objects
        print("Saving merged objects...")
        adata.var = adata.var[['gene_ids', 'gene_name', 'mt']].copy()
        adata.write(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")

    # Basic dim reduction analysis
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
        exp_config = config[exp_id]
        adata = process_experiment(exp_config)
        