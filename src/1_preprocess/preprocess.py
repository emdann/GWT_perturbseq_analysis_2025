'''
Preprocess and merge Perturb-seq data.
'''
import os,sys
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
from pathlib import Path
import glob
import warnings
from tqdm import tqdm
import xarray
warnings.filterwarnings('ignore')

def _read_cellranger_h5(f):
    try:
        a = sc.read_10x_h5(f, gex_only=False)
    except ValueError:
        a = sc.read_h5ad(f)
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

def process_cellranger_h5(f, datadir, sample_id_mapping):
    '''
    Process single cellranger file
    '''
    datadir = _convert_oak_path(datadir)
    f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
    lane_id = f.split('/')[-2]
    f_sample_name = sample_id_mapping[f_sample_name]
    processed_file = f'{datadir}/tmp/{f_sample_name}.{lane_id}.scRNA.h5ad'

    print(f"Processing new file: {f}")
    gex_a, crispr_a = _read_cellranger_h5(f)
    gex_a.obs['library_id'] = f_sample_name
    gex_a.obs['lane_id'] = lane_id
    crispr_a.obs['library_id'] = f_sample_name
    crispr_a.obs['lane_id'] = lane_id
    gex_a.obs_names = gex_a.obs_names + "_" + gex_a.obs['lane_id'] + "_" + gex_a.obs['library_id'] 
    crispr_a.obs_names = crispr_a.obs_names + "_" + crispr_a.obs['lane_id'] + "_" + crispr_a.obs['library_id']
    
    # Process sgRNA adata
    get_sgrna_qc_metrics(crispr_a, min_sgrna_counts=3, q=0.05)
    crispr_a.write_h5ad(f'{datadir}/{f_sample_name}.{lane_id}.sgRNA.h5ad')

    # Process scRNA adata
    gex_a.layers['counts'] = gex_a.X.copy()
    # rsc.get.anndata_to_GPU(gex_a)
    gex_a = _basic_qc(gex_a)
    sc.pp.normalize_total(gex_a)
    sc.pp.log1p(gex_a)
    gex_a.write_h5ad(processed_file)
    
    return(gex_a)

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
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], log1p=True, inplace=True)

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

    # Get all h5 files
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

    print(f"Processing {len(h5_files)} files...")    
    for f in tqdm(h5_files, desc="Processing h5 files"):
        # Check if output already exists
        f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
        lane_id = f.split('/')[-2]
        f_sample_name = sample_id_mapping[f_sample_name]
        processed_file = f'{datadir}/tmp/{f_sample_name}.{lane_id}.scRNA.h5ad'
        
        if os.path.exists(processed_file):
            print(f"Loading existing processed file: {processed_file}")
            gex_a = sc.read_h5ad(processed_file)
        else:
            print(f"Processing {f}")
            gex_a = process_cellranger_h5(f)
        
def merge_experiment(exp_config, ondisk=True):
    datadir = _convert_oak_path(exp_config['datadir'])
    tmpdir = f'{datadir}/tmp/'
    os.makedirs(tmpdir, exist_ok=True)
    experiment_id = os.path.basename(os.path.normpath(datadir))

    # Read sample metadata
    sample_metadata_path = _convert_oak_path(exp_config['sample_metadata'])
    sample_metadata = pd.read_csv(sample_metadata_path)

    if exp_config['lane_ids'] is not None:
        all_lanes = exp_config['lane_ids']

    all_library_ids = sample_metadata['library_id'].tolist()

    # Create list of expected h5ad files by combining all library IDs with all lane IDs
    expected_h5ad_files = []
    for library_id in all_library_ids:
        for lane_id in all_lanes:
            expected_file = f"{tmpdir}/{library_id}.{lane_id}.scRNA.h5ad"
            expected_h5ad_files.append(expected_file)

    # Find all processed h5ad files
    h5ad_files = glob.glob(f"{tmpdir}/*.scRNA.postQC.h5ad")
    if not h5ad_files:
        raise ValueError(f"No processed .scRNA.h5ad files found in {tmpdir}")
    print(f"Found {len(h5ad_files)}/{len(expected_h5ad_files)} processed h5ad files")

    if len(h5ad_files) < len(expected_h5ad_files):
        missing_files = set(expected_h5ad_files) - set(h5ad_files)
        warnings.warn(f"Missing {len(missing_files)} expected h5ad files: {missing_files}")

    # merge
    if ondisk:
        # h5ad_files = h5ad_files[0:4]
        data_dict = {x.split(tmpdir)[-1].split('.h5ad')[0]:x for x in h5ad_files}
        ad.experimental.concat_on_disk(
            data_dict, f"{tmpdir}/{experiment_id}_merged.gex.h5ad"
        )
        adata = ad.experimental.read_lazy(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")
    else:
        raise NotImplementedError('merging in memory is not implemented')
    
    # Merge sample-level metadata
    missing_samples = set(adata.obs['library_id'].values) - set(sample_metadata['library_id'])
    if missing_samples:
        raise ValueError(f"Found samples in data that are missing from metadata: {missing_samples}")
    
    obs_df = adata.obs.copy()
    obs_df = obs_df.set_coords('library_id')
    merged = xarray.merge([
        obs_df, 
        sample_metadata.set_index('library_id').to_xarray()
    ], join='left')
    merged = merged.reset_index('library_id').reset_coords('library_id').drop_dims('library_id')
    adata.obs = merged.to_dataframe()

    # Save merged objects
    print("Saving merged objects...")
    adata.var = adata.var[['gene_ids', 'gene_name', 'mt']].copy()
    adata.write(f"{tmpdir}/{experiment_id}_merged.gex.h5ad")
    return adata 

    return adata

def run_pca_rapids(adata):
    try:
        import rapids_singlecell as rsc
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.highly_variable_genes(adata)
        rsc.pp.pca(adata, n_comps=50)
    except:
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata, n_comps=50)
    

if __name__ == "__main__":
    import argparse
    import yaml
    
    parser = argparse.ArgumentParser(description='Process GWT experiment')
    parser.add_argument('--config', type=str,
                       default='/oak/stanford/groups/pritch/users/emma/bin/GWT_perturbseq_analysis/metadata/experiments_config.yaml',
                       help='Path to experiment config YAML file')
    parser.add_argument('--experiment', type=str,
                       help='Experiment ID to process. If not specified, processes all experiments')
    parser.add_argument('--input_h5', type=str, default=None,
                       help='Process only a specific cellranger output h5 file')
    parser.add_argument('--merge', action='store_true',
                       help='Merge processed files into a single experiment file')
    parser.add_argument('--embedding', action='store_true',
                       help='Run PCA on merged object')
    parser.add_argument('--force', action='store_true',
                       help='Force overwrite of existing files')
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
        if args.input_h5:     
            f = args.input_h5
            # Check if output file exists
            datadir = _convert_oak_path(exp_config['datadir'])
            sample_id_mapping = exp_config['sample_id_mapping']
            f_sample_name = f.split('/')[-1].split('_sample_filtered_feature_bc_matrix')[0]
            lane_id = f.split('/')[-2]
            f_sample_name = sample_id_mapping[f_sample_name]
            processed_file = f'{datadir}/tmp/{f_sample_name}.{lane_id}.scRNA.h5ad'
            
            if os.path.exists(processed_file) and not args.force:
                print(f"Skipping {args.input_h5} - output file {processed_file} already exists")
                continue
            else:
                process_cellranger_h5(f, exp_config['datadir'], exp_config['sample_id_mapping'])
        else:
            if not args.merge:
                process_experiment(exp_config, force=args.force)
            else:
                adata = merge_experiment(exp_config, force=args.force)
                print('Merged object saved.')
            
            if args.embedding:
                run_pca_rapids(adata)
                output_file = f"{datadir}/{args.experiment}_merged.gex.lognorm.h5ad"
                if os.path.exists(output_file) and not args.force:
                    print(f"Skipping embedding - output file {output_file} already exists")
                else:
                    adata.write(output_file)