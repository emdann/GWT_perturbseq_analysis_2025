import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import argparse


def make_pseudobulks(datadir, experiment_name, condition_col='culture_condition', donor_col='donor_id', perturbation_col= 'perturbed_gene_name'):
    """
    Create pseudobulk data from single-cell RNA-seq data.
    
    Parameters:
    -----------
    datadir : str
        Directory containing the data files
    experiment_name : str
        Name of the experiment
    condition_col : str, optional
        Column name for condition information, default is 'culture_condition'
    donor_col : str, optional
        Column name for donor information, default is 'donor_id'
    perturbation_col : str, optional
        Column name for perturbation information, default is 'perturbed_gene_name'
    """
    # Create a subdirectory for perturbation prediction results
    perturbation_prediction_dir = os.path.join(datadir, "perturbation_prediction")
    
    # Check if directory exists, if not create it
    if not os.path.exists(perturbation_prediction_dir):
        os.makedirs(perturbation_prediction_dir)
    # Validate input parameters
    if not os.path.exists(datadir):
        raise ValueError(f"Data directory '{datadir}' does not exist")
    
    # Check if required files exist
    if experiment_name != 'K562_gwps':
        h5ad_file = f'{datadir}/{experiment_name}_merged.gex.h5ad'
        obs_file = f'{datadir}/{experiment_name}_merged.gex.lognorm.postQC_obs.csv'
    
        if not os.path.exists(h5ad_file):
            raise FileNotFoundError(f"Required file not found: {h5ad_file}")
        
        if not os.path.exists(obs_file):
            raise FileNotFoundError(f"Required file not found: {obs_file}")
        
        # Validate column names will be present in the data
        try:
            obs_df = pd.read_csv(obs_file, compression='gzip', index_col=0, nrows=1)
            missing_cols = []
            for col in [condition_col, donor_col, perturbation_col, 'QC_mask', 'perturbed_gene_id']:
                if col not in obs_df.columns and col != 'perturbed_gene_id':  # perturbed_gene_id might be created later
                    missing_cols.append(col)
            
            if missing_cols:
                raise ValueError(f"The following required columns are missing in the observation data: {', '.join(missing_cols)}")
        except Exception as e:
            raise ValueError(f"Error validating observation data: {str(e)}")
    
        print('Loading...')
        adata = sc.read_h5ad(h5ad_file, backed=True)
        adata.obs = pd.read_csv(obs_file, compression='gzip', index_col=0)
        adata = adata[adata.obs['QC_mask']].to_memory()
        adata = adata[adata.obs['guide_id'] != 'multi_sgRNA'].to_memory()
         # Create a new column joining donor and condition with "-" as separator
        adata.obs['pbulk_sample'] = adata.obs[donor_col].astype(str) + "-" + adata.obs[condition_col].astype(str)
        sample_col = 'pbulk_sample'
    else:
        h5ad_file = '{datadir}/K562_gwps_full.h5ad'
        adata = sc.read_h5ad(h5ad_file)
        sample_col = 'cell_line'

    # Normalize
    print('Normalizing...')
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    for s in adata.obs[sample_col].unique():
        adata_s = adata[adata.obs[sample_col] == s].copy()
        adata_bulk = sc.get.aggregate(adata_s, by=perturbation_col, func=['mean', 'count_nonzero'])
        adata_bulk.X = adata_bulk.layers['mean']
        del adata_bulk.layers['mean']
        adata_bulk.write_h5ad(f"{perturbation_prediction_dir}/{experiment_name}.prediction_pseudobulk.{s}.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('--datadir', type=str, required=True, help='Directory containing the data files')
    parser.add_argument('--experiment_name', type=str, required=True, help='Name of the experiment')
    parser.add_argument('--condition_col', type=str, default='culture_condition', help='Column name for condition information')
    parser.add_argument('--donor_col', type=str, default='donor_id', help='Column name for donor information')
    parser.add_argument('--perturb_col', type=str, default='perturbed_gene_name', help='Column name for perturbation information')
    
    args = parser.parse_args()
    
    make_pseudobulks(
        args.datadir,
        args.experiment_name,
        args.condition_col,
        args.donor_col,
        args.perturb_col,
    )