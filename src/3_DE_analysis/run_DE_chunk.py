import os,sys
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import shutil

sys.path.append(os.path.abspath('../'))
from utils import _convert_oak_path

from copy import deepcopy
import argparse

from MultiStatePerturbSeqDataset import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('--config', type=str, required=True, help='Path to config YAML file')
    parser.add_argument('--test_chunk', type=str, required=True, help='Chunk of targets to test (integer or comma-separated list of integers)')
    parser.add_argument('--culture_condition', type=str, required=True, help='Which condition to test in')
    parser.add_argument('--n_cpus', type=int, default=10, help='N cpus')
    args = parser.parse_args()

    n_cpus = args.n_cpus

    # Load configuration from YAML file
    with open(args.config, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    # Extract parameters from config
    datadir = _convert_oak_path(config['datadir'])
    experiment_name = config['experiment_name']
    run_name = config.get('run_name', 'default')
    
    datadir = f'{datadir}/{experiment_name}'
    chunk_ix = args.test_chunk
    cond = args.culture_condition

    # Create directories for results
    de_results_dir = f'{datadir}/DE_results_{run_name}/'
    de_results_tmp_dir = f'{de_results_dir}/tmp/'
    os.makedirs(de_results_dir, exist_ok=True)
    os.makedirs(de_results_tmp_dir, exist_ok=True)

    print('Reading pseudobulk data...')
    pbulk_adata = anndata.read_h5ad(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad')
    pbulk_adata = pbulk_adata[(pbulk_adata.obs['keep_for_DE']) & (pbulk_adata.obs['culture_condition'] == cond)].copy() 

    # Read list of genes for DE testing
    try:
        with open(f'{datadir}/DE_test_genes.{cond}.txt', 'r') as f:
            de_test_genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(de_test_genes)} genes for DE testing")
    except FileNotFoundError:
        raise(FileNotFoundError, f"Warning: DE test genes file not found at {datadir}/DE_test_genes.txt - run feature selection first.")

    pbulk_adata = pbulk_adata[:, de_test_genes].copy()

    # Get targets to test
    target_chunk_matrix = pd.read_csv(f'{datadir}/DE_target2chunk.{cond}.csv.gz', compression='gzip', index_col=0)
    # Allow chunk_ix to be a comma-separated list of integers
    if isinstance(chunk_ix, str):
        chunk_ixs = [int(x) for x in chunk_ix.split(",")]
    elif isinstance(chunk_ix, int):
        chunk_ixs = [chunk_ix]
    elif isinstance(chunk_ix, list):
        chunk_ixs = chunk_ix
    else:
        raise ValueError("chunk_ix must be int, str, or list of ints")

    # Get DE parameters from config
    run_de_params = config.get('run_DE_params', {})
    design_formula = run_de_params.get('design_formula', '~ log10_n_cells + target')
    min_counts_per_gene = run_de_params.get('min_counts_per_gene', 10)

    # Copy the config file to the results directory for reproducibility (do once)
    print('Copying config file to results directory...')
    shutil.copy(args.config, os.path.join(de_results_dir, 'DE_config.yaml'))

    for chunk_ix in chunk_ixs:
        print(f'Processing chunk {chunk_ix}')
        test_targets = target_chunk_matrix.index[target_chunk_matrix[f'chunk_{chunk_ix}'] == 1].tolist()
        
        # Run DE analysis
        ms_perturb_data = MultistatePerturbSeqDataset(
            pbulk_adata,
            sample_cols = ['cell_sample_id'],
            perturbation_type = 'CRISPRi',
            target_col = 'perturbed_gene_id',
            sgrna_col = 'guide_id',
            state_col = 'culture_condition',
            control_level = 'NTC'
            )
        
        print(f"Running DE analysis on pseudobulk data with shape {pbulk_adata.shape}")
        print(f"Testing condition: {cond}")
        print(f"Number of targets to test: {len(test_targets)} (chunk {chunk_ix})")

        model, results = ms_perturb_data.run_target_DE(
            design_formula = design_formula,
            test_state = [cond], test_targets=test_targets,
            min_counts_per_gene = min_counts_per_gene,
            return_model = True,
            n_cpus=n_cpus
            )
        
        print('Saving results...')
        results.to_csv(de_results_tmp_dir + f"DE_results.{cond}.chunk_{chunk_ix}.csv.gz", compression='gzip')

        # Test effect of confounder covariates
        covs = design_formula.replace('~', '').strip().split("+")
        covs = [x.strip() for x in covs]
        covs = [x for x in covs if x != 'target']

        confounder_results_df = pd.DataFrame()
        if 'donor_id' in covs:
            all_donors = model.dds.obs['donor_id'].unique().tolist()[1:]
            for d in all_donors:
                conf_contrast = model.cond(donor_id=d)
                conf_res_df = model.test_contrasts(conf_contrast)        
                conf_res_df['contrast'] = f'donor_id_{d}'
                confounder_results_df = pd.concat([confounder_results_df, conf_res_df])
                
        if 'log10_n_cells' in covs:
            conf_contrast = model.cond(log10_n_cells=1)
            conf_res_df = model.test_contrasts(conf_contrast)        
            conf_res_df['contrast'] = f'log10_n_cells'
            confounder_results_df = pd.concat([confounder_results_df, conf_res_df])

        confounder_results_df.to_csv(de_results_tmp_dir + f"DE_results_confounders.{cond}.chunk_{chunk_ix}.csv.gz", compression='gzip')    
