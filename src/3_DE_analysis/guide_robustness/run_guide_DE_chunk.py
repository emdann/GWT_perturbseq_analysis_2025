import os
# Set threading limits before importing any scientific libraries
n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
os.environ['OMP_NUM_THREADS'] = str(n_cpus)
os.environ['MKL_NUM_THREADS'] = str(n_cpus)
os.environ['NUMBA_NUM_THREADS'] = str(n_cpus)
os.environ['OPENBLAS_NUM_THREADS'] = str(n_cpus)

# CRITICAL: Loky-specific limits
os.environ['LOKY_MAX_CPU_COUNT'] = str(n_cpus)
print(f'Using {n_cpus} CPUs')

import os,sys
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import shutil

sys.path.append(os.path.abspath('../../'))
from utils import _convert_oak_path

from copy import deepcopy
import argparse

from MultiStatePerturbSeqDataset import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DE test by guide')
    parser.add_argument('--config', type=str, required=True, help='Path to config YAML file')
    parser.add_argument('--test_chunk', type=int, required=True, help='Chunk of targets to test')
    parser.add_argument('--culture_condition', type=str, required=True, help='Which condition to test in')
    args = parser.parse_args()

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

    # pbulk_adata = sc.read_h5ad(f'{datadir}/{experiment_name}_merged_by_guide.DE_pseudobulk.h5ad')
    # pbulk_adata = pbulk_adata[(pbulk_adata.obs['keep_for_DE']) & (pbulk_adata.obs.culture_condition == cond)].copy()
    pbulk_adata = anndata.experimental.read_lazy(f'{datadir}/{experiment_name}_merged_by_guide.DE_pseudobulk.h5ad')
    pbulk_adata = pbulk_adata[(pbulk_adata.obs['keep_for_DE']) & (pbulk_adata.obs['culture_condition'] == cond)].to_memory() 
    pbulk_adata.obs['guide_test'] = np.where(pbulk_adata.obs['guide_id'].str.startswith('NTC'), 'NTC', pbulk_adata.obs['guide_id'])

    # Read list of genes for DE testing
    try:
        with open(f'{datadir}/DE_single_guide_test_genes.{cond}.txt', 'r') as f:
            de_test_genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(de_test_genes)} genes for DE testing")
    except FileNotFoundError:
        raise(FileNotFoundError, f"Warning: DE test genes file not found at {datadir}/DE_test_genes.txt - run feature selection first.")

    pbulk_adata = pbulk_adata[:, de_test_genes].copy()

    # Get targets to test
    target_chunk_matrix = pd.read_csv(f'{datadir}/DE_single_guide_target2chunk.{cond}.csv.gz', compression='gzip', index_col=0)
    test_targets = target_chunk_matrix.index[target_chunk_matrix[f'chunk_{chunk_ix}'] == 1].tolist()
    
    # Run DE analysis
    ms_perturb_data = MultistatePerturbSeqDataset(
        pbulk_adata,
        sample_cols = ['cell_sample_id'],
        perturbation_type = 'CRISPRi',
        target_col = 'guide_test',
        sgrna_col = 'guide_id',
        state_col = 'culture_condition',
        control_level = 'NTC'
        )
    
    # Get DE parameters from config
    run_de_params = config.get('run_DE_params', {})
    design_formula = run_de_params.get('design_formula', '~ log10_n_cells + target')
    min_counts_per_gene = run_de_params.get('min_counts_per_gene', 10)
    
    print(f"Running DE analysis on pseudobulk data with shape {pbulk_adata.shape}")
    print(f"Testing condition: {cond}")
    print(f"Number of targets to test: {len(test_targets)}")

    model, results = ms_perturb_data.run_target_DE(
        design_formula = design_formula,
        test_state = [cond], test_targets=test_targets,
        min_counts_per_gene = min_counts_per_gene,
        return_model = True,
        n_cpus = n_cpus
        )
    
    print('Saving results...')
    results.to_csv(de_results_tmp_dir + f"DE_results_by_guide.{cond}.chunk_{chunk_ix}.csv.gz", compression='gzip')
    
    # Copy the config file to the results directory for reproducibility
    print('Copying config file to results directory...')
    shutil.copy(args.config, os.path.join(de_results_dir, 'DE_config.yaml'))
