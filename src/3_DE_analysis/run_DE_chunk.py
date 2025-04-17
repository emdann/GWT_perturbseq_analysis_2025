import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

from copy import deepcopy
import argparse

from MultiStatePerturbSeqDataset import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('--datadir', type=str, required=True, help='Directory containing the data files')
    parser.add_argument('--experiment_name', type=str, required=True, help='Name of the experiment')
    parser.add_argument('--test_chunk', type=int, required=True, help='Chunk of targets to test')
    parser.add_argument('--culture_condition', type=str, required=True, help='Which condition to test in')
    args = parser.parse_args()

    datadir = f'{args.datadir}/{args.experiment_name}'
    experiment_name = args.experiment_name
    chunk_ix = args.test_chunk
    cond = args.culture_condition

    # Create directories for results
    de_results_dir = f'{datadir}/DE_results/'
    de_results_tmp_dir = f'{de_results_dir}/tmp/'
    os.makedirs(de_results_dir, exist_ok=True)
    os.makedirs(de_results_tmp_dir, exist_ok=True)

    pbulk_adata = sc.read_h5ad(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad')
    pbulk_adata = pbulk_adata[pbulk_adata.obs['keep_for_DE']].copy()

    # Read list of genes for DE testing
    try:
        with open(f'{datadir}/DE_test_genes.txt', 'r') as f:
            de_test_genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(de_test_genes)} genes for DE testing")
    except FileNotFoundError:
        raise(FileNotFoundError, f"Warning: DE test genes file not found at {datadir}/DE_test_genes.txt - run feature selection first.")
    
    pbulk_adata = pbulk_adata[:, de_test_genes].copy()

    # Get targets to test
    target_chunk_matrix = pd.read_csv(f'{datadir}/DE_target2chunk.csv.gz', compression='gzip', index_col=0)
    test_targets = target_chunk_matrix.index[target_chunk_matrix[f'chunk_{chunk_ix}'] == 1].tolist()
    
    # Run DE analysis
    ms_perturb_data = MultistatePerturbSeqDataset(pbulk_adata)
    
    results = ms_perturb_data.run_target_DE(
        design_formula = '~ donor_id + target',
        test_state = [cond], test_targets=test_targets)
    
    print('Saving results...')
    results.to_csv(de_results_tmp_dir + f"DE_results.{cond}.chunk_{chunk_ix}.csv.gz", compression='gzip')
    
    