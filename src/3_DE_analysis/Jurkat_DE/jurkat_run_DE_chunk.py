'''
To run on comino:
for chunk in {0..67}; do sbatch --mem=16G --wrap="PYTHONNOUSERSITE=1 /home/emmadann/miniforge3/envs/gwt_de/bin/python Jurkat_DE/jurkat_run_DE_chunk.py --chunk_ix ${chunk}"; done
'''
import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scipy

import argparse
sys.path.append('../')
from MultiStatePerturbSeqDataset import MultistatePerturbSeqDataset

def parse_args():
    parser = argparse.ArgumentParser(description='Run DE analysis on Jurkat pseudobulk data chunk')
    parser.add_argument('--chunk_ix', type=int, required=True, help='Chunk index to process')
    return parser.parse_args()

def main():
    args = parse_args()
    chunk_ix = args.chunk_ix

    pbulk_adata = sc.read_h5ad('/mnt/oak/users/emma/data/GWT/Jurkat_DE_analysis/Jurkat_essential_filt_pseudobulk.for_DE.h5ad')
    target_chunk_matrix = pd.read_csv('/mnt/oak/users/emma/data/GWT/Jurkat_DE_analysis/Jurkat_essential_filt_pseudobulk.target2chunk.csv', index_col=0)
    outdir = '/mnt/oak/users/emma/data/GWT/Jurkat_DE_analysis/'
    de_results_tmp_dir = f'{outdir}/tmp/'

    pbulk_adata.X = pbulk_adata.layers['sum'].copy()
    pbulk_adata.X.data = np.round(pbulk_adata.X.data).astype(np.int32)
    del pbulk_adata.layers['sum']
    pbulk_adata.obs['culture_condition'] = 'jurkat'
    pbulk_adata.obs['guide_id'] = pbulk_adata.obs['gene'].copy()
    pbulk_adata.obs['gem_group'] = pbulk_adata.obs['gem_group'].astype(str)

    # Run DE analysis
    ms_perturb_data = MultistatePerturbSeqDataset(
        pbulk_adata,
        sample_cols = ['gem_group'],
        perturbation_type = 'CRISPRi',
        target_col = 'gene',
        sgrna_col = 'guide_id',
        state_col = 'culture_condition',
        control_level = 'non-targeting'
        )

    test_targets = target_chunk_matrix.index[target_chunk_matrix[f'chunk_{chunk_ix}'] == 1].tolist()

    design_formula = '~ target'
    model, results = ms_perturb_data.run_target_DE(
            design_formula = design_formula,
            test_targets=test_targets,
            min_counts_per_gene = 10,
            return_model = True
            )
    print('Saving results...')
    results.to_csv(de_results_tmp_dir + f"DE_results.chunk_{chunk_ix}.csv.gz", compression='gzip')

if __name__ == '__main__':
    main()
