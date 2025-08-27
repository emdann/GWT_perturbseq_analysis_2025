import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm.notebook import tqdm
import warnings
import argparse

from copy import deepcopy

# Add the parent directory to the path to import from sibling directory
sys.path.append(os.path.abspath('../../1_preprocess/'))
sys.path.append(os.path.abspath('../../'))
sys.path.append(os.path.abspath('../'))
from MultiStatePerturbSeqDataset import *

def _run_DE_test(pbulk_adata):
    pbulk_adata.obs['log10_n_cells'] = np.log10(pbulk_adata.obs['n_cells'])

    n_donors = pbulk_adata.obs['donor_id'].nunique()
    if n_donors > 1:
        design_formula = '~ log10_n_cells + donor_id + target'
    else:
        design_formula = '~ log10_n_cells + target'

    # pbulk_adata = pbulk_adata[:, de_test_genes].copy()
    min_counts_per_gene = 10
    
    ms_perturb_data = MultistatePerturbSeqDataset(
        pbulk_adata,
        sample_cols = ['cell_sample_id'],
        perturbation_type = 'CRISPRi',
        target_col = 'perturbed_gene_id',
        sgrna_col = 'guide_id',
        state_col = 'culture_condition',
        control_level = 'NTC'
        )

    results = ms_perturb_data.run_target_DE(
        design_formula = design_formula,
        test_state = ['Rest'],
        min_counts_per_gene = min_counts_per_gene,
        return_model = False,
        n_cpus=12
        )

    n_cells_target = ms_perturb_data.adata.obs.groupby('target')['n_cells'].sum().reset_index()
    results = pd.merge(results.rename({'contrast':'target'}, axis=1), n_cells_target)
    results['n_donors'] = n_donors
    results['donors'] = '_'.join(pbulk_adata.obs['donor_id'].unique().tolist())
    results['signif'] = results['adj_p_value'] < 0.1
    return(results)


def main(datadir = '.'):
    parser = argparse.ArgumentParser()
    parser.add_argument('k', help='power analysis key')
    parser.add_argument('--depth_perc', type=int, default = 100, help='testing params')
    parser.add_argument('--test', action='store_true', help='testing params')
    args = parser.parse_args()
    k = args.k

    k_adata = sc.read_h5ad(f'{datadir}/power_analysis.{k}.h5ad')
    if args.depth_perc == 100:
        if 'sum' in k_adata.layers:
            k_adata.X = k_adata.layers['sum'].copy()
        k_adata.layers.clear()
    else:
        k_adata.X = k_adata.layers[f'sum_downsampled_{args.depth_perc}perc'].copy()
        k_adata.layers.clear()
    donor_comb_cols = [x for x in k_adata.obs.columns if x.startswith('donor_comb')]
    if args.test:
        donor_comb_cols = donor_comb_cols[0:3]
    all_results = []
    for comb in donor_comb_cols: 
        try:
            # Convert the specific warning to an exception
            with warnings.catch_warnings():
                warnings.filterwarnings("error", 
                                        message="Every gene contains at least one zero, cannot compute log geometric means.*",
                                        category=UserWarning,
                                        module="pydeseq2.dds")

                results = _run_DE_test(k_adata[k_adata.obs[comb]])
                results['donor_comb'] = comb
                results['power_analysis_key'] = k
                all_results.append(results)
        except UserWarning as e:
            print(f"Caught warning for iteration {comb}: {e}")
            continue
    pd.concat(all_results).to_csv(f'{datadir}/power_analysis.{k}_{args.depth_perc}perc.DE_results.csv')

if __name__ == '__main__':
    main()
