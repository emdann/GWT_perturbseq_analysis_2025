import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc
import glob
from tqdm import tqdm
import scipy
from itertools import combinations

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from copy import deepcopy

# Add the parent directory to the path to import from sibling directory
sys.path.append(os.path.abspath('../'))
from MultiStatePerturbSeqDataset import *

def _run_DE_test(pbulk_adata, test_state = 'Rest'):
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
        test_state = [test_state],
        min_counts_per_gene = min_counts_per_gene,
        return_model = False,
        n_cpus=3
        )

    n_cells_target = ms_perturb_data.adata.obs.groupby('target')['n_cells'].sum().reset_index()
    results = pd.merge(results.rename({'contrast':'target'}, axis=1), n_cells_target)
    results['n_donors'] = n_donors
    results['donors'] = '_'.join(pbulk_adata.obs['donor_id'].unique().tolist())
    results['signif'] = results['adj_p_value'] < 0.1
    return(results)



import argparse

parser = argparse.ArgumentParser(description='Run donor robustness DE analysis.')
parser.add_argument('--condition', type=str, required=True, help='Culture condition to test')
parser.add_argument('--target', type=str, required=True, help='Perturbed gene(s), comma-separated if multiple')
parser.add_argument('--datadir', type=str, default='/mnt/oak/users/emma/data/GWT/CD4i_final/', help='Directory containing data')
parser.add_argument('--experiment_name', type=str, default='CD4i_final', help='Name of the experiment')
args = parser.parse_args()

cond = args.condition
selected_perturbed_genes = args.target.split(',')

datadir = args.datadir
experiment_name = args.experiment_name
results_dir = datadir + 'donor_robustness_analysis/'
if not os.path.exists(results_dir):
    os.mkdir(results_dir)


pbulk_adata = anndata.experimental.read_lazy(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad')
keep = (pbulk_adata.obs['perturbed_gene_name'].isin(selected_perturbed_genes + ['NTC'])) & (pbulk_adata.obs['culture_condition'] == cond)
pbulk_adata_test = pbulk_adata[keep].to_memory()
all_donors = pbulk_adata_test.obs['donor_id'].unique().tolist()

# Generate all possible combinations of 2 donors
donor_combinations = []
for train_donors in combinations(all_donors, 2):
    donor_combinations.append(
        list(train_donors)
    )

donor_robustness_results = pd.DataFrame()
all_results_df = pd.DataFrame()
for train_ds in donor_combinations:
    # Run DE tests
    results_df = _run_DE_test(pbulk_adata_test[pbulk_adata_test.obs['donor_id'].isin(train_ds)], test_state=cond)
    all_results_df = pd.concat([all_results_df, results_df])

for t in all_results_df.target.unique():
    all_results_df[all_results_df['target'] == t].to_parquet(f'{results_dir}/DE_donor_robustness.{t}_{cond}.parquet')