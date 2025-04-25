import os,sys
import numpy as np
import anndata
import pandas as pd
import mudata as md
import scanpy as sc
import yaml
import argparse

import matplotlib.pyplot as plt
import seaborn as sns

# Add the parent directory to the path to import from sibling directory
sys.path.append(os.path.abspath('../'))
from utils import feature_selection, _convert_oak_path

# Parse command line arguments
parser = argparse.ArgumentParser(description='Prepare data for differential expression analysis')
parser.add_argument('--config', type=str, required=True, help='Path to config YAML file')
args = parser.parse_args()

# Load configuration from YAML file
with open(args.config, 'r') as config_file:
    config = yaml.safe_load(config_file)

# Extract parameters from config
datadir = _convert_oak_path(config['datadir'])
experiment_name = config['experiment_name']
datadir = f'{datadir}/{experiment_name}'

# Get parameters with defaults
min_replicates = config.get('min_replicates', 3)
min_cells_per_guide = config.get('min_cells_per_guide', 5)
n_hvgs = config.get('n_hvgs', 10000)
chunk_size = config.get('chunk_size', 50)
chunk_split_seed = config.get('chunk_split_seed', 1423)

# Feature selection parameters
feature_selection_params = config.get('feature_selection', {})
highx_min_mean_counts = feature_selection_params.get('highx_min_mean_counts', 2000)
highx_min_pct_dropouts = feature_selection_params.get('highx_min_pct_dropouts', 0.5)
lowx_max_pct_dropouts = feature_selection_params.get('lowx_max_pct_dropouts', 99.9)

# File paths
file_paths = config.get('file_paths', {})
# Define these params directly instead of reading from config
obs_file = f"{experiment_name}_merged.gex.lognorm.postQC_obs.csv"
pseudobulk_file = f"{experiment_name}_merged.DE_pseudobulk.h5ad"
no_effect_guides_file = file_paths.get('no_effect_guides_file', 'no_effect_guides.txt')
de_test_genes_file = file_paths.get('de_test_genes_file', 'DE_test_genes.txt')
target2chunk_file = file_paths.get('target2chunk_file', 'DE_target2chunk.csv.gz')

# Load data using configured file paths
obs_df = pd.read_csv(f'{datadir}/{obs_file}', compression='gzip', index_col=0)
pbulk_adata = sc.read_h5ad(f'{datadir}/{pseudobulk_file}')

## -- Add info on number of cells aggregated -- ##
n_cells = obs_df.groupby(['donor_id', 'culture_condition', 'guide_id']).size().reset_index()
n_cells.columns = ['donor_id', 'culture_condition', 'sgrna', 'n_cells']
n_cells['log10_n_cells'] = np.log10(n_cells['n_cells'])
pbulk_adata.obs = pd.merge(pbulk_adata.obs, n_cells)
pbulk_adata.obs['total_counts'] = np.array(pbulk_adata.X.sum(1)).flatten()

# Exclude samples with very low counts compared to others (threshold by condition)
tot_counts_lower_q = pbulk_adata.obs.groupby('culture_condition')['total_counts'].quantile(config.get('min_counts_quantile', 0.005)).to_dict()

# Load list of ineffective guides
no_effect_guides_file = os.path.join(datadir, 'no_effect_guides.txt')
no_effect_guides = []
try:
    with open(no_effect_guides_file, 'r') as f:
        no_effect_guides = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(no_effect_guides)} ineffective guides from {no_effect_guides_file}")
except FileNotFoundError:
    raise FileNotFoundError(f"No effect guides file not found at {no_effect_guides_file} - run src/1_preprocess/estimate_guide_effect.ipynb")

# Load info on perturbed genes to test
# Create separate boolean columns for each filtering condition
pbulk_adata.obs['keep_min_cells'] = pbulk_adata.obs['n_cells'] >= config.get('min_cells_per_guide', 5)  # Mark guides with at least n cells per condition & sample
pbulk_adata.obs['keep_not_multi_sgrna'] = pbulk_adata.obs['sgrna'] != 'multi_sgRNA'  # Don't count multi-sgRNA samples
pbulk_adata.obs['keep_effective_guides'] = ~pbulk_adata.obs['sgrna'].isin(no_effect_guides)  # Don't count ineffective guides as useful replicates

# Filter counts outliers
pbulk_adata.obs['keep_total_counts'] = True  # Initialize as all True
for c, thresh in tot_counts_lower_q.items():
    mask = pbulk_adata.obs['culture_condition'] == c
    pbulk_adata.obs.loc[mask, 'keep_total_counts'] = pbulk_adata.obs.loc[mask, 'total_counts'] >= thresh

# Combine all conditions into a single keep_for_DE column
pbulk_adata.obs['keep_for_DE'] = (
    pbulk_adata.obs['keep_min_cells'] & 
    pbulk_adata.obs['keep_not_multi_sgrna'] & 
    pbulk_adata.obs['keep_effective_guides'] & 
    pbulk_adata.obs['keep_total_counts']
)

guide_cell_counts = pbulk_adata.obs[['sgrna', 'target', 'culture_condition', 'n_cells', 'keep_for_DE']]
guide_cell_counts = guide_cell_counts[~guide_cell_counts['sgrna'].str.startswith('NTC-')] 

# Count number of replicates per perturbed gene in each condition
genes_replicates = guide_cell_counts.groupby(['target', 'culture_condition'])['keep_for_DE'].sum()\
    .sort_values()\
    .reset_index()\
    .pivot(index='target', columns='culture_condition', values='keep_for_DE')\
    .fillna(0)

genes2test = genes_replicates >= min_replicates
genes2test_dict = {col: genes2test.index[genes2test[col]].tolist() for col in genes2test.columns}


## --- FILTER SAMPLES FOR DE ANALYSIS --- ##
# Add name of perturbed gene
guide_to_gene_id = dict(zip(obs_df['guide_id'], obs_df['perturbed_gene_name']))
pbulk_adata.obs['perturbed_gene_name'] = pbulk_adata.obs['sgrna'].map(guide_to_gene_id)

# Test perturbed genes passing filters in all conditions
test_genes = np.intersect1d(*[x for x in genes2test_dict.values()]).tolist()
test_genes.append('NTC')
# Create a specific column for the test genes condition
pbulk_adata.obs['keep_test_genes'] = pbulk_adata.obs['target'].isin(test_genes)

# Aggregate with the existing keep_for_DE condition
pbulk_adata.obs['keep_for_DE'] = pbulk_adata.obs['keep_for_DE'] & pbulk_adata.obs['keep_test_genes']

pbulk_adata.write_h5ad(f'{datadir}/{pseudobulk_file}')
pbulk_adata = pbulk_adata[pbulk_adata.obs['keep_for_DE']].copy()

## --- FEATURE SELECTION --- ##
feature_selection_var = feature_selection(
    pbulk_adata,
    n_hvgs = n_hvgs,
    subset_adata=False,
    highx_min_mean_counts = highx_min_mean_counts,
    highx_min_pct_dropouts_by_counts = highx_min_pct_dropouts,
    lowx_max_pct_dropouts_by_counts = lowx_max_pct_dropouts,
    return_all=True
    )
# Remove PuroR gene if it exists in the feature selection results
if 'CUSTOM001_PuroR' in feature_selection_var.index:
    feature_selection_var = feature_selection_var.drop('CUSTOM001_PuroR')

all_targets = pbulk_adata.obs['target'].unique()
feature_selection_var['is_target'] = feature_selection_var.index.isin(all_targets)

DE_test_genes = feature_selection_var[feature_selection_var['highly_variable'] | feature_selection_var['is_target']].index.tolist()
# Save DE_test_genes to a text file
with open(f'{datadir}/{de_test_genes_file}', 'w') as f:
    for gene in DE_test_genes:
        f.write(f"{gene}\n")

## --- SPLIT PERTURBATIONS INTO CHUNKS --- ##
all_targets = pbulk_adata.obs['target'].unique().tolist()
all_targets.remove('NTC')

# Randomize targets before splitting (without replacement)
np.random.seed(chunk_split_seed)
np.random.shuffle(all_targets)

# Split all_targets into groups based on chunk_size
target_chunks = [all_targets[i:i+chunk_size] for i in range(0, len(all_targets), chunk_size)]

# Initialize a binary matrix with zeros
target_chunk_matrix = pd.DataFrame(0, 
                                  index=all_targets, 
                                  columns=[f'chunk_{i}' for i in range(len(target_chunks))])

# Fill the matrix with 1s for each target in its respective chunk
for chunk_idx, chunk in enumerate(target_chunks):
    target_chunk_matrix.loc[chunk, f'chunk_{chunk_idx}'] = 1

target_chunk_matrix.to_csv(f'{datadir}/{target2chunk_file}', compression='gzip')