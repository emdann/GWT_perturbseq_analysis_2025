import time
start = time.time()
import numpy as np
nptime = time.time()
print('numpy:',nptime-start)
import pandas as pd
pdtime = time.time()
print('pandas:',pdtime-nptime)
import matplotlib.pyplot as plt
plttime = time.time()
print('plt:',plttime-pdtime)
import anndata 
anntime = time.time()
print('anndata:',anntime-plttime)
import scanpy as sc
sctime = time.time()
print('scanpy:',sctime-anntime)
import seaborn as sns
snstime = time.time()
print('sns:',snstime-sctime)
import os
from scipy.stats import median_abs_deviation
from sys import exit
from pydeseq2.ds import DeseqStats
import pickle

import process as process
import plot as plot

def run_deseq_analysis(dds, covar, output_suffix, figdir, resultsdir):
    contrast = np.zeros(len(dds.obsm['design_matrix'].columns), dtype=int)
    contrast[np.where(dds.obsm['design_matrix'].columns == covar)[0][0]] = 1

    dds.refit()
    ds = DeseqStats(
        dds,
        contrast=contrast,
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
    )
    ds.run_wald_test()
    if ds.cooks_filter:
        ds._cooks_filtering()
    
    if ds.independent_filter:
        ds._independent_filtering()
    else:
        ds._p_value_adjustment()
    
    ds.summary()
    ds.results_df.to_csv(resultsdir+'DE_'+covar+f'_{output_suffix}.csv')
    results_df = ds.results_df

    pvalues = results_df["pvalue"].dropna()
    plot.plot_pval_distribution(pvalues, covar+f'_{output_suffix}', figdir)
    
    return results_df

h5ad_file = '/mnt/oak/users/emma/data/cxg_datasets/Yazar2022.h5ad'
datadir = '/mnt/oak/users/emma/data/GWT/OneK1K_analysis/'
figdir = 'figures/'
resultsdir = 'results/'

load_processed_data = True
ct = 'CD8T'

h5ad_file_processed = datadir + f"Yazar2022_{ct}_processed.h5ad"

print('loading processed 1k1k dataset')
adata_1k1k = anndata.experimental.read_lazy(h5ad_file_processed)
adata = anndata.AnnData(
        obs=adata_1k1k.obs.to_dataframe(),
        var=adata_1k1k.var.to_dataframe(),
        X = adata_1k1k.X
    )
adata = adata.to_memory()
adata.var['highly_variable'] = adata.var_names.isin(adata.var['dispersions_norm'].nlargest(10000).index)

print('counting cells')
cell_type_counts_full, total_counts_sum_full = process.read_total_counts(h5ad_file, datadir)

adata.X.data = adata.X.data.astype(float)

# process by cell type
cell_type_counts, total_counts_sum = process.compute_cell_type_counts(adata)

counts_sum_sparse, counts_sum, mean_var_expression_donor = process.gene_expression_sum_mean_var(adata)
metadata = process.calculate_metadata_by_donor_celltype(adata)
B_T_pct_df = process.calculate_B_T_pct(cell_type_counts_full)
metadata = metadata.reset_index().merge(B_T_pct_df, on="donor_id", how="left").set_index('donor+cell_type')

metadata_donors = process.calculate_metadata_by_donor(adata, cell_type_counts, cell_type_counts_full, total_counts_sum_full, B_T_pct_df)

adata_sums = anndata.AnnData(
	X=counts_sum_sparse,
	obs=metadata,	
	var=adata.var.loc[counts_sum.columns][['feature_name']]  
)

adata_sums = process.dimensionality_reduction_summed(adata_sums, figdir)
adata_sums.write_h5ad(datadir + f"Yazar2022_{ct}_processed.pbulk.h5ad")
gene_names_dict = dict(zip(adata.var.index, adata.var['feature_name']))
chromosome_dict = process.create_chromosome_dict(adata.var.index)

# add top 5 PCs as covars to account for pool differences
if ct == 'CD4T':
    adata_pcs = process.pool_differences_PCA(adata_sums, chromosome_dict, figdir)
    pcs_df = pd.DataFrame(adata_pcs.obsm["X_pca"][:, :5], index=adata_pcs.obs.index, columns=[f"PC{i+1}" for i in range(5)])
    pcs_df["donor_id"] = adata_pcs.obs["donor_id"]
else:
    metadata_df = pd.read_csv(datadir + 'Yazar2022_CD4T_processed.metadata.csv')
    pcs_df = metadata_df[['donor_id'] + [f'PC{i}' for i in range(1,6)]].drop_duplicates()

metadata = metadata.merge(pcs_df, on="donor_id", how="left").set_index(metadata.index)

# DEseqs
print('testing')
metadata = metadata.loc[:,~metadata.columns.str.startswith('pool')]
metadata.to_csv(datadir + f'Yazar2022_{ct}_processed.metadata.csv')

# Split donors into train (80%) and validation (20%) sets with stratification by age_cat
unique_donors = metadata['donor_id'].unique()
donor_age_cats = metadata.groupby('donor_id')['age_cat'].first()

train_donors = []
val_donors = []

for age_cat in donor_age_cats.unique():
    age_cat_donors = donor_age_cats[donor_age_cats == age_cat].index
    n_train = int(0.8 * len(age_cat_donors))
    train_age_cat = np.random.choice(age_cat_donors, size=n_train, replace=False)
    val_age_cat = np.array([d for d in age_cat_donors if d not in train_age_cat])
    
    train_donors.extend(train_age_cat)
    val_donors.extend(val_age_cat)

train_donors = np.array(train_donors)
val_donors = np.array(val_donors)

# Create train and validation metadata/counts
train_metadata = metadata[metadata['donor_id'].isin(train_donors)]
val_metadata = metadata[metadata['donor_id'].isin(val_donors)]
train_counts = counts_sum.loc[train_metadata.index]
val_counts = counts_sum.loc[val_metadata.index]

covars = list(metadata.columns.drop(['age','donor_id','predicted.celltype.l2','total_counts','cell_counts','avg_count_per_cell','log_counts_per_cell']))

# Fit DESeq for training and validation sets
train_dds = process.fit_DEseql(train_metadata, train_counts, cols=covars)
val_dds = process.fit_DEseql(val_metadata, val_counts, cols=covars)

for covar in covars:
    if covar.startswith('PC') or covar.startswith('CD'): continue
    print(covar)
    
    train_results = run_deseq_analysis(train_dds, covar, f'{ct}_train', figdir, resultsdir)
    val_results = run_deseq_analysis(val_dds, covar, f'{ct}_validation', figdir, resultsdir)
