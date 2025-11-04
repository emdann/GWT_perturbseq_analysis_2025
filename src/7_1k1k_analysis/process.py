import anndata
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix
from adjustText import adjust_text

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def load_data(h5ad_file):
	print('Loading dataset')
	return anndata.read_h5ad(h5ad_file)

def filter_cell_types(adata, keep_cell_types):
	CD4_pos = adata.obs['cell_type'].isin(keep_cell_types)
	return adata[CD4_pos]

def normalize_counts(adata):
	print('Normalizing counts')
	scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
	adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
	return adata

def plot_normalization(adata, figdir, counts_col='nCount_RNA'):
	plt.clf()
	fig, axes = plt.subplots(1, 2, figsize=(10, 5))
	sns.histplot(adata.obs[counts_col], bins=100, kde=False, ax=axes[0])
	axes[0].set_title("Total counts")
	sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
	axes[1].set_title("Shifted logarithm")
	# plt.savefig(figdir+'shifted_log_normalization.pdf')

def perform_qc(adata, figdir, gene_name_col='feature_name', counts_col='nCount_RNA'):
	print('Performing QC')
	adata.var["mt"] = adata.var[gene_name_col].str.startswith("MT-")
	adata.var["ribo"] = adata.var[gene_name_col].str.startswith(("RPS", "RPL"))
	adata.var["hb"] = adata.var[gene_name_col].str.contains("^HB[^(P)]")
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
	plt.clf()
	fig, axs = plt.subplots(2, 2, figsize=(15, 12))
	sns.histplot(adata.obs[counts_col], bins=100, kde=False, ax=axs[0,0])
	sc.pl.violin(adata, "pct_counts_mt", ax=axs[1,0], show=False)
	sc.pl.scatter(adata, counts_col, "n_genes_by_counts", color="pct_counts_mt", ax=axs[0,1], show=False)
	plt.tight_layout()
	# plt.savefig(figdir+'qc.png', dpi=300)
	return adata

def is_outlier(adata, metric: str, nmads: int):
	M = adata.obs[metric]
	return (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)

def filter_outliers(adata):
	adata.obs["outlier"] = (
		is_outlier(adata, "log1p_total_counts", 5) |
		is_outlier(adata, "log1p_n_genes_by_counts", 5) |
		is_outlier(adata, "pct_counts_in_top_20_genes", 5)
	)
	adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)
	print(f"Total number of cells: {adata.n_obs}")
	adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
	print(f"Number of cells after filtering of low-quality cells: {adata.n_obs}")
	return adata

def select_highly_variable_genes(adata, figdir):
	print('Selecting highly variable genes')
	sc.pp.highly_variable_genes(adata, layer="log1p_norm", n_top_genes=5000)
	plt.clf()
	ax = sns.scatterplot(data=adata.var, x="means", y="dispersions", hue='highly_variable', s=5)
	ax.set_xlim(None, 1.5)
	ax.set_ylim(None, 3)
	# plt.savefig(figdir+'variable_genes.pdf')

def perform_dimensionality_reduction(adata, figdir, color_by='age'):
	print('Performing dimensionality reduction')
	adata.layers["counts"] = adata.X.copy()
	adata.X = adata.layers['log1p_norm'].copy()
	sc.pp.pca(adata, svd_solver="arpack", mask_var='highly_variable')

	plt.clf()
	sc.pl.pca_scatter(adata, color=color_by)
	# plt.savefig(figdir+'pca.png', dpi=300)

	sc.pp.neighbors(adata, n_pcs=40)
	sc.tl.umap(adata)
	sc.pl.umap(adata, color=color_by)
	plt.tight_layout()
	# plt.savefig(figdir+'umap_'+color_by+'.png', dpi=300)
	return adata

def save_processed_data(adata, output_file):
	adata.write(output_file, compression="gzip")
	print('Wrote processed data')

def read_total_counts(h5ad_file, datadir):
	adata_full = anndata.read_h5ad(h5ad_file, backed='r')
	
	donor_ids_full = adata_full.obs['donor_id']
	cell_types_full = adata_full.obs['predicted.celltype.l2']
	total_counts_full = adata_full.obs['nCount_RNA']
	
	df = pd.DataFrame({'donor_id': donor_ids_full, 'cell_type': cell_types_full, 'total_counts': total_counts_full})
	
	cell_type_counts_full = pd.crosstab(df['donor_id'], df['cell_type'])
	total_counts_sum_full = pd.crosstab(df['donor_id'], df['cell_type'], values=df['total_counts'], aggfunc=sum)
	
	cell_type_counts_full.to_csv(datadir+'cell_type_counts_full.csv', index=True)
	total_counts_sum_full.to_csv(datadir+'total_counts_sum_full.csv', index=True)

	return cell_type_counts_full, total_counts_sum_full

def compute_cell_type_counts(adata):
	df = pd.DataFrame({
		'donor_id': adata.obs['donor_id'], 
		'cell_type': adata.obs['predicted.celltype.l2'], 
		'total_counts': adata.obs['total_counts']
	})
	cell_type_counts = pd.crosstab(df['donor_id'], df['cell_type'])
	cell_type_counts['Total'] = cell_type_counts.sum(axis=1)
	total_counts_sum = pd.crosstab(df['donor_id'], df['cell_type'], values=df['total_counts'], aggfunc=sum)
	total_counts_sum['Total'] = total_counts_sum.sum(axis=1)
	return cell_type_counts, total_counts_sum


def gene_expression_sum_mean_var(adata, donor_col='donor_id', celltype_col='predicted.celltype.l2', additional_col=None):
	###### Group by donor, cell type, and optionally additional column
	if additional_col:
		group_keys = adata.obs[donor_col].str.cat([adata.obs[celltype_col], adata.obs[additional_col]], sep=', ')
	else:
		group_keys = adata.obs[donor_col].str.cat(adata.obs[celltype_col], sep=', ')
	
	unique_groups, group_indices = np.unique(group_keys, return_inverse=True)
	
	hvg_sparse_matrix = adata[:, adata.var['highly_variable']].X  # Sparse matrix of HVGs
	
	n_groups = len(unique_groups)
	n_samples = adata.shape[0]
	group_indicator = csr_matrix((np.ones(n_samples), (group_indices, np.arange(n_samples))), shape=(n_groups, n_samples))
	
	# Perform sparse matrix multiplication to sum rows by group
	counts_sum_sparse = group_indicator @ hvg_sparse_matrix  # Sparse matrix with grouped sums
	
	# Convert the result to a DataFrame
	counts_sum = pd.DataFrame.sparse.from_spmatrix(
		counts_sum_sparse,
		index=unique_groups,
		columns=adata.var.index[adata.var['highly_variable']]
	)
	
	### Counts mean and variance
	n_cells_per_group = group_indicator.sum(axis=1).A1
	gene_means = counts_sum_sparse / n_cells_per_group[:, None]
	
	gene_means_df = pd.DataFrame.sparse.from_spmatrix(
		gene_means,
		index=unique_groups,
		columns=adata.var.index[adata.var['highly_variable']]
	)
	
	gene_means_df['cell type'] = gene_means_df.index.str.split(', ').str[1]
	
	cell_types_list = gene_means_df['cell type'].unique()
	mean_var_expression_donor = pd.DataFrame(index=gene_means_df.columns)
	for ct in cell_types_list:
		means = gene_means_df.loc[gene_means_df['cell type'] == ct].drop(columns=['cell type']).astype(float)
		mean_var_expression_donor['Mean, ' + ct] = means.mean(axis=0)
		mean_var_expression_donor['Variance, ' + ct] = means.var(axis=0)
	mean_var_expression_donor = mean_var_expression_donor.drop('cell type')
	
	return counts_sum_sparse, counts_sum, mean_var_expression_donor


def calculate_metadata_by_donor_celltype(adata):
	adata.obs['male'] = adata.obs['sex']=='male'
	
	# Create binary columns for each cell type level except the last one
	cell_types = adata.obs['predicted.celltype.l2'].unique()
	for ct in cell_types[:-1]:
		col_name = ct.replace(' ', '_')
		adata.obs[col_name] = adata.obs['predicted.celltype.l2']==ct
	
	# Get list of binary cell type columns
	cell_type_cols = [ct.replace(' ', '_') for ct in cell_types[:-1]]
	
	# Create metadata matrix with all relevant columns
	metadata_cols = ['donor_id', 'predicted.celltype.l2', 'male', 'age', 'total_counts', 'pool_number'] + cell_type_cols
	metadata_matrix = adata.obs[metadata_cols]
	metadata_matrix['cell_counts'] = 1
	
	metadata_matrix.index = metadata_matrix['donor_id'].str.cat(metadata_matrix['predicted.celltype.l2'], sep=', ')
	
	metadata = metadata_matrix.groupby(metadata_matrix.index).first().drop(columns=['total_counts','cell_counts'])
	metadata[['total_counts','cell_counts']] = metadata_matrix[['total_counts','cell_counts']].groupby(metadata_matrix.index).sum()[['total_counts','cell_counts']]
	metadata.index.rename('donor+cell_type', inplace=True)
	metadata['avg_count_per_cell'] = metadata['total_counts'] / metadata['cell_counts']
	metadata['log_counts_per_cell'] = np.log10(metadata['avg_count_per_cell'])
	metadata["age_cat"] = metadata["age"].apply(categorize_age)

	pool_one_hot = pd.get_dummies(metadata['pool_number'], prefix='pool')
	metadata = pd.concat([metadata, pool_one_hot], axis=1)
	return metadata

def calculate_metadata_by_donor(adata, cell_type_counts_full, total_counts_sum_full, B_T_pct_df):
	adata.obs['male'] = adata.obs['sex']=='male'
	metadata_donors = adata.obs.groupby(adata.obs['donor_id']).first()[['male','age', 'pool_number']]
	metadata_donors["age_cat"] = metadata_donors["age"].apply(categorize_age)
	pool_one_hot = pd.get_dummies(metadata_donors['pool_number'], prefix='pool')
	metadata = pd.concat([metadata_donors, pool_one_hot], axis=1)
	
	metadata_donors["total_cells"] = cell_type_counts_full.sum(axis=1)
	metadata_donors["total_counts"] = total_counts_sum_full.sum(axis=1)
	metadata_donors['avg_count_per_cell'] = metadata_donors['total_counts'] / metadata_donors['total_cells']
	metadata_donors['log_counts_per_cell'] = np.log10(metadata_donors['avg_count_per_cell'])
	metadata_donors = metadata_donors.reset_index().merge(B_T_pct_df, on="donor_id", how="left").set_index('donor_id')

	return metadata_donors

def calculate_metadata_by_target_condition(adata):
	adata.obs['Stimulated'] = adata.obs['condition']=='Stimulated'
	adata.obs['Teff'] = adata.obs['cell_type']=='Teff'
	adata.obs['individualA'] = adata.obs['individual']=='A'
	
	metadata_matrix = adata.obs[['sgrna','target','hash_ID','Stimulated','Teff','individualA','S_Score','G2M_Score','total_counts','individual']]
	metadata_matrix['cell_counts'] = 1
	target_dummies = pd.get_dummies(metadata_matrix['target'], prefix='target').astype(int)
	target_dummies.columns = target_dummies.columns.str.replace('-', '_')
	target_dummies = target_dummies.sub(target_dummies['target_NO_TARGET'], axis=0)
	metadata_matrix = pd.concat([metadata_matrix, target_dummies], axis=1) 
	
	metadata_matrix.index = metadata_matrix['sgrna'].str.cat(metadata_matrix['hash_ID'], sep=', ').str.cat(metadata_matrix['individual'], sep=', ')
	
	metadata = metadata_matrix.groupby(metadata_matrix.index).first().drop(columns=['total_counts','cell_counts'])
	metadata[['total_counts','cell_counts']] = metadata_matrix[['total_counts','cell_counts']].groupby(metadata_matrix.index).sum()[['total_counts','cell_counts']]
	metadata.index.rename('target+condition', inplace=True)
	metadata['avg_count_per_cell'] = metadata['total_counts'] / metadata['cell_counts']
	metadata['log_counts_per_cell'] = np.log10(metadata['avg_count_per_cell'])

	return metadata


def calculate_B_T_pct(cell_type_counts_full):
	# calculate B and T percentage
	b_cell_types = ["B intermediate", "B memory", "B naive"]
	t_cell_types = ["CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "Treg", "dnT", "gdT", "MAIT"]
	
	total_cells_per_donor = cell_type_counts_full.sum(axis=1)
	total_b_cells = cell_type_counts_full[b_cell_types].sum(axis=1)
	total_t_cells = cell_type_counts_full[t_cell_types].sum(axis=1)
	
	percent_b_cells = (total_b_cells / total_cells_per_donor) * 100
	percent_t_cells = (total_t_cells / total_cells_per_donor) * 100
	
	percentages_df = pd.DataFrame({"percent_B_cells": percent_b_cells, "percent_T_cells": percent_t_cells})
	percentages_df = percentages_df.reset_index()

	return percentages_df

def dimensionality_reduction_summed(adata_sums, figdir, plot=False):
	scales_counts = sc.pp.normalize_total(adata_sums, target_sum=None, inplace=False)
	adata_sums.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
	sc.pp.pca(adata_sums, layer='log1p_norm') 

	if plot:
		plt.close()
		for col in adata_sums.obs.columns:
			plt.clf()
			sc.pl.pca_scatter(adata_sums, color=col, annotate_var_explained=True, components=['1,2', '2,3', '3,4','4,5'], size=20000/len(adata_sums))
			#sc.pl.pca_scatter(adata_sums, color=col, annotate_var_explained=True, components=['2,3'], size=35000/len(adata_sums))
			#plt.title('Cell Type / State')
			plt.tight_layout()
			# plt.savefig(figdir+'umaps_on_summed/pca_'+col+'.png', dpi=300)

		# UMAP
		sc.pp.neighbors(adata_sums, n_pcs=20)
		sc.tl.umap(adata_sums)
		
		for col in adata_sums.obs.columns:
			plt.clf()
			sc.pl.umap(adata_sums, color=col)
			plt.tight_layout()
			# plt.savefig(figdir+'umaps_on_summed/umap_'+col+'.png', dpi=300)
		
	return adata_sums
	
	## t-SNE
	#sc.tl.tsne(adata_sums, use_rep="X_pca", n_pcs=20)
	#
	#plt.close()
	#for col in adata_sums.obs.columns:
	#	plt.clf()
	#	sc.pl.tsne(adata_sums, color=col)
	#	plt.tight_layout()
	#	plt.savefig(figdir+'umaps_on_summed/tsne_'+col+'.png', dpi=300)
	#
	#return adata_sums

def pool_differences_PCA(adata_sums, chromosome_dict, figdir):
	if 'predicted.celltype.l2' in adata_sums.obs:
		adata_pcs = adata_sums[adata_sums.obs['predicted.celltype.l2']=='CD4 Naive']
	else:
		adata_pcs = adata_sums.copy()
	adata_pcs.var['chromosome'] = adata_pcs.var.index.map(chromosome_dict)
	adata_pcs = adata_pcs[:, ~adata_pcs.var["chromosome"].isin(["X", "Y", "MT"])].copy()

	scales_counts = sc.pp.normalize_total(adata_pcs, target_sum=None, inplace=False)
	adata_pcs.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
	sc.pp.pca(adata_pcs, layer='log1p_norm') 

	plt.close()
	for col in adata_pcs.obs.columns:
		plt.clf()
		sc.pl.pca_scatter(adata_pcs, color=col, annotate_var_explained=True, components=['1,2', '2,3', '3,4','4,5'])
		plt.tight_layout()
		# plt.savefig(figdir+'umaps_covars/pca_'+col+'.png', dpi=300)
	
	return adata_pcs


def fit_DEseql(metadata, counts_sum, cols=['male', 'CD4_Naive', 'CD4_TCM', 'age_cat', 'percent_B_cells', 'percent_T_cells']):
	covars = metadata[cols]
	
	design = '~' + ' + '.join(covars.columns)
	dds = DeseqDataSet(counts=counts_sum, metadata=covars, design=design)
	
	dds.fit_size_factors()
	dds.fit_genewise_dispersions()
	dds.fit_dispersion_trend()
	dds.fit_dispersion_prior()
	print(
		f"logres_prior={dds.uns['_squared_logres']}, sigma_prior={dds.uns['prior_disp_var']}"
	) 
	dds.fit_MAP_dispersions()
	dds.fit_LFC()
	return dds

def categorize_age(age):
	if age < 40:
		return 1
	elif 40 <= age < 60:
		return 2
	elif 60 <= age < 70:
		return 3
	elif 70 <= age < 80:
		return 4
	else:
		return 5


def create_chromosome_dict(gene_list):
	from pyensembl import EnsemblRelease
	gene_chromosome_dict = {}

	# Load human Ensembl dataset (GRCh38)
	data = EnsemblRelease(release=104, species='homo_sapiens')
	
	for gene in gene_list:
		# try:
		chrom = data.gene_by_id(gene).contig
		previous_chrom = chrom
		gene_chromosome_dict[gene] = chrom
		# except:
		# 	gene_chromosome_dict[gene] = previous_chrom  # they seem to be in order

	return gene_chromosome_dict


def compute_target_counts(adata):
	df = pd.DataFrame({
		'target': adata.obs['target'], 
		'condition': adata.obs['hash_ID'], 
		'total_counts': adata.obs['total_counts']
	})
	cell_type_counts = pd.crosstab(df['target'], df['condition'])
	cell_type_counts['Total'] = cell_type_counts.sum(axis=1)
	total_counts_sum = pd.crosstab(df['target'], df['condition'], values=df['total_counts'], aggfunc=sum)
	total_counts_sum['Total'] = total_counts_sum.sum(axis=1)
	return cell_type_counts, total_counts_sum

