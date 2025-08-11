import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sys import exit
import re
from scipy.stats import pearsonr

def parse_filename_de(filename):
	match = re.match(r"DE_(.*?)\.csv", filename)
	if match:
		var = match.groups()
		return var[0]
	return None

def parse_filename_perturb(filename):
	match = re.match(r"DE_(Resting|Stimulated)-(Teff|Treg)_(?:target_)?(\w+).csv", filename)
	if match:
		condition, celltype, target = match.groups()
		return condition, celltype, target
	return None, None, None

def load_data_perturb(directory):
	data = {}
	for filename in os.listdir(directory):
		if filename.endswith(".csv"):
			condition, celltype, target = parse_filename_perturb(filename)
			if condition and celltype and target:
				filepath = os.path.join(directory, filename)
				df = pd.read_csv(filepath, index_col=0)
				
				if '-'.join([condition,celltype]) not in data:
					data['-'.join([condition,celltype])] = {}
				data['-'.join([condition,celltype])][target] = df
	return data
	
def load_data_de(directory):
	data = {}
	for filename in os.listdir(directory):
		if filename.endswith(".csv"):
			var = parse_filename_de(filename)
			if var:
				filepath = os.path.join(directory, filename)
				df = pd.read_csv(filepath, index_col=0)
				
				if var not in data:
					data[var] = {}
				data[var] = df
	return data

def plot_correlation(data, condition1, target1, condition2, target2, figdir, plot='automatic'):
	df1 = data[condition1][target1].copy()
	df2 = data[condition2][target2].copy()
	
	df1['z'] = df1['log2FoldChange'] / df1['lfcSE']
	df2['z'] = df2['log2FoldChange'] / df2['lfcSE']
	
	# Drop NaN and infinite values
	df1 = df1.replace([np.inf, -np.inf], np.nan).dropna(subset=['z'])
	df2 = df2.replace([np.inf, -np.inf], np.nan).dropna(subset=['z'])

	case1 = f"{condition1}_{target1}"
	case2 = f"{condition2}_{target2}"

	# Merge on gene identifiers (index)
	merged_df = df1[['z']].rename(columns={'z': case1}).merge(
		df2[['z']].rename(columns={'z': case2}),
		left_index=True, right_index=True
	).dropna()  # Drop any remaining NaNs after merging

	# Ensure there are valid values for correlation
	if len(merged_df) < 2 or merged_df[case1].std() == 0 or merged_df[case2].std() == 0:
		print("Insufficient data for correlation calculation.")
		return np.nan, np.nan

	# Compute correlation and p-value
	correlation, p_value = pearsonr(merged_df[case1], merged_df[case2])

	if plot == True or (plot == 'automatic' and np.abs(correlation) > 0.15):
		# Compute quadrants
		q1 = ((merged_df[case1] > 0) & (merged_df[case2] > 0)).sum()
		q2 = ((merged_df[case1] < 0) & (merged_df[case2] > 0)).sum()
		q3 = ((merged_df[case1] < 0) & (merged_df[case2] < 0)).sum()
		q4 = ((merged_df[case1] > 0) & (merged_df[case2] < 0)).sum()
		total = len(merged_df)

		# Percentages
		q1_pct, q2_pct, q3_pct, q4_pct = (q1 / total * 100, q2 / total * 100, q3 / total * 100, q4 / total * 100)

		# Plot
		plt.figure(figsize=(6, 6))
		plt.axhline(0, color='gray', linestyle='--')
		plt.axvline(0, color='gray', linestyle='--')
		sns.scatterplot(data=merged_df, x=case1, y=case2, alpha=0.7, s=4)
		sns.regplot(data=merged_df, x=case1, y=case2, scatter=False, color='red', line_kws={'linewidth': 1})  # Best fit line

		# Label quadrants
		plt.text(0.8, 0.9, f"{q1_pct:.1f}%", transform=plt.gca().transAxes, fontsize=9, color='k')
		plt.text(0.1, 0.9, f"{q2_pct:.1f}%", transform=plt.gca().transAxes, fontsize=9, color='k')
		plt.text(0.1, 0.1, f"{q3_pct:.1f}%", transform=plt.gca().transAxes, fontsize=9, color='k')
		plt.text(0.8, 0.1, f"{q4_pct:.1f}%", transform=plt.gca().transAxes, fontsize=9, color='k')

		plt.title(f"{condition1} {target1} and {condition2} {target2}\nCorrelation: {correlation:.2f}, p-value: {p_value:.2e}")
		plt.xlabel(f"{case1} Fold Change z-score")
		plt.ylabel(f"{case2} Fold Change z-score")
		plt.savefig(f"{figdir}correlation_{case1}_{case2}.png", dpi=400)

	return correlation, p_value

#def heatmap(correlations, p_values):
#	correlations = correlations.fillna(0).astype(float)
#	p_values = p_values.fillna(1).astype(float) 
#	
#	# Create annotation labels with stars for significance
#	stars = np.where(p_values < 0.01, '*', '')  # Only annotate with stars
#	
#	fig, ax = plt.subplots(figsize=(12, 10))
#	g = sns.clustermap(
#		correlations, cmap="bwr", annot=stars, fmt="", linewidths=0.5, 
#		center=0, row_cluster=True, col_cluster=True
#	)
#	return fig, g

def heatmap(correlations, p_values, match_order=False):
	correlations = correlations.fillna(0).astype(float)
	p_values = p_values.fillna(1).astype(float)
	
	# Create annotation labels with stars for significance
	stars = np.where(p_values < 0.01, '*', '')  # Only annotate with stars
	
	fig, ax = plt.subplots(figsize=(12, 10))
	g = sns.clustermap(
		correlations, cmap="bwr", annot=stars, fmt="", linewidths=0.5, 
		center=0, row_cluster=True, col_cluster=True, annot_kws={"size": 13}
	)
	
	# Match clustering order if match_order is True
	if match_order:
		order = g.dendrogram_row.reordered_ind  # Ensure row and column order are the same
		correlations = correlations.iloc[order, order]
		stars = stars[order][:, order]
		
		plt.close(g.fig)  # Close previous figure to avoid duplication
		g = sns.clustermap(
			correlations, cmap="bwr", annot=stars, fmt="", linewidths=0.5,
			center=0, row_cluster=False, col_cluster=False, annot_kws={"size": 13}
		)
	
	return fig, g


def heatmap_coefs(df):
	fig, ax = plt.subplots(figsize=(10, 8))
	g = sns.clustermap(df, cmap="bwr", annot=False, linewidths=0.5, center=0, row_cluster=True, col_cluster=False)
	return fig, g

def create_chromosome_dict(gene_list):
	from pyensembl import EnsemblRelease
	gene_chromosome_dict = {}

	# Load human Ensembl dataset (GRCh38)
	data = EnsemblRelease(release=104, species='homo_sapiens')
	
	for gene in gene_list:
		try:
			chrom = data.gene_by_id(gene).contig
			previous_chrom = chrom
			gene_chromosome_dict[gene] = chrom
		except:
			gene_chromosome_dict[gene] = previous_chrom  # they seem to be in order

	return gene_chromosome_dict

def create_gene_name_dict(gene_list):
	from pyensembl import EnsemblRelease
	gene_name_dict = {}

	# Load human Ensembl dataset (GRCh38)
	ensembl = EnsemblRelease(release=104, species='homo_sapiens')
	
	for gene in gene_list:
		name = ensembl.gene_by_id(gene).gene_name
		previous_chrom = name
		gene_name_dict[gene] = name

	return gene_name_dict



def plot_pearson_correlations(corr_df, se_df):
	"""
	Plots a bar chart of Pearson correlations with standard errors.
	
	Parameters:
	corr_df (pd.Series): A pandas Series containing correlation values.
	se_df (pd.Series): A pandas Series containing standard errors.
	"""
	# Combine data into a DataFrame and sort by correlation values
	df = pd.DataFrame({'correlation': corr_df, 'std_error': se_df})
	df = df.sort_values(by='correlation', ascending=False)
	
	# Define colors: highlight 'all_targets' in red
	colors = ['red' if gene == 'all_targets' else 'royalblue' for gene in df.index]
	
	# Plot settings
	fig, ax = plt.subplots(figsize=(10,6))
	ax.bar(df.index, df['correlation'], yerr=df['std_error'], color=colors, capsize=5)
	
	# Formatting
	plt.xticks(rotation=45, ha='right')
	ax.set_ylabel('Pearson Correlation')
	return fig
