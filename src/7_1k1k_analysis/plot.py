import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
import seaborn as sns
import pandas as pd
import anndata 
import scanpy as sc
import matplotlib.patches as mpatches
from sklearn.preprocessing import MinMaxScaler
from adjustText import adjust_text

def plot_umap_by_obs(adata, figdir):
	for col in adata.obs.columns:
		if len(adata.obs[col].unique())>1:
			plt.clf()
			if col=='sex':
				sc.pl.umap(adata, color=col, title='Sex', palette={'female': 'violet', 'male': 'blue'})
			elif col=='predicted.celltype.l2':
				sc.pl.umap(adata, color=col, title='Cell Type')
			else: continue
			plt.tight_layout()
			plt.savefig(figdir+'umaps/umap_'+col+'.png', dpi=300)


def plot_cell_type_counts(cell_type_counts, figdir):
	fig, axes = plt.subplots(nrows=1, ncols=cell_type_counts.shape[1], figsize=(15, 5), sharey=True)
	for i, cell_type in enumerate(cell_type_counts.columns):
		axes[i].hist(cell_type_counts[cell_type], bins=np.arange(0,cell_type_counts.max().median(),100), color='skyblue', edgecolor='black')
		axes[i].set_title(f'{cell_type}')
		axes[i].set_xlabel('Number of cells')
		axes[i].set_ylabel('Number of Donors/Perturbations')
	plt.tight_layout()
	plt.savefig(figdir+'cell_type_counts_each_donor.pdf')

def analyze_total_counts(metadata_donors, figdir):
	metadata_donors["CD4 Naive %"] = (metadata_donors["CD4 Naive"] / metadata_donors["Total"]) * 100
	metadata_donors["CD4 TCM %"] = (metadata_donors["CD4 TCM"] / metadata_donors["Total"]) * 100
	plt.figure(figsize=(8, 6))
	plt.scatter(metadata_donors["CD4 Naive %"], metadata_donors["CD4 TCM %"], alpha=0.7)
	plt.xlabel("CD4 Naive (%)", color='darkviolet')
	plt.ylabel("CD4 TCM (%)", color='darkturquoise')
	plt.title("CD4 Naive vs CD4 TCM Percentages")
	plt.grid(True)
	plt.savefig(figdir+'naive_tcm_percentages.pdf')

def calculate_r_squared(x, y, coeffs):
	p = np.poly1d(coeffs)
	y_pred = p(x)
	ss_res = np.sum((y - y_pred) ** 2)
	ss_tot = np.sum((y - np.mean(y)) ** 2)
	return 1 - (ss_res / ss_tot)

def cell_counts_boxplot(cell_type_counts_full, figdir):
	# Compute median count per cell type and sort in descending order
	medians = cell_type_counts_full.median().sort_values(ascending=False)

	# Select only valid cell types
	plot_data = cell_type_counts_full[medians[medians>1].index]
	plot_data.columns = plot_data.columns.set_categories(list(plot_data.columns), ordered=True)

	# Explicitly pass sorted order to sns.boxplot()
	plt.figure(figsize=(12, 6))
	sns.boxplot(data=plot_data, palette='crest') 
	plt.xticks(fontsize=12, rotation=60)
	plt.xlabel("Cell Type", fontsize=15)
	plt.ylabel("Number of Cells per Donor", fontsize=15)
	plt.title("Number of Cells per Donor Across Cell Types", fontsize=20)
	plt.tight_layout()
	plt.savefig(figdir + 'total_cells_each_donor.pdf')

def mean_vs_variance_expression(mean_var_expression_donor, figdir):
	# Identify cell types dynamically
	cell_types = set()
	for col in mean_var_expression_donor.columns:
		if col.startswith("Mean, "):
			cell_types.add(col.replace("Mean, ", ""))
	
	# Assign colors
	colors = ["blue", "green", "red", "purple"]  # Extend if necessary
	color_map = {cell_type: colors[i % len(colors)] for i, cell_type in enumerate(cell_types)}
	
	# Plot log-log scatter
	plt.clf()
	plt.figure(figsize=(8, 6))
	for cell_type in cell_types:
		mean_col = f"Mean, {cell_type}"
		var_col = f"Variance, {cell_type}"
		
		if mean_col in mean_var_expression_donor.columns and var_col in mean_var_expression_donor.columns:
			mean_values = mean_var_expression_donor[mean_col]
			var_values = mean_var_expression_donor[var_col]
			
			# Filter out zero values to avoid log(0) issues
			nonzero_mask = (mean_values > 0) & (var_values > 0)
			mean_values = mean_values[nonzero_mask]
			var_values = var_values[nonzero_mask]
			
			plt.scatter(np.log10(mean_values), np.log10(var_values), 
						label=f'one gene, {cell_type}', 
						alpha=0.7, color=color_map[cell_type], s=1)
	
	plt.title("Mean vs Variance by Cell Type", fontsize=14)
	plt.xlabel("Log10(Mean(each gene))", fontsize=12)
	plt.ylabel("Log10(Variance(each gene))", fontsize=12)
	plt.legend(title="Cell Type")
	plt.grid(True, which="both", linestyle="--", linewidth=0.5)
	plt.savefig(figdir + 'mean_vs_variance_each_gene.png', dpi=300)


#def mean_vs_variance_expression(mean_var_expression_donor, figdir):
#	# Plot log-log scatter
#	plt.clf()
#	colors = ["blue", "green", "red"]  # Assign a color to each cell type
#	plt.figure(figsize=(8, 6))
#	for cell_type, color in zip(['CD4 Naive', 'CD4 TCM', 'CD4 TEM'], colors):
#		mean_col = f"Mean, {cell_type}"
#		var_col = f"Variance, {cell_type}"
#		mean_values = mean_var_expression_donor[mean_col]
#		var_values = mean_var_expression_donor[var_col]
#		
#		# Filter out zero values to avoid log(0) issues
#		nonzero_mask = (mean_values > 0) & (var_values > 0)
#		mean_values = mean_values[nonzero_mask]
#		var_values = var_values[nonzero_mask]
#		
#		plt.scatter(np.log10(mean_values), np.log10(var_values), label='one gene, '+cell_type, alpha=0.7, color=color, s=1)
#	
#	plt.title("Mean vs Variance Across Donors by Cell Type", fontsize=14)
#	plt.xlabel("Log10(Mean(each gene,across donors))", fontsize=12)
#	plt.ylabel("Log10(Variance(each gene,across donors))", fontsize=12)
#	plt.legend(title="Cell Type")
#	plt.grid(True, which="both", linestyle="--", linewidth=0.5)
#	plt.savefig(figdir+'mean_vs_variance_across_donors_each_gene.png', dpi=300)

def total_cells_vs_total_counts(metadata_donors, figdir):
	plt.figure(figsize=(8, 6))
	plt.scatter(metadata_donors["total_cells"], metadata_donors["total_counts"], 
				alpha=0.7, s=3, c=metadata_donors["color"], label=metadata_donors["sex"])
	plt.yscale('log')
	plt.xscale('log')
	plt.title("Total Counts vs Total Cells for Each Donor", fontsize=14)
	plt.xlabel("Log10(Total Cells)", fontsize=12)
	plt.ylabel("Log10(Total Counts)", fontsize=12)
	plt.grid(True, which="both", linestyle="--", linewidth=0.5)
	legend_elements = [Line2D([0], [0], marker='o', color='w', markersize=5, markerfacecolor='red', label='Female'), Line2D([0], [0], marker='o', color='w', markersize=5, markerfacecolor='blue', label='Male')]
	plt.legend(handles=legend_elements)
	plt.savefig(figdir + 'total_cells_vs_total_counts_by_sex.png', dpi=300)

def pct_BT_cells_by_age(metadata_donors, figdir, var='T'):
	colors = {'male': 'blue', 'female': 'violet'}
	metadata_donors["color"] = metadata_donors["sex"].map(colors)

	plt.figure(figsize=(8, 6))
	plt.scatter(metadata_donors["age"], metadata_donors[f"percent_{var}_cells"],
				alpha=0.7, s=7, c=metadata_donors["color"], label=metadata_donors["sex"])

	# Separate the data by sex
	males = metadata_donors[metadata_donors["sex"] == "male"]
	females = metadata_donors[metadata_donors["sex"] == "female"]

	# Best fit line for males
	male_slope, male_intercept, male_r_value, _, _ = linregress(males["age"], males[f"percent_{var}_cells"])
	male_r_squared = male_r_value ** 2
	male_line_x = np.linspace(min(males["age"]), max(males["age"]), 100)
	male_line_y = male_slope * male_line_x + male_intercept
	plt.plot(male_line_x, male_line_y, color='blue', label=f'Male fit (R² = {male_r_squared:.2f})')

	# Best fit line for females
	female_slope, female_intercept, female_r_value, _, _ = linregress(females["age"], females[f"percent_{var}_cells"])
	female_r_squared = female_r_value ** 2
	female_line_x = np.linspace(min(females["age"]), max(females["age"]), 100)
	female_line_y = female_slope * female_line_x + female_intercept
	plt.plot(female_line_x, female_line_y, color='violet', label=f'Female fit (R² = {female_r_squared:.2f})')

	# Add titles and labels
	plt.title(f"Percent {var} cells by Age", fontsize=19)
	plt.xlabel("Age", fontsize=14)
	plt.ylabel(f"Percent {var} cells", fontsize=14)
	plt.grid(True, which="both", linestyle="--", linewidth=0.5)

	legend_elements = [
		Line2D([0], [0], marker='o', color='w', markersize=5, markerfacecolor=colors['female'], label='Female'),
		Line2D([0], [0], marker='o', color='w', markersize=5, markerfacecolor=colors['male'], label='Male')
	]
	plt.legend(handles=legend_elements, loc='upper left', fontsize=14)
	plt.savefig(figdir + f'pct_{var}_cells_by_age_by_sex.pdf')
	plt.close()

def avg_counts_per_cell(metadata, figdir):
	plt.clf()
	plt.figure()
	plt.hist(metadata['log_counts_per_cell'], bins=100, color="royalblue", edgecolor="black", alpha=0.75)
	plt.xlabel("Average Counts per Cell", fontsize=12)
	plt.ylabel("# Donors", fontsize=12)
	plt.title("Distribution of Counts per Cell", fontsize=14)
	plt.savefig(figdir + 'avg_counts_per_cell_hist.pdf')


def sex_pcts(metadata_donors, figdir):
	male_percentage = 100* (metadata_donors['sex']=='male').sum() / len(metadata_donors)
	female_percentage = 100 - male_percentage
	
	plt.figure(figsize=(6, 6))
	plt.pie(
		[male_percentage, female_percentage],
		labels=['Male', 'Female'],
		autopct='%d%%',
		colors=['skyblue', 'plum'], 
		startangle=90,
		textprops={'fontsize': 19}
	)
	#plt.title('Percentage of Male vs Female', fontsize=16)
	plt.tight_layout()
	plt.savefig(figdir+'percent_male_female.pdf')

def age_distribution(metadata_donors, figdir):
	plt.figure(figsize=(8, 6))
	plt.hist(metadata_donors['age'], bins=20, color='mediumslateblue', edgecolor='black', alpha=0.7)
	plt.title('Distribution of Age', fontsize=19)
	plt.xlabel('Age', fontsize=15)
	plt.ylabel('Number People', fontsize=15)
	plt.grid(axis='y', linestyle='--', linewidth=0.5)
	plt.savefig(figdir+'age_distribution.pdf')

def round_to_half(number):
	cap = round(number * 2) / 2
	return cap

#def plot_volcano_plot(ds, name, gene_names_dict, gene_chromosome_dict, figdir, padj_threshold=1e-3, log2fc_threshold=0.25, gate='or', info=None, adjust=False):
#	log2fc = ds.results_df["log2FoldChange"].values
#	padj = ds.results_df["padj"].values
#	ensembl_ids = ds.results_df.index
#	gene_names = [gene_names_dict[ensembl] for ensembl in ensembl_ids]
#	neg_log10_padj = -np.log10(padj)
#
#	if gate == 'or':
#		significant = (padj < padj_threshold) | (abs(log2fc) > log2fc_threshold)
#	else:
#		significant = (padj < padj_threshold) & (abs(log2fc) > log2fc_threshold)
#
#	# Assign colors based on chromosome location
#	colors = []
#	for gene, sig in zip(ensembl_ids, significant):
#		if not sig:
#			colors.append("gray")  # Non-significant genes in gray
#		else:
#			chrom = gene_chromosome_dict[gene]
#			if chrom == "X":
#				colors.append("red")
#			elif chrom == "Y":
#				colors.append("blue")
#			else:
#				colors.append("green")
#
#	# Cap log2fc values
#	cap_value = round_to_half(np.sort(np.abs(log2fc[~np.isnan(log2fc)]))[-3])
#	log2fc_clipped = np.clip(log2fc, -1 * cap_value, cap_value)
#	neg_log10_padj_clipped = np.clip(neg_log10_padj, 0, 50)
#
#	plt.figure(figsize=(8, 6))
#	plt.scatter(log2fc_clipped, neg_log10_padj_clipped, c=colors, alpha=0.7, s=6)
#
#	legend_patches = [
#		mpatches.Patch(color='blue', label='Y'),
#		mpatches.Patch(color='red', label='X'),
#		mpatches.Patch(color='green', label='Autosome')
#	]
#
#	plt.legend(handles=legend_patches)
#
#	# Draw threshold lines
#	plt.axvline(0, linestyle="--", color="black", linewidth=0.8)
#
#	# Label significant genes with improved text positioning
#	# Separate text annotations for high and low p-values
#	text_annotations_pos = []
#	text_annotations_neg = []
#	for i, gene in enumerate(gene_names):
#		if significant[i]:
#			if neg_log10_padj_clipped[i] >= 20:
#				# High p-values: Place label without adjusting
#				plt.text(log2fc_clipped[i], neg_log10_padj_clipped[i], gene, 
#						 fontsize=8, ha='right', va="bottom")
#			else:
#
#				if log2fc_clipped[i] > 0:
#					# Low p-values: Add to list for adjustment
#					text_pos = plt.text(log2fc_clipped[i], neg_log10_padj_clipped[i], gene, 
#									fontsize=8, ha='center', va="center")
#					text_annotations_pos.append(text_pos)
#				else:
#					# Low p-values: Add to list for adjustment
#					text_neg = plt.text(log2fc_clipped[i], neg_log10_padj_clipped[i], gene, 
#									fontsize=8, ha='center', va="center")
#					text_annotations_pos.append(text_neg)
#
#
#	plt.xlabel("Log2 Fold Change", fontsize=12)
#	plt.ylabel("-Log10 Adjusted p-value", fontsize=12)
#	plt.title(f"Volcano Plot of DESeq2 Results for {name} {info}", fontsize=14)
#
#	plt.grid(True, linestyle="--", linewidth=0.5)
#	plt.tight_layout()
#	print(text_annotations_pos)
#
#	# Adjust only the low p-value labels (neg_log10_padj < 20)
#	adjust_text(text_annotations_pos,
#				#avoid_self=True,
#				#force_explode=(2, 0),  
#				#expand=(2,0),
#				only_move='x',
#				#max_move=None,
#				arrowprops=dict(arrowstyle="-", color="gray", lw=0.3, alpha=0.5))
#	#adjust_text(text_annotations_neg,
#				#avoid_self=True,
#				#force_explode=(2, 0),  
#				#expand=(2,0),
#				#only_move='x-',
#				#max_move=None,
#				#arrowprops=dict(arrowstyle="-", color="gray", lw=0.3, alpha=0.5))
#
#	plt.savefig(figdir + f'deseqs_volcano_plot_{name}.pdf')

def plot_volcano_plot(results_df, name, gene_names_dict, gene_chromosome_dict, figdir, padj_threshold=1e-3, log2fc_threshold=0.25, gate='or', info=None):
	results_df = results_df.dropna()
	log2fc = results_df["log2FoldChange"].values
	padj = results_df["padj"].values
	ensembl_ids = results_df.index
	gene_names = [gene_names_dict[ensembl] for ensembl in ensembl_ids]
	neg_log10_padj = -np.log10(padj)

	if gate=='or':
		significant = (padj < padj_threshold) | (abs(log2fc) > log2fc_threshold)
	else:
		significant = (padj < padj_threshold) & (abs(log2fc) > log2fc_threshold)

	# Assign colors based on chromosome location
	colors = []
	for gene, sig in zip(ensembl_ids, significant):
		if not sig:
			colors.append("gray")  # Non-significant genes in gray
		else:
			chrom = gene_chromosome_dict[gene]
			if chrom == "X":
				colors.append("red") 
			elif chrom == "Y":
				colors.append("blue") 
			else:
				colors.append("green") 

	# Cap log2fc values
	cap_value = round_to_half(np.sort(np.abs(log2fc))[-3])
	if cap_value==0: cap_value = np.sort(np.abs(log2fc))[-1]
	log2fc_clipped = np.clip(log2fc, -1*cap_value, cap_value)
	neg_log10_padj_clipped = np.clip(neg_log10_padj, 0, 50)

	plt.figure(figsize=(8, 6))
	plt.scatter(log2fc_clipped, neg_log10_padj_clipped, c=colors, alpha=0.7, s=6)

	legend_patches = [
		mpatches.Patch(color='blue', label='Y'),
		mpatches.Patch(color='red', label='X'),
		mpatches.Patch(color='green', label='Autosome')
	]
	
	plt.legend(handles=legend_patches)

	# Draw threshold lines
	plt.axvline(0, linestyle="--", color="black", linewidth=0.8)

	# Label significant genes
	for i, gene in enumerate(gene_names):
		if significant[i]:
			plt.text(log2fc_clipped[i], neg_log10_padj_clipped[i], gene, fontsize=9, ha="right", va="bottom", rotation=-30)

	plt.xlabel("Log2 Fold Change", fontsize=12)
	plt.ylabel("-Log10 Adjusted p-value", fontsize=12)
	plt.title(f"Volcano Plot of DESeq2 Results for {name} {info}", fontsize=14)

	plt.grid(True, linestyle="--", linewidth=0.5)
	plt.tight_layout()

	plt.savefig(figdir + f'deseqs_volcano_plot_{name}.svg')


#def plot_volcano_plot(ds, name, gene_names_dict, chromosome_dict, figdir):
#	log2fc = ds.results_df["log2FoldChange"].values
#	padj = ds.results_df["padj"].values
#	ensembl_ids = ds.results_df.index
#	gene_names = [gene_names_dict[ensembl] for ensembl in ensembl_ids]
#	neg_log10_padj = -np.log10(padj)
#
#	padj_threshold = 1e-30
#	log2fc_threshold = 0.15  # Adjust as needed
#
#	significant = (padj < padj_threshold) | (abs(log2fc) > log2fc_threshold)
#	colors = np.where(significant, "red", "gray")
#
#	# Cap log2fc values at -1 and 1 for display
#	log2fc_clipped = np.clip(log2fc, -1, 1)
#	neg_log10_padj_clipped = np.clip(neg_log10_padj, 0, 50)
#
#	plt.figure(figsize=(8, 6))
#	plt.scatter(log2fc_clipped, neg_log10_padj_clipped, c=colors, alpha=0.7, s=5)
#
#	# Draw threshold lines
#	#plt.axhline(-np.log10(padj_threshold), linestyle="--", color="black", linewidth=0.8)  
#	plt.axvline(0, linestyle="--", color="black", linewidth=0.8)  
#	#plt.axvline(-log2fc_threshold, linestyle="--", color="black", linewidth=0.8)  
#	#plt.axvline(log2fc_threshold, linestyle="--", color="black", linewidth=0.8)
#
#	# Label significant genes
#	for i, gene in enumerate(gene_names):
#		if significant[i]:
#			plt.text(log2fc_clipped[i], neg_log10_padj_clipped[i], gene, fontsize=4, ha="right", va="bottom", rotation=-30)
#
#	plt.xlabel("Log2 Fold Change", fontsize=12)
#	plt.ylabel("-Log10 Adjusted p-value", fontsize=12)
#	plt.title(f"Volcano Plot of DESeq2 Results for {name}", fontsize=14)
#
#	plt.grid(True, linestyle="--", linewidth=0.5)
#	plt.tight_layout()
#
#	plt.savefig(figdir + f'deseqs_volcano_plot_{name}.pdf')


def plot_pval_distribution(pvalues, covar, figdir):
	plt.figure(figsize=(8, 5))
	sns.histplot(pvalues, bins=np.arange(0,1,0.02), kde=True, color="blue", edgecolor="black")
	plt.xlabel("P-value", fontsize=12)
	plt.ylabel("Frequency", fontsize=12)
	plt.title("Histogram of P-values, "+covar, fontsize=14)
	plt.savefig(figdir+'deseqs_pvalue_hist_'+covar+'.pdf')

def scale_dataframe(df):
	"""
	Scales each column of a Pandas DataFrame to the range [0, 1].

	Args:
		df (pd.DataFrame): The DataFrame to scale.

	Returns:
		 pd.DataFrame: A new DataFrame with scaled columns.
	"""
	scaler = MinMaxScaler()
	scaled_values = scaler.fit_transform(df)
	scaled_df = pd.DataFrame(scaled_values, columns=df.columns, index=df.index)
	return scaled_df

def plot_design_matrix(design_matrix, figdir):
	#design_matrix = design_matrix.drop(columns=['Intercept'])
	design_matrix = scale_dataframe(design_matrix)
	design_matrix['Intercept'] = 1
	sort_columns = ['Intercept', 'CD4_Naive', 'CD4_TCM', 'male', 'age_cat', 'percent_T_cells', 'percent_B_cells'] + \
				[col for col in design_matrix.columns if 'PC' in col]
	sorted_design_matrix = design_matrix.sort_values(by=sort_columns)

	plt.figure(figsize=(12, 6))
	sns.heatmap(sorted_design_matrix[sort_columns], cmap="viridis", xticklabels=True, yticklabels=False)
	plt.xticks(rotation=45, ha='right', fontsize=10)

	# Label the axes
	plt.xlabel("Features")
	plt.ylabel("Samples")
	plt.title("Design Matrix Heatmap")
	plt.tight_layout()
	plt.savefig(figdir+'design_matrix_heatmap.pdf')


def plot_design_matrix_kos(design_matrix, condition, figdir):
	design_matrix = scale_dataframe(design_matrix)
	design_matrix['Intercept'] = 1
	sort_columns = ['Intercept'] + \
				[col for col in design_matrix.columns if 'target' in col] + \
				['individualA', 'S_Score', 'G2M_Score']
	sorted_design_matrix = design_matrix.sort_values(by=sort_columns)

	plt.figure(figsize=(12, 6))
	sns.heatmap(sorted_design_matrix[sort_columns], cmap="viridis", xticklabels=True, yticklabels=False)
	plt.xticks(rotation=45, ha='right', fontsize=10)

	# Label the axes
	plt.xlabel("Features")
	plt.ylabel("Samples")
	plt.title("Design Matrix Heatmap for "+condition, fontsize=16)
	plt.tight_layout()
	plt.savefig(figdir+f'design_matrix_heatmap_{condition}.pdf')

