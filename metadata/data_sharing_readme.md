# Data artifacts

## Cell-level data

Filenames: `D*_*.assigned_guide.h5ad`

Each AnnData object contains cell expression profiles for cells from one donor (D1, D2, D3, D4) and culture condition (Rest, Stim8hr, Stim48hr). Cells from different 10X lanes are concatenated. Each observation represents a cell. Each variable is a measured gene in the transcriptome.

### Observation Metadata (`.obs`)
Annotations for each single cell:

- **`lane_id`**: 10X lane identifier (corresponds to one cellranger output)
- **`n_genes_by_counts`**: Number of genes with non-zero counts detected in the cell
- **`total_counts`**: Total UMI counts in the cell
- **`pct_counts_mt`**: Percentage of counts mapping to mitochondrial genes
- **`top_guide_UMI_counts`**: UMI counts for the most abundant guide RNA in the cell
- **`guide_id`**: Unique identifier for the guide RNA detected in the cell (if more than one guide was detected, we annotate as "multi-guide")
- **`perturbed_gene_name`**: Name of the gene perturbed by the detected guide (before target curation)
- **`perturbed_gene_id`**: Ensembl gene ID of the perturbed gene (before target curation)
- **`guide_type`**: Type of guide (e.g., targeting, non-targeting)
- **`PuroR`**: Puromycin resistance marker expression level
- **`guide_group`**: Group classification for the guide 
- **`low_quality`**: Boolean flag indicating low-quality cells to be filtered

### Variable Metadata (`.var`)
Annotations for each measured gene:

- **`gene_ids`**: Ensembl gene identifiers
- **`feature_types`**: Type of feature (e.g., Gene Expression)
- **`genome`**: Reference genome used for alignment
- **`gene_name`**: Gene symbols
- **`mt`**: Boolean flag indicating mitochondrial genes

### Expression Matrix (`.X`)
Single-cell gene expression data:

- **Content**: UMI counts for each gene in each cell
- **Data type**: Sparse matrix (likely CSR format)

## Pseudobulk-level data

Filename: `GWCD4i.pseudobulk_merged.h5ad`

This AnnData object contains pseudobulk expression profiles. Each observation represents a pseudobulk (aggregated by guide, donor and culture condition). Each variable is a measured gene in the transcriptome (`n_vars = 18,129`).

### Observation Metadata (`.obs`)
Annotations for each pseudobulk sample:

- **`10xrun_id`**: processing batch identifier (R1 or R2)
- **`donor_id`**: Donor identifier
- **`culture_condition`**: Culture condition (Rest, Stim8hr, Stim48hr)
- **`guide_id`**: Unique guide identifier
- **`perturbed_gene_name`**: Name of the gene perturbed by the guide (note that the annotated gene in the guide identifier doesn't always match because we did some post-hoc curation of the target gene)
- **`perturbed_gene_id`**: Ensembl gene ID of the perturbed gene
- **`guide_type`**: Type of guide (e.g., targeting, non-targeting)
- **`n_cells`**: Number of cells aggregated in this pseudobulk sample
- **`total_counts`**: Total UMI counts across all cells in this pseudobulk
- **`log10_n_cells`**: Log10-transformed number of cells
- **`keep_min_cells`**: Boolean flag indicating sample passes minimum cell count threshold to be used for DE analysis
- **`keep_effective_guides`**: Boolean flag indicating guide was considered effective (t-test significant) to be used for DE analysis
- **`keep_total_counts`**: Boolean flag indicating sample passes total counts threshold to be used for DE analysis
- **`keep_for_DE`**: Boolean flag indicating sample is suitable for differential expression analysis
- **`keep_test_genes`**: Boolean flag indicating whether the perturbed gene passes criteria for differential expression analysis

### Variable Metadata (`.var`)
Annotations for each measured gene:

- **`gene_ids`**: Ensembl gene identifiers
- **`gene_name`**: Gene symbols

### Expression Matrix (`.X`)
Sum of UMI counts across cells for each gene in each pseudobulk sample

## Differential Expression Results 

Filename: `GWCD4i.DE_stats.h5ad`

This AnnData object contains genome-wide differential expression results from a perturb-seq experiment in CD4+ T cells. Each observation represents a single perturbation (perturbed gene) tested in a specific culture condition (`n_obs = 33,983`). Each variable is a measured gene in the transcriptome (`n_vars = 10,282`).

### Observation Metadata (`.obs`)
Annotations for each perturbation-condition pair:

- **`target_contrast_gene_name`**: Name of the perturbed gene  
- **`culture_condition`**: culture condition (Rest, Stim8hr, Stim48hr)  
- **`target_contrast`**: Unique identifier for the perturbed gene  
- **`chunk`**: differential expression processing group identifier  
- **`n_cells_target`**: Number of cells with targeting guide for the perturbed gene  
- **`n_up_genes`**: Count of significantly upregulated genes (10% FDR)  
- **`n_down_genes`**: Count of significantly downregulated genes (10% FDR)  
- **`n_total_de_genes`**: Total number of significantly differentially expressed genes (10% FDR)  
- **`ontarget_effect_size`**: Effect size of the perturbation on its intended target gene  
- **`ontarget_significant`**: Boolean indicating whether on-target knockdown was significant (10% FDR)  
- **`target_baseMean`**: Mean baseline expression of the target gene  
- **`offtarget_flag`**: Flag indicating potential off-target effects (TSS within 10 kb with significant down-regulation)  
- **`n_total_genes_category`**: Category based on number of trans-effects  
- **`ontarget_effect_category`**: Category based on on-target / off-target effects  
- **`n_downstream`**: Number of genes significantly affected by this perturbation, excluding on-target efffect (incoming trans-effects)

### Variable Metadata (`.var`)
Annotations for each measured gene:

- **`gene_ids`**: Gene identifiers (e.g., Ensembl IDs)
- **`gene_name`**: Gene symbols

### Variable Matrices (`.varm`)
Summary statistics for measured genes across conditions:

- **`measured_genes_stats_Stim8hr`**: Gene-level statistics for 8-hour stimulation condition
- **`measured_genes_stats_Stim48hr`**: Gene-level statistics for 48-hour stimulation condition
- **`measured_genes_stats_Rest`**: Gene-level statistics for resting/unstimulated condition

### Data Layers (`.layers`)
Differential expression statistics for each perturbation-gene pair (from DESeq2):

- **`log_fc`**: Log2 fold change
- **`p_value`**: Raw p-values from differential expression testing
- **`adj_p_value`**: FDR-adjusted p-values
- **`baseMean`**: Mean normalized expression of the gene across cells
- **`lfcSE`**: Standard error of log fold change
- **`zscore`**: Z-scores for differential expression (logFC / lfcSE)


## Supplementary tables

### Sample metadata

Filename: `sample_metadata.suppl_table.csv`

This supplementary table contains experimental metadata for all samples in the perturb-seq screen. Each row represents a unique biological sample with information about the experimental setup, library preparation, sequencing details, and donor demographics.

- **`cell_sample_id`**: Unique identifier for the biological sample
- **`10xrun_id`**: Unique identifier for run/batch (R1 or R2)
- **`donor_id`**: Donor identifier
- **`culture_condition`**: Culture condition applied to the cells (Rest, Stim8hr, Stim48hr)
- **`library_id`**: Unique identifier for the sequencing library (matches cellranger outputs)
- **`library_prep_kit`**: Library preparation kit used for sample processing (e.g., GEMX_flex_v2)
- **`probe_hyb_loading`**: Probe hybridization loading information (cell count and probe details)
- **`GEM_loading`**: GEM loading information for 10x Genomics workflow
- **`sequencing_platform`**: Sequencing platform used (e.g., Ultima)
- **`age`**: Donor age in years
- **`sex`**: Donor sex (Male/Female)
- **`ethnicity`**: Donor ethnicity
- **`weight_kg`**: Donor weight in kilograms
- **`height_cm`**: Donor height in centimeters
- **`smoker`**: Smoking status (Yes/No)
- **`blood type`**: Donor blood type
- **`anticoagulant`**: Anticoagulant used for blood collection
- **`harvest_date`**: Date of blood sample collection

### Differential expression statistics for each perturbation-condition pair

Filename: `DE_stats.suppl_table.csv`

See `.obs` of "Differential expression results"

### Guide library metadata

Filename: `sgrna_library_metadata.suppl_table.csv`

Contains metadata for the sgRNA guide library used in the genome-wide CRISPR perturbation screen. Each row represents a single guide RNA with its genomic targeting information, design details, and potential off-target considerations.

- **`sgRNA`**: Unique identifier for the guide RNA
- **`chromosome`**: Chromosome of the target site
- **`pos`**: Genomic position of the guide target site
- **`strand`**: DNA strand orientation of the target site (+ or -)
- **`seq`**: Full guide RNA sequence
- **`seq_last19bp`**: Last 19 base pairs of the guide sequence
- **`PAM`**: boolean flag for presence of Protospacer Adjacent Motif sequence
- **`note`**: Additional notes about the guide design
- **`flag`**: Quality control or classification flag
- **`target_gene_name_from_sgRNA`**: Target gene name derived from the sgRNA identifier
- **`designed_target_gene_id`**: Ensembl gene ID of the intended target gene (as designed)
- **`designed_target_gene_name`**: Gene name of the intended target gene (as designed)
- **`target_gene_id`**: Ensembl gene ID of the actual/validated target gene
- **`target_gene_name`**: Gene name of the actual/validated target gene
- **`distance_to_closest_target_tss`**: Distance (in base pairs) from guide to the closest transcription start site (TSS) of the target gene
- **`nearby_gene_within_2kb`**: Boolean or count indicating genes within 2 kb of the guide target site
- **`nearby_gene_within_30kb`**: Boolean or count indicating genes within 30 kb of the guide target site
- **`nearest_within2kb_gene_id`**: Ensembl gene ID of the nearest gene within 2 kb
- **`nearest_within2kb_gene_name`**: Gene name of the nearest gene within 2 kb
- **`nearest_within2kb_gene_dist`**: Distance to the nearest gene within 2 kb
- **`nearest_within2kb_nontarget_gene_id`**: Ensembl gene ID of the nearest non-target gene within 2 kb
- **`nearest_within2kb_nontarget_gene_name`**: Gene name of the nearest non-target gene within 2 kb
- **`nearest_within2kb_nontarget_gene_dist`**: Distance to the nearest non-target gene within 2 kb
- **`putative_bidirectional_promoter`**: Flag indicating potential bidirectional promoter region (may affect multiple genes)
- **`other_alignment_chromosome`**: Chromosome with potential off-target alignment
- **`other_alignment_pos`**: Genomic position of potential off-target alignment

### Guide knockdown efficiency

Filename: `guide_kd_efficiency.suppl_table.csv`

Summary statistics on knockdown efficiency of each sgRNA guide across three culture conditions.

- **`index`**: sgRNA ID
- **`guide_mean_expr`**: Mean log-normalized expression of the target gene in cells carrying this guide
- **`guide_std_expr`**: Standard deviation of log-normalized target gene expression in cells carrying this guide (set to 0.01 for guides with zero variance, 100 for guides with only one cell)
- **`guide_n`**: Number of cells carrying this guide
- **`ntc_mean_expr`**: Mean log-normalized expression of the target gene in non-targeting control cells
- **`ntc_std_expr`**: Standard deviation of log-normalized target gene expression in non-targeting control cells
- **`ntc_n`**: Total number of non-targeting control cells across all samples
- **`t_statistic`**: Welch's t-test statistic comparing guide expression vs NTC expression (negative values indicate knockdown)
- **`p_value`**: Nominal p-value from Welch's t-test
- **`adj_p_value`**: Benjamini-Hochberg FDR-adjusted p-value (minimum value capped at 1e-16)
- **`signif_knockdown`**: Boolean indicating significant knockdown (adj_p_value < 0.1 AND t_statistic < 0)
- **`perturbed_gene_id`**: Ensembl gene ID of the target gene
- **`rank`**: Rank of the target gene based on mean expression in NTC cells (1 = lowest expressed)
- **`high_confidence_no_effect_guides`**: Boolean indicating guides with high confidence of having no knockdown effect (criteria: non-significant knockdown, >10 cells with guide, target expression in NTCs >0.001)
- **`culture_condition`**: Culture condition for this measurement (Rest, Stim8hr, or Stim48hr)

### CD4+ T cell aging signature differential expression results

Filename: `CD4T_aging_signature_DE_results_full.suppl_table.csv`

Full differential expression results for DE analysis of age-associated changes in CD4+ T cells across all cohorts.

- **`variable`**: Ensembl gene ID of the measured gene
- **`gene_name`**: Gene symbol
- **`baseMean`**: Mean baseline expression of the gene
- **`log_fc`**: Log2 fold change
- **`lfcSE`**: Standard error of log fold change
- **`stat`**: Test statistic
- **`p_value`**: Raw p-value from differential expression testing
- **`adj_p_value`**: FDR-adjusted p-value
- **`contrast`**: comparison cohort
- **`zscore`**: Z-score for differential expression (log_fc / lfcSE)

### Th2/Th1 polarization signature differential expression results

Filename: `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv`


Full differential expression results for DE analysis of Th2 vs Th1 changes in CD4+ T cells across all cohorts.


- **`variable`**: Gene symbol
- **`baseMean`**: Mean baseline expression of the gene
- **`log_fc`**: Log2 fold change
- **`lfcSE`**: Standard error of log fold change
- **`stat`**: Test statistic
- **`p_value`**: Raw p-value from differential expression testing
- **`adj_p_value`**: FDR-adjusted p-value
- **`contrast`**: comparison cohort
- **`zscore`**: Z-score for differential expression (log_fc / lfcSE)

### Cluster autoimmune disease enrichment results

Filename: `cluster_autoimmune_enrichment_results.suppl_table.csv`

Enrichment analysis results for autoimmune disease-associated genes within perturbation effect clusters.

- **`cluster`**: Cluster identifier
- **`disease`**: Disease category (autoimmune disease)
- **`gene_set`**: Gene set being tested (downstream effects by condition)
- **`odds_ratio`**: Odds ratio from Fisher's exact test
- **`ci_low`**: Lower bound of 95% confidence interval for odds ratio
- **`ci_high`**: Upper bound of 95% confidence interval for odds ratio
- **`p_value`**: Raw p-value from Fisher's exact test
- **`p_adj_fdr`**: FDR-adjusted p-value
- **`cluster_size`**: Number of genes in the cluster
- **`in_cluster_in_disease`**: Count of genes both in cluster and associated with disease
- **`in_cluster_not_disease`**: Count of genes in cluster but not associated with disease
- **`not_cluster_in_disease`**: Count of disease-associated genes not in cluster
- **`not_cluster_not_disease`**: Count of genes neither in cluster nor associated with disease
- **`intersecting_genes`**: List of genes that overlap between cluster and disease association
- **`negative_control_disease`**: Boolean flag indicating if this is a negative control disease category

### Aging prediction regulator coefficients

Filename: `aging_prediction_condition_comparison_regulator_coefficients.csv`

Model coefficients from linear models predicting the CD4+ T cell aging signature across different datasets (perturb-seq in CD4+ T cells vs K562 cells).

- **`coef_mean`**: Mean coefficient value for the regulator across model fits
- **`coef_sem`**: Standard error of the mean for the coefficient
- **`coef_rank`**: Rank of the regulator coefficient (0-1 scale, higher = stronger effect)
- **`regulator`**: Gene symbol of the regulator
- **`known_regulators`**: Boolean indicating if this is a known regulator of aging
- **`dataset_key`**: Dataset identifier for model comparison (e.g., CD4T_K562)
- **`regulator_type`**: Type/category of regulator
- **`celltype`**: Cell type or condition context (K562, Rest, Stim8hr, Stim48hr)
- **`signature`**: Signature being predicted (CD4T)

### Polarization prediction regulator coefficients

Filename: `polarization_prediction_condition_comparison_regulator_coefficients.csv`

Model coefficients from linear models predicting T cell activation and polarization signatures across different culture conditions.

- **`coef_mean`**: Mean coefficient value for the regulator across model fits
- **`coef_sem`**: Standard error of the mean for the coefficient
- **`coef_rank`**: Rank of the regulator coefficient (0-1 scale, higher = stronger effect)
- **`regulator`**: Gene symbol of the regulator
- **`known_regulators`**: Boolean indicating if this is a known regulator of the signature
- **`dataset_key`**: Dataset identifier for model comparison (e.g., activation_Rest, polarization_Stim8hr)
- **`regulator_type`**: Type/category of regulator
- **`celltype`**: Culture condition context (Rest, Stim8hr, Stim48hr)
- **`signature`**: Signature being predicted (activation or polarization)
