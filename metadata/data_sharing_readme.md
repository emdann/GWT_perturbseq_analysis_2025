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

This supplementary table contains experimental metadata for all samples in the perturb-seq screen. Each row represents a unique biological sample with information about the experimental setup, library preparation, and sequencing details.

- **`10xrun_id`**: Unique identifier for run/batch (R1 or R2)
- **`cell_sample_id`**: Unique identifier for the biological sample 
- **`donor_id`**: Donor identifier for the donor
- **`culture_condition`**: Culture condition applied to the cells (Rest, Stim8hr, Stim48hr)
- **`library_id`**: Unique identifier for the sequencing library (matches cellranger outputs)
- **`library_prep_kit`**: Library preparation kit used for sample processing
- **`probe_hyb_loading`**: Probe hybridization loading information 
- **`GEM_loading`**: GEM loading information for 10x Genomics workflow
- **`sequencing_platform`**: Sequencing platform used

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