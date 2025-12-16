# Guide Knockdown Efficiency Supplementary Table

## Description
Summary statistics on knockdown efficiency of each sgRNA guide across three culture conditions.

- **index**: sgRNA ID
- **guide_mean_expr**: Mean log-normalized expression of the target gene in cells carrying this guide
- **guide_std_expr**: Standard deviation of log-normalized target gene expression in cells carrying this guide (set to 0.01 for guides with zero variance, 100 for guides with only one cell)
- **guide_n**: Number of cells carrying this guide
- **ntc_mean_expr**: Mean log-normalized expression of the target gene in non-targeting control cells
- **ntc_std_expr**: Standard deviation of log-normalized target gene expression in non-targeting control cells
- **ntc_n**: Total number of non-targeting control cells across all samples
- **t_statistic**: Welch's t-test statistic comparing guide expression vs NTC expression (negative values indicate knockdown)
- **p_value**: Nominal p-value from Welch's t-test
- **adj_p_value**: Benjamini-Hochberg FDR-adjusted p-value (minimum value capped at 1e-16)
- **signif_knockdown**: Boolean indicating significant knockdown (adj_p_value < 0.1 AND t_statistic < 0)
- **perturbed_gene_id**: Ensembl gene ID of the target gene
- **rank**: Rank of the target gene based on mean expression in NTC cells (1 = lowest expressed)
- **high_confidence_no_effect_guides**: Boolean indicating guides with high confidence of having no knockdown effect (criteria: non-significant knockdown, >10 cells with guide, target expression in NTCs >0.001)
- **culture_condition**: Culture condition for this measurement (Rest, Stim8hr, or Stim48hr)