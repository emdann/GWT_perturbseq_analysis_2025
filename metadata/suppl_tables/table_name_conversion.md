# Conversion lookup table

Mapping from source file (this directory) to the renamed file in
`../../../../5_manuscripts/supplementary_table_renamed/`

## Source -> renamed

sample_metadata.suppl_table.csv -> Supplementary_Table_1_Sample_metadata.csv
stabl_constructs.csv -> Supplementary_Table_2_Plasmid_constructs.csv
sgrna_library_metadata.suppl_table.csv -> Supplementary_Table_3_gRNA_guide_library_metadata.csv
QC_summaries_per_sample_lane.csv -> Supplementary_Table_4_QC_summaries_per_sample_lane.csv
guide_kd_efficiency.suppl_table.csv -> Supplementary_Table_5_guide_kd_efficiency.csv
guide_offtarget_analysis_results.csv -> Supplementary_Table_6_guide_off_target.csv
DE_stats.suppl_table.csv -> Supplementary_Table_7_DE_stats.csv
K562_comparison.suppl_table.csv -> Supplementary_Table_8_K562_comparison.csv
IL10IL21bulkRNAseq_DESeq2_results_byguide.csv -> Supplementary_Table_9_DE_byguide_IL10IL21_validation.csv
IL10IL21bulkRNAseq_DESeq2_results_bygene.csv -> Supplementary_Table_10_DE_bygene_IL10IL21_validation.csv
IL10_IL21_arrayed_validation_combined.csv -> Supplementary_Table_11_IL10IL21_flow_cytometry_validation.csv
proliferation_stats_IL10IL21.csv -> Supplementary_Table_12_Proliferation_IL10IL21_validation.csv
clustering_results_and_annotations.csv -> Supplementary_Table_13_clustering_results_and_annotations.csv
clustering_downstream_genes.csv.gz -> Supplementary_Table_14_clustering_downstream_genes.csv.gz
clusterTCR_deseq2_results.csv -> Supplementary_Table_15_DE_TCRcluster_validation.csv
Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv -> Supplementary_Table_16_Th2_Th1_polarization_signature_DE_results_full.csv
polarization_prediction_condition_comparison_regulator_coefficients.csv -> Supplementary_Table_17_polarization_prediction_condition_comparison_regulator_coefficients.csv
Th1Th2_validation_summary.suppl_table.csv -> Supplementary_Table_18_Th1Th2_arrayed_validation_summary.csv
Th1Th2bulkRNAseq_DESeq2_results.csv.gz -> Supplementary_Table_19_DE_Th1Th2_validation_bulkRNAseq.csv.gz
CD4T_aging_signature_DE_results_full.suppl_table.csv -> Supplementary_Table_20_CD4T_aging_signature_DE_results_full.csv
aging_prediction_condition_comparison_regulator_coefficients.csv -> Supplementary_Table_21_aging_prediction_condition_comparison_regulator_coefficients.csv
cluster_autoimmune_enrichment_results.suppl_table.csv -> Supplementary_Table_22_cluster_autoimmune_enrichment_results.csv
validation_donor_metadata.csv -> Supplementary_Table_23_followup_donor_metadata.csv

# Instruction
Output to ../../../../5_manuscripts/supplementary_table_renamed with file name renamed according to conversion lookup table. Any compressed file (e.g. .csv.gz) should be decompressed in this output directory.
Output to ../../../GWT_perturbseq_analysis_2025/metadata/suppl_table with the original name, for files larger than 50Mb, compress
