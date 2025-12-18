
module load system cairo/1.14.10
module load python/3.12.1

##python3
import warnings

warnings.filterwarnings("ignore")
import anndata as ad

adata = ad.read_h5ad("CD4i_final.merged_DE_results.h5ad")

import pandas as pd

logfc=adata.layers["log_fc"]
df = pd.DataFrame(logfc)
df.to_csv("CD4i_final_merged.DE_pseudobulk_logFC.csv")

Pval=adata.layers["p_value"]
df = pd.DataFrame(Pval)
df.to_csv("CD4i_final_merged.DE_pseudobulk_pvalue.csv")

lfcse=adata.layers["lfcSE"]
df = pd.DataFrame(lfcse)
df.to_csv("CD4i_final_merged.DE_pseudobulk_lfcSE.csv")

adata.obs.to_csv("CD4i_final_merged.DE_pseudobulk_metadata.csv")
adata.var.to_csv("CD4i_final_merged.DE_pseudobulk_genes.csv")


##put Backman_LymphocyteCount_fullFeatures.per_gene_estimates.tsv in input/burden

mkdir -p BurdenRegCor/GeneLevel/

file_list=$( ls input/burden )

for FILE in $file_list; do
sbatch  --error="Regulator_burden_cor.err"  --output="Regulator_burden_cor.log" --time=24:00:00  --mem=64GB  R_4.2_execute.sh  Regulator_burden_correlation_GWT.R ${FILE} 
done

file_list2=$( ls | grep limma )
for FILE in $file_list; do
	for FILE2 in $file_list2; do
sbatch  --error="Regulator_burden_cor.err"  --output="Regulator_burden_cor.log" --time=24:00:00  --mem=64GB  R_4.2_execute.sh  Regulator_burden_correlation_CelLine.R ${FILE} ${FILE2} 
done
done