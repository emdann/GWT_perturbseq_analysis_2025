# GWT_perturbseq_analysis
Analysis of genome-wide perturb-seq screen on primary T cells (see our [manuscript]())

## Contents

- `src` - analysis code
    - `1_preprocess/` - ingest and preprocess new experiments from cellranger outputs
    - `2_embedding/` - cell state embedding
    - `3_DE_analysis/` - differential expression analysis
    - `4_polarization_signatures/` - analysis of polarization signatures
    - `5_cytokine_regulators/` - analysis of cytokine regulators
    - `6_functional_interaction/` - functional interaction analysis
    - `7_1k1k_analysis/` - 1k1k dataset analysis
    - `8_lymphocyte_counts_LoF/` - lymphocyte counts loss-of-function analysis
    - `_misc/` - miscellaneous utility scripts
- `metadata` - sample and experimental metadata, configs, gene annotations etc

Please refer to the [figure map](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/figure_map.md) to find which scripts were used to generate a specific figure in the manuscript.

## Set-up compute environment

```
conda env create -f environment.yaml
conda activate gwt-env
```

## Data pointers

- Processed data (cell-level count matrices, pseudobulk-level count matrices and differential expression estimates, analysis results) are available via the [Biohub Virtual Cells Platform](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq). A detailed description of shared files can be found [here](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/data_sharing_readme.md).
- Additional supplementary tables and metadata are available [here](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata)
- Raw sequencing data and cellranger outputs will be made available through SRA/GEO (accession: SRP643211 / GSE314342)


## Citation

If you use this data or code in your work, please cite

Zhu R., Dann E. et al. (2025) Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits. _bioRxiv_

## Contact

For any questions, please post an [issue](https://github.com/emdann/GWT_perturbseq_analysis_2025/issues?q=sort%3Aupdated-desc+is%3Aissue+is%3Aopen) in this repository, or contact by email `emmadann<at>stanford.edu` or `ronghui.zhu<at>gladstone.ucsf.edu`. 