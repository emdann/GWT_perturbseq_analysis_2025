# GWT_perturbseq_analysis
Analysis of genome-wide perturb-seq screen on primary T cells (see our [manuscript](https://urldefense.com/v3/__https://kwnsfk27.r.eu-west-1.awstrack.me/L0/https:*2F*2Fauthors.elsevier.com*2Fsd*2Farticle*2FS0092-8674(26)00929-3/1/010201a048e35ea6-6688df17-b320-4de3-8aba-096d4f5be26f-000000/rURbivJIptX9Fm8r7lDPxZmCeow=473__;JSUlJSU!!G92We9drHetJ8EofZw!fAtxxwDXcuX14X51Gp86Giy2pLU--OiX92LO_zkFiL7vOpdBUxqMO0ZtwUv28RxazWc_YK3CpWt0N4iW_erRsas$))

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
    - `9_ChIP_analysis/` - ChIP-seq analysis of direct/indirect trans-effects
    - `_misc/` - miscellaneous utility scripts
- `metadata` - sample and experimental metadata, configs, gene annotations etc

Please refer to the [figure map](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/figure_map.md) to find which scripts were used to generate a specific figure in the manuscript.

## Set-up compute environment

```
conda env create -f environment.yaml
conda activate gwt-env
```

## Data pointers

Processed data (cell-level count matrices, pseudobulk-level count matrices and differential expression estimates, analysis results) are available via the [Biohub Virtual Cells Platform](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq). Run the following AWS CLI command in your terminal to explore the available data:
```
aws s3 ls --no-sign-request s3://genome-scale-tcell-perturb-seq/marson2025_data/
``` 
A detailed description of shared files can be found [here](https://github.com/emdann/GWT_perturbseq_analysis_2025/blob/master/metadata/data_sharing_readme.md).

Additional supplementary tables and metadata are available [here](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata)
Raw sequencing data and cellranger outputs will be made available through SRA/GEO (accession: SRP643211 / GSE314342)

## Citation

If you use this data or code in your work, please cite

Zhu R., Dann E. et al. (2026) Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits. _Cell_

## Contact

For any questions, please post an [issue](https://github.com/emdann/GWT_perturbseq_analysis_2025/issues?q=sort%3Aupdated-desc+is%3Aissue+is%3Aopen) in this repository, or contact by email `emmadann<at>stanford.edu` or `ronghui.zhu<at>gladstone.ucsf.edu`. 
