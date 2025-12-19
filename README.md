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
    - `_misc/` - miscellaneous / old code
- `metadata` - sample and experimental metadata, configs, gene annotations etc

## Set-up compute environment

```
conda env create -f environment.yaml
conda activate gwt-env
```

## Citation

If you use this data or code in your work, please cite

Zhu R., Dann E. et al. (2025) Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits. _bioRxiv_

## Contact

For any questions, please post an [issue](https://github.com/emdann/GWT_perturbseq_analysis_2025/issues?q=sort%3Aupdated-desc+is%3Aissue+is%3Aopen) in this repository, or contact by email `emmadann<at>stanford.edu` or `ronghui.zhu<at>gladstone.ucsf.edu`. 