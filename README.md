# GWT_perturbseq_analysis
Analysis of genome-wide perturb-seq screen on primary T cells

## Contents

- `src` - analysis code
    - `1_preprocess/` - ingest and preprocess new experiments from cellranger outputs
    - `2_embedding` - cell state embedding
    - `3_DE_analysis` - differential expression analysis
    - `4_polarization_analysis` - analysis of polarization signatures
    - `_misc` - miscellaneous / old code
- `metadata` - sample and experimental metadata, configs, gene annotations etc

## Set-up compute environment

To install required packages (including [perturbseq_tools module](https://github.com/emdann/perturbseq_tools))

```
conda env create -f environment.yaml
conda activate gwt-env
```