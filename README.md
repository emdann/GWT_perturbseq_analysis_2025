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

## Contributing code

Copy directory in your local environment 
```
git clone git@github.com:emdann/GWT_perturbseq_analysis.git
```

To add new code/notebooks, add them to `src` in an appropriate directory, then run:
```
git add src/path/to/new_script.py
git commit -m 'adding script for analysis x'
git push origin master
```

To edit existing code/notebooks without creating conflicts, make edits, then run:
```
git checkout -b new-branch-name
git add .
git commit -m 'esiting script x y and z'
git push origin new-branch-name
```

## Processed data

### Count matrices / AnnData objects

These objects are saved and processed by 10X lane. For now they are stored under two separate experiment names, one for each data drop (`CD4iR1_Psomagen` and `CD4iR2_Psomagen`)

Base path:
- Dropbox `GRNPerturbSeq/3_expts/processed_data/{experiment_name}` 
- oak `/oak/stanford/groups/pritch/users/emma/data/GWT/{experiment_name}`

Files:
- `tmp/{sample}.{lane}scRNA.postQC.h5ad` - count matrices and annotations after QC with sgRNA assignment
- `QC_summary_stats.csv` - summary of QC metric statistics for each sample and lane
- `perturbation_counts.csv` - count of number of cells per perturbation for each sample and lane

### Guide effect estimates 

Base path:
- Dropbox: `GRNPerturbSeq/3_expts/processed_data/CD4i_final/` 
- oak: `/oak/stanford/groups/pritch/users/emma/data/GWT/CD4i_final/`

Files:
- `CD4i_final.guide_effect.{culture_condition}.csv` - summary stats to assess sgRNA effect on target gene compared to expression of gene in NTC controls
- `no_effect_guides.txt` - guides with no significant effect in any condition

### Differential expression analysis files

Base path:
- Dropbox: `GRNPerturbSeq/3_expts/processed_data/CD4i_final/` 
- oak: `/oak/stanford/groups/pritch/users/emma/data/GWT/CD4i_final/`

Files:
- `CD4i_final_merged.DE_pseudobulk.h5ad` - pseudobulked gene expression counts per guide+sample+condition (summing expression profile)
- `DE_results_all_confounders/CD4i_final_merged.DE_results.h5ad` - DE analysis results (obs are perturbations x condition, vars are transcriptome genes)
- `DE_results_all_confounders/DE_summary_stats_per_target.csv` - Summary of on-target effects and overall effect for each perturbation and condition


To sync processed data with Dropbox, see [example script](https://github.com/emdann/GWT_perturbseq_analysis/blob/master/src/_misc/sync2dropbox.sh).  
