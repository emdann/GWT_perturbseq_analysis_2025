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

- On Dropbox: `GRNPerturbSeq/3_expts/processed_data/{experiment_name}` 
- On oak: `/oak/stanford/groups/pritch/users/emma/data/GWT/{experiment_name}`

### Count matrices / AnnData objects
- `{sample}.{lane}scRNA.postQC_obs.h5ad` - count matrices and annotations after QC with sgRNA assignment

### QC stats
- `QC_summary_stats.csv` - summary of QC metric statistics for each sample and lane
- `perturbation_counts.csv` - count of number of cells per perturbation for each sample and lane
- `{EXPERIMENT_NAME}.guide_effect.{culture_condition}.csv` - summary stats to assess sgRNA effect on target gene compared to expression of gene in NTC controls
- `no_effect_guides.txt` - guides with no significant effect in any condition

### Differential expression analysis files
- `{experiment_name}_merged.DE_pseudobulk.h5ad` - pseudobulked gene expression counts per guide+sample+condition (summing expression profile)
- `DE_results_{run_name}/{experiment_name}.gex.lognorm.h5ad` - DE analysis results (obs are perturbations x condition, vars are transcriptome genes)
- `DE_results_{run_name}/DE_summary_stats_per_target.csv` - Summary of on-target effects and overall effect for each perturbation and condition


To sync processed data with Dropbox

```bash
# DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
DATADIR=/mnt/oak/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR1_Psomagen
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/
DROPBOX_PATH=GRNPerturbSeq/3_expts/processed_data/

# Define list of files to copy
FILES_TO_COPY=(
    "${EXPERIMENT_NAME}_merged.gex.lognorm.postQC_obs.csv"
    "${EXPERIMENT_NAME}_merged.DE_pseudobulk.h5ad"
    "DE_results/${EXPERIMENT_NAME}.merged_DE_results.h5ad"
    "knockdown_efficacy_simple.csv"
    "guide_ontarget_effect_simple.csv"
    "${EXPERIMENT_NAME}_merged.gex.h5ad"
)

rclone mkdir dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/

# Copy each file to Dropbox
for f in "${FILES_TO_COPY[@]}"; do 
    rclone copy ${EXPDIR}/${f} dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/ --checksum --ignore-times
done

# Sync raw outputs
DATADIR=/mnt/oak/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR1_Psomagen
DROPBOX_PATH=GRNPerturbSeq/3_expts/
LOCAL_DATA_DIR="$DATADIR/czi-psomagen"
rclone mkdir dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/

find "$LOCAL_DATA_DIR" -name "web_summary.html" -type f | while read file; do
    relative_path=${file#$LOCAL_DATA_DIR/}
    relative_dir=$(dirname "$relative_path")
    rclone copy "$file" "dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/${relative_dir}/" --checksum --ignore-times
done
```





