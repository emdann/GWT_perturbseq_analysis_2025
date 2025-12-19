# Data ingestion and preprocessing

## Set up 

1. Update experiment config file: add new experiment in `metadata/experiments_config.yaml` - new entry should be called as EXPERIMENT_NAME
2. Make folder structure
```bash
EXPERIMENT_NAME=CD4iR2_Psomagen
python make_GWT_directories.py $EXPERIMENT_NAME
```
3. Process sample metadata 
```bash
# first download GWT_sample_metadata.xlsx 
python process_sample_metadata.py --experiment_name $EXPERIMENT_NAME --datadir /mnt/oak/users/emma/data/GWT/
```
4. Make sgRNA library metadata
```bash
sgRNA_annotation/prep_sgrna_library_metadata.ipynb
```

## Run preprocessing workflow locally

No SLURM (suitable only for small experiments)

```bash
# Ingest cellranger outputs and merge
python preprocess.py --config ../../metadata/experiments_config.yaml --experiment $EXPERIMENT_NAME --merge --embedding

# Guide assignments
python sgrna_assignment.py $EXPERIMENT_NAME --config ../../metadata/experiments_config.yaml

# QC analysis
# see notebook qc_final.ipynb 
```

## Run preprocessing workflow with SLURM

This runs preprocessing and QC on each .h5 cellranger output in parallel.

1. Ingest and basic preprocessing of cellranger outputs
```bash
conda activate gwt-env
EXPERIMENT_NAME=CD4iR2_Psomagen
DATADIR=/path/to/data/GWT/
for H5_FILE in $(ls ${DATADIR}/${EXPERIMENT_NAME}/cellranger_outs/*/*); do
  sbatch \
      --partition=pritch \
      --job-name=preprocess_${EXPERIMENT_NAME} \
      --output=$GROUP_SCRATCH/emma/slurm-process_%j.out \
      --error=$GROUP_SCRATCH/emma/slurm-process_%j.err \
      --mem=100G  \
      --time=01:00:00 \
      --wrap="python preprocess.py --experiment ${EXPERIMENT_NAME} --input_h5 ${H5_FILE} --force"
done
```

2. Guide RNA assignment in chunks
```bash
EXPDIR=/path/to/data/GWT/${EXPERIMENT_NAME}/
for h5ad_file in $(ls $EXPDIR*.sgRNA.h5ad); do
  ./submit_sgrna_assignment.sh $EXPERIMENT_NAME $h5ad_file
done
```
3. Merge to cell-level guide assignment for each sample 
```bash
python sgrna_assignment.py $EXPERIMENT_NAME --merge
```

4. Compute QC stats and exclude low quality cells
```bash
H5AD_FILES=$(ls ${DATADIR}/${EXPERIMENT_NAME}/tmp/*.scRNA.h5ad)
for f in $H5AD_FILES; do
  SAMPLE_NAME=$(basename ${f} .scRNA.h5ad)
  INPUT_CSV="${DATADIR}/${EXPERIMENT_NAME}/${SAMPLE_NAME}.sgrna_assignment.csv"
  OUTPUT_H5AD="${DATADIR}/${EXPERIMENT_NAME}/tmp/${SAMPLE_NAME}.scRNA.postQC.h5ad"
  
  if [ ! -f "${OUTPUT_H5AD}" ] && [ -f "${INPUT_CSV}" ]; then
    echo $SAMPLE_NAME
    sbatch \
          --partition=pritch \
          --job-name=qc_${EXPERIMENT_NAME}_${SAMPLE_NAME} \
          --output=$GROUP_SCRATCH/emma/slurm-qc_%j.out \
          --error=$GROUP_SCRATCH/emma/slurm-qc_%j.err \
          --mem=24G  \
          --time=01:00:00 \
          --wrap="python qc_samples.py --experiment_name=${EXPERIMENT_NAME} --sample_id=${SAMPLE_NAME}"
  fi
done
```

## QC analysis

- `qc_final.ipynb` - Plots on quality control stats for each experiment
- `estimate_guide_effect.ipynb` - Estimate KD effect of each guide (to exclude ineffective from DE analysis)

## Outputs
- `QC_summary_stats.csv` - summary of QC metric statistics for each sample and lane
- `perturbation_counts.csv` - count of number of cells per perturbation for each sample and lane
- `{EXPERIMENT_NAME}.guide_effect.{culture_condition}.csv` - summary stats to assess sgRNA effect on target gene compared to expression of gene in NTC controls
- `no_effect_guides.txt` - guides with no significant effect in any condition