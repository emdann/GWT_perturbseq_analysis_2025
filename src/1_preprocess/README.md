## Moving cellranger outs from Dropbox to Sherlock

1. Update experiment config file: add new experiment in `metadata/experiments_config.yaml` - new entry should be called as EXPERIMENT_NAME

```bash
# Setup folders in Sherlock
EXPERIMENT_NAME=CD4iR1_Psomagen
python make_GWT_directories.py $EXPERIMENT_NAME

# DROPBOX_PATH=GRNPerturbSeq/3_expts/CRiCD4IL2_Illumina/lane13_10x/
DROPBOX_PATH=correct_spikein_sequence/lane5_10x_puro_sequence_corrected
rclone copy dropbox:"${DROPBOX_PATH}" "${EXPDIR}/cellranger_outs/" --progress
# rclone copy dropbox:${DROPBOX_PATH}/ ${EXPDIR}/cellranger_outs/ --include "_filtered_feature_bc_matrix.h5"
done
```

## Processing data from a new experiment

1. Update experiment config file: add new experiment in `metadata/experiments_config.yaml` - new entry should be called as EXPERIMENT_NAME

2. Download sample metadata 
```bash
# Often easier to first download GWT_sample_metadata.xlsx separately:
# rclone copy gdrive:GWT_perturbseq_analysis/metadata/GWT_sample_metadata.xlsx /oak/stanford/groups/pritch/users/emma/data/GWT/
python process_sample_metadata.py --experiment_name $EXPERIMENT_NAME --datadir /mnt/oak/users/emma/data/GWT/
```

3. Ingest and basic preprocessing of cellranger outputs

```bash
## All at once (works on small experiments)
python preprocess.py --config ../../metadata/experiments_config.yaml --experiment $EXPERIMENT_NAME

## Process each h5 file in parallel
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
for H5_FILE in $(ls ${DATADIR}/${EXPERIMENT_NAME}/cellranger_outs/*/*); do
  BASENAME=$(basename ${H5_FILE} .h5)
  SAMPLE_NAME=$(echo ${BASENAME} | sed 's/_sample_filtered_feature_bc_matrix//')
  LANE_ID=$(basename $(dirname ${H5_FILE}))
  OUTPUT_FILE="${DATADIR}/${EXPERIMENT_NAME}/tmp/${SAMPLE_NAME}_CD4i_R1_Ultima.${LANE_ID}.scRNA.h5ad"
  # Check if output file already exists
  if [ ! -f "${OUTPUT_FILE}" ]; then
    sbatch \
        --partition=pritch \
        --job-name=preprocess_${EXPERIMENT_NAME} \
        --output=$GROUP_SCRATCH/emma/slurm-process_%j.out \
        --error=$GROUP_SCRATCH/emma/slurm-process_%j.err \
        --mem=100G  \
        --time=01:00:00 \
        --wrap="python preprocess.py --experiment ${EXPERIMENT_NAME} --input_h5 ${H5_FILE} --force"
    else 
      echo "${OUTPUT_FILE} exists" 
  fi
done
```

<!-- Merging
```bash
conda activate pertpy-milo
sbatch \
    --partition=pritch \
    --job-name=merge_${EXPERIMENT_NAME} \
    --output=$GROUP_SCRATCH/emma/slurm-merge_%j.out \
    --error=$GROUP_SCRATCH/emma/slurm-merge_%j.err \
    --mem=24G  \
    --time=12:00:00 \
    --wrap="python merge_samples.py"
``` -->

4. Guide RNA assignment 
```bash
# Compute all assignments (for small experiments)
python sgrna_assignment.py $EXPERIMENT_NAME --config ../../metadata/experiments_config.yaml

# Compute assignment for each sample-lane in parallel
conda activate perturb-vs-tissue-env
EXPDIR=/oak/stanford/groups/pritch/users/emma/data/GWT/${EXPERIMENT_NAME}/
for h5ad_file in $(ls $EXPDIR*.sgRNA.h5ad); do
  SAMPLE_NAME=$(basename $h5ad_file .sgRNA.h5ad)
  if [ ! -f "$EXPDIR/${SAMPLE_NAME}.sgrna_assignment.csv" ]; then
        ./submit_sgrna_assignment.sh $EXPERIMENT_NAME $h5ad_file
    else
      echo "Skipping ${SAMPLE_NAME} - assignment file already exists"
  fi
done

# Merge to cell-level assignment for each sample
python sgrna_assignment.py $EXPERIMENT_NAME --merge
```

5. Compute QC stats and exclude low quality cells
```bash
# Quality control filtering
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR1_Psomagen
H5AD_FILES=$(ls ${DATADIR}/${EXPERIMENT_NAME}/tmp/*.scRNA.h5ad)
for f in $H5AD_FILES; do
  SAMPLE_NAME=$(basename ${f} .scRNA.h5ad)
  INPUT_CSV="${DATADIR}/${EXPERIMENT_NAME}/${SAMPLE_NAME}.sgrna_assignment.csv"
  OUTPUT_H5AD="${DATADIR}/${EXPERIMENT_NAME}/tmp/${SAMPLE_NAME}.scRNA.postQC.h5ad"
  
  # if [ ! -f "${OUTPUT_H5AD}" ] && [ -f "${INPUT_CSV}" ]; then
    echo $SAMPLE_NAME
    sbatch \
          --partition=pritch \
          --job-name=qc_${EXPERIMENT_NAME}_${SAMPLE_NAME} \
          --output=$GROUP_SCRATCH/emma/slurm-qc_%j.out \
          --error=$GROUP_SCRATCH/emma/slurm-qc_%j.err \
          --mem=24G  \
          --time=00:30:00 \
          --wrap="python qc_samples.py --experiment_name=${EXPERIMENT_NAME} --sample_id=${SAMPLE_NAME}"
  # fi
done
```

```bash
# Copy notebook to run analysis
cp qc_PilotD2Redo_Lane2.ipynb qc_${EXPERIMENT_NAME}.ipynb

# To convert to report after editing
jupyter nbconvert qc_${EXPERIMENT_NAME}.ipynb --to html
```

6. Estimate effect of each guide (to exclude ineffective from DE analysis) - see `estimate_guide_effect.ipynb`


## Outputs
- `QC_summary_stats.csv` - summary of QC metric statistics for each sample and lane
- `perturbation_counts.csv` - count of number of cells per perturbation for each sample and lane
- `{EXPERIMENT_NAME}.guide_effect.{culture_condition}.csv` - summary stats to assess sgRNA effect on target gene compared to expression of gene in NTC controls
- `no_effect_guides.txt` - guides with no significant effect in any condition

<!-- - `{EXPERIMENT_NAME}_merged.gex.h5ad` - merged raw count matrices (output of `preprocess.py`)
- `{EXPERIMENT_NAME}_merged.gex.lognorm.postQC_obs.csv` - cell-level metadata after QC, with sgRNA assignment and mask for high quality cells (`QC_mask`) (output of `qc_{EXPERIMENT_NAME}.ipynb`)
- `knockdown_efficacy_simple.csv` - coarse estimate of knock-down efficiency for each perturbed gene 
- `guide_ontarget_effect_simple.csv` - estimate of knock-down efficiency for each guide -->
