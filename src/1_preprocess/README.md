## Moving cellranger outs from Dropbox to Sherlock

```bash
# Setup folders in Sherlock
EXPERIMENT_NAME=CRiCD4_Run1_Illumina
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

3. Preprocess and merge cellranger outputs

```bash
python preprocess.py --config ../../metadata/experiments_config.yaml --experiment $EXPERIMENT_NAME
```

4. Guide RNA assignment 
```bash
# Compute all assignments 
python sgrna_assignment.py $EXPERIMENT_NAME --config ../../metadata/experiments_config.yaml

# To submit with SLURM (parallelize)
conda activate perturb-vs-tissue-env
EXPDIR=/oak/stanford/groups/pritch/users/emma/data/GWT/${EXPERIMENT_NAME}/
for h5ad_file in $(ls $EXPDIR*.sgRNA.h5ad); do
  ./submit_sgrna_assignment.sh $EXPERIMENT_NAME $h5ad_file
done 

# Merge to cell-level assignment for each sample
python sgrna_assignment.py $EXPERIMENT_NAME --merge

```

5. Basic QC plots and analysis
```bash
# Copy notebook to run analysis
cp qc_PilotD2Redo_Lane2.ipynb qc_${EXPERIMENT_NAME}.ipynb

# To convert to report after editing
jupyter nbconvert qc_${EXPERIMENT_NAME}.ipynb --to html
```

## Outputs

- `{EXPERIMENT_NAME}_merged.gex.lognorm.h5ad` - merged count matrices (output of `preprocess.py`)
- `{EXPERIMENT_NAME}_merged.gex.lognorm.postQC.h5ad` - merged count matrices after filtering of low quality cells, clustering and sgRNA assignment (output of `qc_{EXPERIMENT_NAME}.ipynb`)
- `knockdown_efficacy_simple.csv` - coarse estimate of knock-down efficiency for each perturbed gene 
- `sgRNA_assignment.csv` - estimate of knock-down efficiency for each perturbed gene 

## Sharing outputs on Dropbox

From Comino:
```bash
DROPBOX_PATH=GRNPerturbSeq/3_expts/processed_data/
rclone mkdir dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/
for f in ${EXPERIMENT_NAME}_merged.gex.lognorm.postQC.h5ad ${EXPERIMENT_NAME}_merged.gex.counts.postQC.h5ad knockdown_efficacy_simple.csv; do 
    rclone copy ${EXPDIR}/${f} dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/ --checksum --ignore-times
    done
```