## Processing data from a new experiment

1. Update experiment config file: add new experiment in `metadata/experiments_config.yaml`

2. Download sample metadata 
```bash
# Often easier to first download GWT_sample_metadata.xlsx separately:
# rclone copy gdrive:GWT_perturbseq_analysis/metadata/GWT_sample_metadata.xlsx /oak/stanford/groups/pritch/users/emma/data/GWT/

EXPERIMENT_NAME=PilotD2Redo_Lane2
python process_sample_metadata.py --experiment_name $EXPERIMENT_NAME --datadir /mnt/oak/users/emma/data/GWT/
```

3. Preprocess and merge cellranger outputs
```bash
EXPERIMENT_NAME=PilotD2Redo_Lane2
DATADIR=/mnt/oak/users/emma/data/GWT/
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/

# Store h5 files in cellranger_outs/ directory
mkdir ${EXPDIR}/cellranger_outs/
mv /path/to/cellranger_outputs/*.h5 $EXPDIR/cellranger_outs/

# Run preprocessing
python preprocess.py --config ../../metadata/experiments_config.yaml --experiment $EXPERIMENT_NAME
```

4. Guide RNA assignment 
```bash
python sgrna_assignment.py --config ../../metadata/experiments_config.yaml --experiment $EXPERIMENT_NAME
```

5. Basic QC plots and analysis
```bash
# Copy notebook to run analysis
EXPERIMENT_NAME=PilotD2Redo_Lane2
cp qc_PilotD2Redo.ipynb qc_${EXPERIMENT_NAME}.ipynb

# To convert to report after editing
jupyter nbconvert qc_${EXPERIMENT_NAME}.ipynb --to html
```