'''
Sync analysis outputs and intermediate files to Dropbox
'''
DATADIR=/mnt/oak/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR2_Psomagen
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/
DROPBOX_PATH=GRNPerturbSeq/3_expts/processed_data/

rclone mkdir dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/

# Sync raw outputs
LOCAL_DATA_DIR="$DATADIR/czi-psomagen"
find "$LOCAL_DATA_DIR" -name "web_summary.html" -type f | while read file; do
    relative_path=${file#$LOCAL_DATA_DIR/}
    relative_dir=$(dirname "$relative_path")
    rclone copy "$file" "dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/${relative_dir}/" --checksum --ignore-times
done

# Copy AnnData objects x sample-lane
H5AD_FILES=$(ls ${DATADIR}/${EXPERIMENT_NAME}/tmp/*.scRNA.postQC.h5ad)
rclone sync ${DATADIR}/${EXPERIMENT_NAME}/tmp/ dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/ --include "*.scRNA.postQC.h5ad" --checksum -v


# Summary files
FILES_TO_COPY=(
    # "${EXPERIMENT_NAME}_merged.gex.lognorm.postQC_obs.csv"  
    "DE_results_all_confounders/DE_summary_stats_per_target.csv"
    # "QC_summary_stats.csv"
    # "perturbation_counts.csv"
    # "${EXPERIMENT_NAME}.guide_effect.Stim8hr.csv"
    # "${EXPERIMENT_NAME}.guide_effect.Stim48hr.csv"
    # "${EXPERIMENT_NAME}.guide_effect.Rest.csv"
    # "no_effect_guides.txt"
    # "${EXPERIMENT_NAME}_merged.DE_pseudobulk.h5ad"
    # "DE_results_all_confounders/${EXPERIMENT_NAME}.merged_DE_results.h5ad"
)

# Copy each file to Dropbox
for f in "${FILES_TO_COPY[@]}"; do 
    rclone sync ${EXPDIR} dropbox:${DROPBOX_PATH}${EXPERIMENT_NAME}/ --include "${f}" --checksum --ignore-times -v
done

## Sharing raw counts 
DATADIR=/mnt/oak/users/emma/data/GWT/
DROPBOX_PATH=GWT-counts-Christina

EXPERIMENT_NAME=CD4iR1_Psomagen
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/
rclone sync ${DATADIR}/${EXPERIMENT_NAME}/tmp/ dropbox:${DROPBOX_PATH}/ --include "*.scRNA.h5ad" --checksum -v
rclone sync ${DATADIR}/${EXPERIMENT_NAME}/ dropbox:${DROPBOX_PATH}/ --include "*.sgrna_assignment.csv" --checksum -v

EXPERIMENT_NAME=CD4iR2_Psomagen
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/
rclone sync ${DATADIR}/${EXPERIMENT_NAME}/tmp/ dropbox:${DROPBOX_PATH}/ --include "*.scRNA.postQC.h5ad" --checksum -v
rclone sync ${DATADIR}/${EXPERIMENT_NAME}/ dropbox:${DROPBOX_PATH}/ --include "*.sgrna_assignment.csv" --checksum -v


DATADIR=/mnt/oak/users/emma/data/GWT/
EXPERIMENT_NAME=CD4i_final
EXPDIR=${DATADIR}/${EXPERIMENT_NAME}/
DROPBOX_PATH=GRNPerturbSeq/3_expts/processed_data/

rclone sync dropbox:${DROPBOX_PATH}/analysis_largefiles/ ${DATADIR}/${EXPERIMENT_NAME}/ --include "nde50ntotal100_varfiltered_clustering_downstream_genes.csv" --checksum -v

## Reload from backup
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
DROPBOX_PATH=GRNPerturbSeq/3_expts/processed_data/CD4i_final/DE_results_all_confounders/
rclone copy "dropbox:${DROPBOX_PATH}" "${DATADIR}/CD4i_final/DE_results_all_confounders/" \
  --checksum -v
