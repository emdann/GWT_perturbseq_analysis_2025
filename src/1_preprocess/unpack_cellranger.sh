#!/bin/bash 

if [ "$#" -lt 1 ]; then
   echo "Usage: $0 CELLRANGER_DIR [TARGET_DIR]"
   exit 1
fi

CELLRANGER_DIR=$1
if [ "$#" -eq 2 ]; then
   TARGET_DIR=$2
else
   TARGET_DIR=$(dirname "$CELLRANGER_DIR")/TcellsGW_$(basename "$CELLRANGER_DIR")
fi

mkdir -p "${TARGET_DIR}/cellranger_outs"

for dir in ${CELLRANGER_DIR}/outs/per_sample_outs/*; do
   dirname=$(basename $dir)
   cp "${dir}/count/sample_filtered_feature_bc_matrix.h5" "${TARGET_DIR}/cellranger_outs/${dirname}_sample_filtered_feature_bc_matrix.h5"
done