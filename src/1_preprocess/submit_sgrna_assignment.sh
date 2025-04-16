#!/bin/bash

if [ -z "$1" ]; then
    echo "Error: Experiment name must be provided as first argument"
    echo "Usage: $0 <experiment_name>"
    exit 1
fi

EXPERIMENT_NAME=$1
H5AD_FILE=$2
SAMPLE=$(basename $H5AD_FILE .sgRNA.h5ad)
CHUNK_SIZE=200                          # Number of sgRNAs per chunk
MAX_CORES_PER_NODE=10                   # Maximum CPUs available on your nodes
TOTAL_GUIDES=22000                      # Total number of guides

# Calculate the number of chunks
N_CHUNKS=$(( ($TOTAL_GUIDES + $CHUNK_SIZE - 1) / $CHUNK_SIZE ))

# Submit a job array for all chunks
echo "${EXPERIMENT_NAME} - ${H5AD_FILE}"
echo "Submitting job array for ${N_CHUNKS} chunks..."

# Submit the job array
sbatch \
    --partition=pritch \
    --job-name=sgrna_${SAMPLE}_array \
    --output=$GROUP_SCRATCH/emma/slurm-sgrna_${SAMPLE}_chunk%a_%j.out \
    --error=$GROUP_SCRATCH/emma/slurm-sgrna_${SAMPLE}_chunk%a_%j.err \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=$MAX_CORES_PER_NODE \
    --mem=24G \
    --time=1:00:00 \
    --array=0-$(($N_CHUNKS-1)) \
    --export=ALL,EXPERIMENT_NAME=$EXPERIMENT_NAME,H5AD_FILE=$H5AD_FILE,CHUNK_SIZE=$CHUNK_SIZE,TOTAL_GUIDES=$TOTAL_GUIDES,MAX_CORES_PER_NODE=$MAX_CORES_PER_NODE \
    --wrap="i=\$SLURM_ARRAY_TASK_ID; \
           start_idx=\$(( i * CHUNK_SIZE )); \
           end_idx=\$(( (i + 1) * CHUNK_SIZE )); \
           if (( end_idx > TOTAL_GUIDES )); then \
               end_idx=\$TOTAL_GUIDES; \
           fi; \
           echo \"Processing chunk \$i (guides \$start_idx-\$end_idx)...\"; \
           python sgrna_assignment.py \$EXPERIMENT_NAME --crispr_h5ad \$H5AD_FILE --n_cores \$MAX_CORES_PER_NODE --n_guides_parallel \$CHUNK_SIZE --start_idx \$start_idx --end_idx \$end_idx --chunk_id \$i"

echo "All jobs submitted."