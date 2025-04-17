#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <experiment_name> <culture_condition>"
    echo "Example: $0 experiment1 baseline"
    exit 1
fi

EXPERIMENT_NAME=$1
CONDITION=$2
DATADIR=${3:-"/oak/stanford/groups/pritch/users/emma/data/GWT/"}  # Default path or override with arg

# Check if the chunk file exists
CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_target2chunk.csv.gz"
if [ ! -f "$CHUNK_FILE" ]; then
    echo "Error: Chunk file not found at $CHUNK_FILE"
    exit 1
fi

# Get number of chunks from the header of the gzipped CSV file
N_CHUNKS=$(zcat "$CHUNK_FILE" | head -n 1 | tr ',' '\n' | grep -c "chunk_")
if [ $N_CHUNKS -eq 0 ]; then
    echo "Error: No chunks found in $CHUNK_FILE"
    exit 1
fi

conda activate pertpy-milo

echo "Submitting $N_CHUNKS DE analysis jobs for experiment $EXPERIMENT_NAME, condition $CONDITION"

# Create a job array for all chunks (0-based indexing for chunks)
sbatch \
    --partition=pritch \
    --job-name=DE_${EXPERIMENT_NAME}_${CONDITION} \
    --output=$GROUP_SCRATCH/emma/slurm-DE_%A_%a.out \
    --error=$GROUP_SCRATCH/emma/slurm-DE_%A_%a.err \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=32G \
    --time=2:00:00 \
    --array=0-$((N_CHUNKS-1)) \
    --wrap="python run_DE_chunk.py \
        --datadir $DATADIR \
        --experiment_name $EXPERIMENT_NAME \
        --test_chunk \$SLURM_ARRAY_TASK_ID \
        --culture_condition $CONDITION"