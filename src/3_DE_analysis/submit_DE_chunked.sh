#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <experiment_name> <culture_condition>"
    echo "Example: $0 experiment1 baseline"
    exit 1
fi

CONFIG=$1
CONDITION=$2
DATADIR=${3:-"/oak/stanford/groups/pritch/users/emma/data/GWT/"}  # Default path or override with arg

# Extract experiment_name and datadir from the config file
EXPERIMENT_NAME=$(grep "experiment_name:" "$CONFIG" | cut -d ":" -f2 | tr -d ' "')
DATADIR_CONFIG=$(grep "datadir:" "$CONFIG" | cut -d ":" -f2 | tr -d ' "')
if [ -z "$EXPERIMENT_NAME" ]; then
    echo "Error: Could not extract experiment_name from config file"
    exit 1
fi
if [ -n "$DATADIR_CONFIG" ]; then
    DATADIR=$DATADIR_CONFIG
fi

# Check if the chunk file exists
CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_target2chunk.${CONDITION}.csv.gz"
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
    --cpus-per-task=3 \
    --mem=50G \
    --time=2:00:00 \
    --array=0-$((N_CHUNKS-1)) \
    --wrap="python run_DE_chunk.py \
        --config $CONFIG \
        --test_chunk \$SLURM_ARRAY_TASK_ID \
        --culture_condition $CONDITION"

# # Create output directory if it doesn't exist
# OUTPUT_DIR="${DATADIR}/${EXPERIMENT_NAME}/DE_results_all_confounders/tmp/"
# # mkdir -p "$OUTPUT_DIR"

# # Find missing chunks
# CONDITION=Rest
# # Check if the chunk file exists
# CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_target2chunk.${CONDITION}.csv.gz"
# if [ ! -f "$CHUNK_FILE" ]; then
#     echo "Error: Chunk file not found at $CHUNK_FILE"
#     exit 1
# fi

# # Get number of chunks from the header of the gzipped CSV file
# N_CHUNKS=$(zcat "$CHUNK_FILE" | head -n 1 | tr ',' '\n' | grep -c "chunk_")
# if [ $N_CHUNKS -eq 0 ]; then
#     echo "Error: No chunks found in $CHUNK_FILE"
#     exit 1
# fi

# MISSING_CHUNKS=()
# for ((i=0; i<N_CHUNKS; i++)); do
#     OUTPUT_FILE="${OUTPUT_DIR}/DE_results.${CONDITION}.chunk_${i}.csv.gz"
#     if [ ! -f "$OUTPUT_FILE" ]; then
#         echo "Missing output file for chunk $i: $OUTPUT_FILE"
#         MISSING_CHUNKS+=("$i")
#     fi
# done

# for i in "${MISSING_CHUNKS[@]}"; do 
#     sbatch \
#     --partition=pritch \
#     --job-name=DE_${EXPERIMENT_NAME}_${CONDITION}_missingchunk \
#     --output=$GROUP_SCRATCH/emma/slurm-DE_%j_missingchunk.out \
#     --error=$GROUP_SCRATCH/emma/slurm-DE_%j_missingchunk.err \
#     --nodes=1 \
#     --ntasks=1 \
#     --cpus-per-task=3 \
#     --mem=50G \
#     --time=3:00:00 \
#     --wrap="python run_DE_chunk.py \
#         --config $CONFIG \
#         --test_chunk ${i} \
#         --culture_condition $CONDITION"
#     done