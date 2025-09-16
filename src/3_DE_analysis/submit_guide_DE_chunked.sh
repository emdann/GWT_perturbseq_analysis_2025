#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <experiment_name> <culture_condition>"
    echo "Example: $0 experiment1 baseline"
    exit 1
fi

CONFIG=$1
CONDITION=$2
DATADIR=${3:-"/mnt/oak/users/emma/data/GWT/"}  # Default path or override with arg

# Extract experiment_name and datadir from the YAML config file
EXPERIMENT_NAME=$(grep "experiment_name:" "$CONFIG" | sed 's/.*experiment_name:[[:space:]]*//' | sed 's/[[:space:]]*$//' | tr -d '"')
DATADIR_CONFIG=$(grep "datadir:" "$CONFIG" | sed 's/.*datadir:[[:space:]]*//' | sed 's/[[:space:]]*$//' | tr -d '"')
if [ -z "$EXPERIMENT_NAME" ]; then
    echo "Error: Could not extract experiment_name from config file"
    exit 1
fi
if [ -n "$DATADIR_CONFIG" ]; then
    DATADIR=$DATADIR_CONFIG
fi

# Check if the chunk file exists
CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_single_guide_target2chunk.${CONDITION}.csv.gz"
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

# Get the run_name from YAML config to construct the correct results directory
RUN_NAME=$(grep "run_name:" "$CONFIG" | sed 's/.*run_name:[[:space:]]*//' | sed 's/[[:space:]]*$//' | tr -d '"')
if [ -z "$RUN_NAME" ]; then
    RUN_NAME="default"
fi
RESULTS_DIR="${DATADIR}/${EXPERIMENT_NAME}/DE_results_${RUN_NAME}/tmp"


echo "Checking for existing results in: $RESULTS_DIR"
echo "Total chunks to process: $N_CHUNKS (0 to $((N_CHUNKS-1)))"

# Find missing chunks
MISSING_CHUNKS=()
for chunk_id in $(seq 0 $((N_CHUNKS-1))); do
    OUTPUT_FILE="${RESULTS_DIR}/DE_results_by_guide.${CONDITION}.chunk_${chunk_id}.csv.gz"
    if [ ! -f "$OUTPUT_FILE" ]; then
        MISSING_CHUNKS+=($chunk_id)
    fi
done

# Check if any chunks are missing
if [ ${#MISSING_CHUNKS[@]} -eq 0 ]; then
    echo "All chunks are already completed for condition: $CONDITION"
    echo "No jobs to submit."
    exit 0
fi

echo "Found ${#MISSING_CHUNKS[@]} missing chunks: ${MISSING_CHUNKS[*]}"

conda activate pertpy-milo

echo "Submitting ${#MISSING_CHUNKS[@]} individual DE analysis jobs for experiment $EXPERIMENT_NAME, condition $CONDITION"
echo "Missing chunks: ${MISSING_CHUNKS[*]}"

# Submit individual jobs for each missing chunk
for chunk_id in "${MISSING_CHUNKS[@]}"; do
    echo "Submitting job for chunk $chunk_id"
    sbatch \
        --job-name=DE_guide_${EXPERIMENT_NAME}_${CONDITION}_chunk${chunk_id} \
        --output=./slurm-DE_guide_${EXPERIMENT_NAME}_${CONDITION}_chunk${chunk_id}.out \
        --error=./slurm-DE_guide_${EXPERIMENT_NAME}_${CONDITION}_chunk${chunk_id}.err \
        --mem=75G \
        --time=2:00:00 \
        --wrap="python run_guide_DE_chunk.py \
                    --config $CONFIG \
                    --test_chunk $chunk_id \
                    --culture_condition $CONDITION"
done

echo "Submitted ${#MISSING_CHUNKS[@]} individual jobs for chunks: ${MISSING_CHUNKS[*]}"

