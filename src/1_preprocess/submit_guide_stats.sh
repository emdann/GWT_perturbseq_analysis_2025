#!/bin/bash
#SBATCH --job-name=guide_stats
#SBATCH --array=1-284
#SBATCH --time=02:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/guide_stats_%A_%a.out
#SBATCH --error=logs/guide_stats_%A_%a.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Combine sample lists
cd /mnt/oak/users/emma/bin/GWT_perturbseq_analysis/src/1_preprocess

# Create combined sample list with experiment names
{
    # CD4iR1_Psomagen samples (92 samples)
    while read sample_id; do
        echo "CD4iR1_Psomagen,$sample_id"
    done < CD4iR1_samples.txt

    # CD4iR2_Psomagen samples (192 samples)
    while read sample_id; do
        echo "CD4iR2_Psomagen,$sample_id"
    done < CD4iR2_samples.txt
} > all_samples.txt

# Get the current sample info based on array index
SAMPLE_INFO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all_samples.txt)
EXPERIMENT_NAME=$(echo $SAMPLE_INFO | cut -d',' -f1)
SAMPLE_ID=$(echo $SAMPLE_INFO | cut -d',' -f2)

echo "Processing experiment: $EXPERIMENT_NAME, sample: $SAMPLE_ID"

# Activate conda environment (adjust path as needed)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rapids_singlecell

# Run the guide group stats computation
python qc_samples.py \
    --experiment_name $EXPERIMENT_NAME \
    --sample_id $SAMPLE_ID \
    --guide_group_stats

echo "Completed processing $SAMPLE_ID"