#!/bin/bash

# Create logs directory if it doesn't exist
mkdir -p logs

# Data directory
datadir="/mnt/oak/users/emma/data/GWT/CD4i_final/"
results_dir="${datadir}donor_robustness_analysis/"

# Find all target files
target_files=($(ls ${results_dir}selected_targets_robustness_*.txt))

# Loop through each target file
for target_file in "${target_files[@]}"; do
    # Extract condition from filename
    condition=$(basename "$target_file" | sed 's/selected_targets_robustness_\(.*\)\.txt/\1/')
    
    # Skip the "all" condition file if it exists
    if [[ "$condition" == "all" ]]; then
        continue
    fi
    
    # Read targets from file and join with commas
    targets=$(paste -sd "," "$target_file")
    
    echo "Submitting job for condition: $condition"
    
    sbatch --job-name="robustness_${condition}" \
           --time=6:00:00 \
           --mem=50G \
           --cpus-per-task=4 \
           --output="logs/robustness_${condition}_%j.out" \
           --error="logs/robustness_${condition}_%j.err" \
           --wrap="python get_donor_robustness.py --condition $condition --target $targets"
done

echo "All jobs submitted!"