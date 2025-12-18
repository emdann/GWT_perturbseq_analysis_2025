#!/bin/bash
#
# Submit SLURM jobs to merge technical replicates for all biological samples
#

# Set up paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHON_SCRIPT="${SCRIPT_DIR}/merge_sample_for_upload.py"

# SLURM parameters
TIME="8:00:00"    
MEM="64G"         
CPUS="2"          

# Log directory
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Define all biological samples (from notebook output: 8 total samples)
# Donors: D1, D2, D3, D4
# Conditions: Rest, Stim8hr
DONORS=("D1" "D2" "D3" "D4")
CONDITIONS=("Rest" "Stim8hr" "Stim48hr")

# Submit jobs for all combinations
echo "Submitting merge jobs for all biological samples..."
echo "========================================"

for donor in "${DONORS[@]}"; do
    for condition in "${CONDITIONS[@]}"; do
        bio_sample="${donor}_${condition}"

        # Define log files
        log_out="${LOG_DIR}/${bio_sample}.out"
        log_err="${LOG_DIR}/${bio_sample}.err"

        # Submit job
        job_id=$(sbatch \
            --time="${TIME}" \
            --mem="${MEM}" \
            --cpus-per-task="${CPUS}" \
            --job-name="merge_${bio_sample}" \
            --output="${log_out}" \
            --error="${log_err}" \
            --wrap="eval \"\$(conda shell.bash hook)\" && conda activate rapids_singlecell && python -u ${PYTHON_SCRIPT} --donor ${donor} --condition ${condition}" \
            | awk '{print $4}')

        echo "Submitted ${bio_sample}: Job ID ${job_id}"
    done
done

echo "========================================"
echo "All jobs submitted!"
echo "Monitor jobs with: squeue -u \$USER"
echo "Check logs in: ${LOG_DIR}"
