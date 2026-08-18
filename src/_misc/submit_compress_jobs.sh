#!/bin/bash
#
# Submit SLURM jobs to recompress the to_share h5ads (gzip + shuffle).
#
# Only run this after compress_h5ad_for_upload.py has been validated on one
# sample -- each job reads ~110-160 GiB over NFS and writes a new file.
#
# Usage:
#   ./submit_compress_jobs.sh              # all 12 samples
#   ./submit_compress_jobs.sh D1_Rest ...  # named samples only
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHON_SCRIPT="${SCRIPT_DIR}/compress_h5ad_for_upload.py"

# I/O bound, single stream at ~75 MiB/s; ~1.5h/file with headroom.
TIME="6:00:00"
MEM="8G"      # 256 MiB slab buffers + gzip state; 8G is generous
CPUS="2"

LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

if [ $# -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for donor in D1 D2 D3 D4; do
        for condition in Rest Stim8hr Stim48hr; do
            SAMPLES+=("${donor}_${condition}")
        done
    done
fi

echo "Submitting compression jobs for ${#SAMPLES[@]} samples..."
echo "========================================"

for sample in "${SAMPLES[@]}"; do
    job_id=$(sbatch \
        --time="${TIME}" \
        --mem="${MEM}" \
        --cpus-per-task="${CPUS}" \
        --job-name="compress_${sample}" \
        --output="${LOG_DIR}/compress_${sample}.out" \
        --error="${LOG_DIR}/compress_${sample}.err" \
        --wrap="eval \"\$(conda shell.bash hook)\" && conda activate rapids_singlecell && python -u ${PYTHON_SCRIPT} --sample ${sample}" \
        | awk '{print $4}')
    echo "Submitted ${sample}: Job ID ${job_id}"
done

echo "========================================"
echo "Monitor jobs with: squeue -u \$USER"
echo "Check logs in: ${LOG_DIR}"
echo
echo "Each job verifies its own output (CRC32 of X + obs/var equality) and exits"
echo "non-zero on mismatch. Confirm all succeeded before deleting any originals:"
echo "  grep -l SUCCESS ${LOG_DIR}/compress_*.out | wc -l"
