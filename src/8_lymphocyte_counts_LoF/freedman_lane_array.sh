conda activate gwt

# Parameters
PERMUTATIONS=1000
SEED=42
MIN_TARGETS=20

TRAIT=lymph
CONDS=Stim8hr,Stim48hr,Stim8hr+Stim48hr,K562,K562+Stim8hr,K562+Stim48hr
BASE_COLOR_MAP=Stim8hr:mediumseagreen,Stim48hr:green,Stim8hr+Stim48hr:darkblue,K562:lightgray,K562+Stim8hr:lightblue,K562+Stim48hr:steelblue

# Directories
ITERATION_DIR=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/permutation_iterations
FINAL_OUTDIR=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/permutation_results

# Input files
COND_MATS_FILE=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/cond_mats.pkl
ENSG2SYM_FILE=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/ensg2sym.pkl
SYM2ENSG_FILE=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/sym2ensg.pkl
LOF_FILE=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/lof.tsv
SHET_FILE=/oak/stanford/groups/pritch/users/rgoto/data/GWT_causal_gene_discovery/CD4i/burden_reg_corr/${TRAIT}/shet.tsv

# Create directories
mkdir -p ${ITERATION_DIR}
mkdir -p ${FINAL_OUTDIR}

echo "=== Freedman-Lane Permutation Analysis ==="
echo "Trait: ${TRAIT}"
echo "Conditions: ${CONDS}"
echo "Permutations: ${PERMUTATIONS}"
echo "Iteration directory: ${ITERATION_DIR}"
echo "Final output directory: ${FINAL_OUTDIR}"
echo "=========================================="

# Step 1: Submit array job for individual iterations
echo "Submitting array job for ${PERMUTATIONS} iterations..."
JOB_ID=$(sbatch \
    --partition=pritch \
    --job-name=freedman_lane_iter_${TRAIT} \
    --output=$SCRATCH/slurm-burd/freedman_lane_iter_%A_%a.out \
    --error=$SCRATCH/slurm-burd/freedman_lane_iter_%A_%a.err \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=16G \
    --array=0-$((PERMUTATIONS-1))%400 \
    --wrap="python /oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/freedman_lane_single_iter.py \
            --cond_mats_file ${COND_MATS_FILE} \
            --ensg2sym_file ${ENSG2SYM_FILE} \
            --sym2ensg_file ${SYM2ENSG_FILE} \
            --lof_file ${LOF_FILE} \
            --shet_file ${SHET_FILE} \
            --conds ${CONDS} \
            --output_dir ${ITERATION_DIR} \
            --seed ${SEED} \
            --min_targets ${MIN_TARGETS} \
            --iteration_idx \$SLURM_ARRAY_TASK_ID" | grep -o '[0-9]\+')

echo "Submitted array job with ID: ${JOB_ID}"

# Step 2: Submit integration job that waits for array job to complete
echo "Submitting integration job..."
INTEGRATION_JOB_ID=$(sbatch \
    --partition=pritch \
    --job-name=freedman_lane_integrate_${TRAIT} \
    --output=$SCRATCH/slurm-burd/freedman_lane_integrate_%j.out \
    --error=$SCRATCH/slurm-burd/freedman_lane_integrate_%j.err \
    --time=1:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=32G \
    --dependency=afterok:${JOB_ID} \
    --wrap="python /oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/freedman_lane_integrate.py \
            --output_dir ${ITERATION_DIR} \
            --base_color_map ${BASE_COLOR_MAP} \
            --permutations ${PERMUTATIONS} \
            --final_output_dir ${FINAL_OUTDIR}" | grep -o '[0-9]\+')

echo "Submitted integration job with ID: ${INTEGRATION_JOB_ID}"
echo ""
echo "Job submission complete!"
echo "Array job ID: ${JOB_ID} (${PERMUTATIONS} iterations)"
echo "Integration job ID: ${INTEGRATION_JOB_ID} (waits for array job)"
echo ""
echo "Monitor progress with:"
echo "  squeue -j ${JOB_ID},${INTEGRATION_JOB_ID}"
echo "  tail -f $SCRATCH/slurm-burd/freedman_lane_iter_${JOB_ID}_*.out"
echo "  tail -f $SCRATCH/slurm-burd/freedman_lane_integrate_${INTEGRATION_JOB_ID}.out"
