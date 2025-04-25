## Preparing for DE analysis

1. Identify putative ineffective guides (no on-target knockdown with sufficient confidence) - see `src/1_preprocess/estimate_guide_effect.ipynb`

2. Pseudobulk dataset by replicate
```bash
# DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
DATADIR=/mnt/oak/users/emma/data/GWT/
EXPERIMENT_NAME=CRiCD4_Run1_Illumina
python make_pseudobulk.py --experiment_name $EXPERIMENT_NAME --datadir ${DATADIR}${EXPERIMENT_NAME}
```

3. Select common features to test and split perturbations to test into chunks - store parameters in config file
```bash
python prep_DE.py --config ${DATADIR}${EXPERIMENT_NAME}/DE_results/DE_config.yaml
```

Output files (required for DE analysis scripts):

- Pseudobulk expression, with annotation of samples to keep for DE analysis in `.obs['keep_for_DE']` (`{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad`)
- List of transcriptome genes to test for DE (`{datadir}/DE_test_genes.txt`)
- Assignment of perturbations to processing chunks (`f'{datadir}/DE_target2chunk.csv.gz'`)

## Running DE analysis 

Submit DE analysis for each chunk (with SLURM)
```bash
for c in Rest Stim8hr; do 
    # ./submit_DE_chunked.sh DE_config.yaml $c
    ./submit_DE_chunked.sh DE_config_nodonor.yaml $c
done
```

Run MASH for groups of chunks
```bash
CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_target2chunk.csv.gz"
N_CHUNKS=$(zcat "$CHUNK_FILE" | head -n 1 | tr ',' '\n' | grep -c "chunk_")
N_GROUPS=$(( (N_CHUNKS + 4) / 5 ))  # Integer division with ceiling

## Locally
for i in $(seq 0 $((${N_GROUPS}-1))); do
    python run_MASH_chunk.py --datadir $DATADIR --experiment_name $EXPERIMENT_NAME --test_group $i
done

## With SLURM
sbatch \
    --partition=pritch \
    --job-name=MASH_${EXPERIMENT_NAME} \
    --output=$GROUP_SCRATCH/emma/slurm-MASH_%A_%a.out \
    --error=$GROUP_SCRATCH/emma/slurm-MASH_%A_%a.err \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=3 \
    --mem=50G \
    --time=0:45:00 \
    --array=0-$((N_GROUPS-1)) \
    --wrap="python run_MASH_chunk.py \
        --config DE_config.yaml \
        --test_group \$SLURM_ARRAY_TASK_ID"
```

Merge outputs in AnnData object
```bash
python merge_DE_results.py --config DE_config.yaml
```

Output files (required for DE analysis scripts):
- DE analysis statistics for each perturbation and condition (`{datadir}/DE_results_{run_name}/{experiment_name}.merged_DE_results.h5ad`)

## Analysing DE results

- `DE_results_analysis.ipynb` - exploratory analysis of DE results 
- `FACS_comparison.ipynb` - comparison of DE results with FACS screens from Marson Lab