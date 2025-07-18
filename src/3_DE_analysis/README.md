## Preparing for DE analysis

1. Identify putative ineffective guides (no on-target knockdown with sufficient confidence) - see `src/1_preprocess/estimate_guide_effect.ipynb`

2. Pseudobulk dataset by replicate
```bash
conda activate pertpy-milo
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR2_Psomagen

H5AD_FILES=$(ls $DATADIR/$EXPERIMENT_NAME/tmp/*.postQC.h5ad)
for f in $H5AD_FILES; do
sbatch \
    --partition=pritch \
    --job-name=pbulk_${f} \
    --output=$GROUP_SCRATCH/emma/slurm-pbulk_%j.out \
    --error=$GROUP_SCRATCH/emma/slurm-pbulk_%j.err \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=20G \
    --time=0:30:00 \
    --wrap="python make_pseudobulk.py aggregate ${f} --sample_metadata_csv ${DATADIR}/sample_metadata/GWT_sample_metadata.${EXPERIMENT_NAME}.csv"
done
```

Merge pseudobulks for each sample
```bash
for SAMPLE_ID in CD4i_R2_D4_Stim8hr_CD4i_R2_Ultima; do
    sbatch \
        --partition=pritch \
        --job-name=pbulk_merge_${s} \
        --output=$GROUP_SCRATCH/emma/slurm-pbulk_%j.out \
        --error=$GROUP_SCRATCH/emma/slurm-pbulk_%j.err \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=20G \
        --time=0:30:00 \
        --wrap="python make_pseudobulk.py merge $SAMPLE_ID"
done
```

3. Select common features to test and split perturbations to test into chunks - store parameters in config file
```bash
prep_DE.ipynb 
```

Output files (required for DE analysis scripts):

- Pseudobulk expression, with annotation of samples to keep for DE analysis in `.obs['keep_for_DE']` (`{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad`)
- List of transcriptome genes to test for DE (`{datadir}/DE_test_genes.{condition}.txt`)
- Assignment of perturbations to processing chunks (`f'{datadir}/DE_target2chunk.{condition}csv.gz'`)

## Running DE analysis 

Submit DE analysis for each chunk (with SLURM)
```bash
for c in Rest Stim8hr Stim48hr; do 
    ./submit_DE_chunked.sh DE_config_full.yaml $c
done

# Temporary
conda activate rpy2-voodoo
CONFIG=DE_config_full.yaml
CONDITION=Stim8hr

# Extract experiment_name and datadir from the config file
EXPERIMENT_NAME=$(grep "experiment_name:" "$CONFIG" | cut -d ":" -f2 | tr -d ' "')
DATADIR=/mnt/oak/users/emma/data/GWT/
RESULTS_DIR=/mnt/oak/users/emma/data/GWT/CD4i_final/DE_results_all_confounders/tmp/
CHUNK_FILE="${DATADIR}/${EXPERIMENT_NAME}/DE_target2chunk.${CONDITION}.csv.gz"
N_CHUNKS=$(zcat "$CHUNK_FILE" | head -n 1 | tr ',' '\n' | grep -c "chunk_")

# Find missing chunks
MISSING_CHUNKS=()
for chunk_id in $(seq 0 $((N_CHUNKS-1))); do
    OUTPUT_FILE="${RESULTS_DIR}/DE_results.${CONDITION}.chunk_${chunk_id}.csv.gz"
    if [ ! -f "$OUTPUT_FILE" ]; then
        MISSING_CHUNKS+=($chunk_id)
    fi
done

# Process missing chunks in groups of 3
for ((i=0; i<${#MISSING_CHUNKS[@]}; i+=3)); do
    # Get up to 3 chunks starting at index i
    CHUNK_GROUP=()
    for ((j=i; j<i+3 && j<${#MISSING_CHUNKS[@]}; j++)); do
        CHUNK_GROUP+=(${MISSING_CHUNKS[j]})
    done
    
    # Convert chunk group to comma-separated string
    ARRAY_SPEC=$(IFS=,; echo "${CHUNK_GROUP[*]}")
    echo $ARRAY_SPEC

    sbatch --mem=75G \
        --output=slurm-DE_%j.out \
        --error=slurm-DE_%j.err \
        --wrap="python run_DE_chunk.py \
            --config $CONFIG \
            --test_chunk $ARRAY_SPEC \
            --culture_condition $CONDITION \
            --n_cpus 3"
done

```bash
for c in Rest Stim8hr Stim48hr; do 
    ./submit_guide_DE_chunked.sh DE_config_by_guide.yaml $c
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
python merge_DE_results.py --config DE_config_CD4iR1_Psomagen.yaml
```

## Analysing DE results

- `DE_results_analysis.ipynb` - exploratory analysis of DE results 
- `FACS_comparison.ipynb` - comparison of DE results with FACS screens from Marson Lab

## Output files 

-  `{experiment_name}_merged.DE_pseudobulk.h5ad` - Pseudobulked data object (sum of counts across donor-condition-guide)
- `{datadir}/DE_results_{run_name}/{experiment_name}.merged_DE_results.h5ad` - DE analysis statistics for each perturbation and condition
- `{datadir}/DE_results_{run_name}/DE_summary_stats_per_target.csv` - Summary of on-target effects and overall effect for each perturbation and condition