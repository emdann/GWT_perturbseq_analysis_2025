## Preparing for DE analysis

1. Identify putative ineffective guides (no on-target knockdown with sufficient confidence) - see `src/1_preprocess/estimate_guide_effect.ipynb`

2. Pseudobulk dataset by replicate
```bash
conda activate pertpy-milo
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
EXPERIMENT_NAME=CD4iR1_Psomagen

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
    --wrap="python make_pseudobulk.py ${f} --sample_metadata_csv ${DATADIR}/sample_metadata/GWT_sample_metadata.${EXPERIMENT_NAME}.csv"
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
for c in Rest Stim8hr; do 
    ./submit_DE_chunked.sh DE_config_CD4iR1_Psomagen.yaml $c
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