## Preparing for DE analysis

1. Identify putative ineffective guides (no on-target knockdown with sufficient confidence) - see `src/1_preprocess/estimate_guide_effect.ipynb`

2. Pseudobulk dataset by replicate
```bash
DATADIR=/oak/stanford/groups/pritch/users/emma/data/GWT/
EXPERIMENT_NAME=CRiCD4_Run1_Illumina
python make_pseudobulk.py --experiment_name $EXPERIMENT_NAME --datadir ${DATADIR}${EXPERIMENT_NAME}
```

3. Select common features to test and split perturbations to test into chunks - store parameters in config file
```bash
python prep_DE.py --config ${DATADIR}/DE_results/DE_config.yaml
```

Output files (required for DE analysis scripts):

- Pseudobulk expression, with annotation of samples to keep for DE analysis in `.obs['keep_for_DE']` (`{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad`)
- List of transcriptome genes to test for DE (`{datadir}/DE_test_genes.txt`)
- Assignment of perturbations to processing chunks (`f'{datadir}/DE_target2chunk.csv.gz'`)

## Running DE analysis 

Submit DE analysis for each chunk (with SLURM)
```bash
for c in Rest Stim8hr; do 
    ./submit_DE_chunked.sh $EXPERIMENT_NAME $c
done
```

Merge outputs in AnnData object
```bash
python merge_DE_results.py --experiment_name $EXPERIMENT_NAME
```

Output files (required for DE analysis scripts):
- DE analysis statistics for each perturbation and condition (`{datadir}/DE_results/{experiment_name}.merged_DE_results.h5ad`)