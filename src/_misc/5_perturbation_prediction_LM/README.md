Make pseudobulks
```bash
conda activate perturb-vs-tissue-env
EXPERIMENT_NAME=CRiCD4_Run1_Illumina
sbatch \
    --partition=pritch \
    --job-name=prediction_pbulk_${EXPERIMENT_NAME} \
    --output=$GROUP_SCRATCH/emma/slurm-pbulk.out \
    --error=$GROUP_SCRATCH/emma/slurm-pbulk.err \
    --mem=75G \
    --time=3:00:00 \
    --wrap="python prediction_pseudobulk.py --datadir /oak/stanford/groups/pritch/users/emma/data/GWT/${EXPERIMENT_NAME}/ --experiment_name ${EXPERIMENT_NAME}"

# Compare to Replogle2022
sbatch \
    --partition=pritch \
    --job-name=prediction_pbulk_K562_gwps \
    --output=$GROUP_SCRATCH/emma/slurm-pbulk-k562.out \
    --error=$GROUP_SCRATCH/emma/slurm-pbulk-k562.err \
    --mem=75G \
    --time=3:00:00 \
    --wrap="python prediction_pseudobulk.py --datadir /oak/stanford/groups/pritch/users/emma/data/wesvae-data/ --experiment_name K562_gwps --perturb_col gene"
```