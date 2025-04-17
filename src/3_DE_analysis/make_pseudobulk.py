from MultiStatePerturbSeqDataset import *
import scanpy as sc
import pandas as pd
import numpy as np
import argparse


def make_pseudobulk(datadir, experiment_name, condition_col='culture_condition', sgrna_col='guide_id'):
    """
    Create pseudobulk data from single-cell RNA-seq data.
    
    Parameters:
    -----------
    datadir : str
        Directory containing the data files
    experiment_name : str
        Name of the experiment
    condition_col : str, optional
        Column name for condition information, default is 'culture_condition'
    sgrna_col : str, optional
        Column name for sgRNA information, default is 'guide_id'
    """
    adata = sc.read_h5ad(f'{datadir}/tmp/{experiment_name}_merged.gex.h5ad', backed=True)
    adata.obs = pd.read_csv(f'{datadir}/{experiment_name}_merged.gex.lognorm.postQC_obs.csv', compression='gzip', index_col=0)
    adata = adata[adata.obs['QC_mask']].to_memory()

    adata.obs['perturbed_gene_id'] = np.where(adata.obs['perturbed_gene_name'] == 'NTC', 'NTC', adata.obs['perturbed_gene_id'])

    mspert_data = MultistatePerturbSeqDataset(
        adata,
        state_col=condition_col,
        sample_cols=['donor_id'],
        target_col='perturbed_gene_id',
        perturbation_type='CRISPRi',
        sgrna_col=sgrna_col,
        control_level='NTC'
    )
    mspert_data.pseudobulk()
    mspert_data.save(f'{datadir}/{experiment_name}_merged.DE_pseudobulk.h5ad')
    
    return mspert_data


def main():
    parser = argparse.ArgumentParser(description='Create pseudobulk data from single-cell RNA-seq data')
    parser.add_argument('--datadir', type=str, required=True, help='Directory containing the data files')
    parser.add_argument('--experiment_name', type=str, required=True, help='Name of the experiment')
    parser.add_argument('--condition_col', type=str, default='culture_condition', help='Column name for condition information')
    parser.add_argument('--sgrna_col', type=str, default='guide_id', help='Column name for sgRNA information')
    
    args = parser.parse_args()
    
    make_pseudobulk(
        args.datadir,
        args.experiment_name,
        args.condition_col,
        args.sgrna_col
    )


if __name__ == "__main__":
    main()