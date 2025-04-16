'''
Script to train contrastiveVI model on perturb-seq datasets.
Saves model and embeddings for downstream analysis and visualization.

To run on comino:
sbatch --gres=gpu:1 --mem=100G -o slurm-%j.out -e slurm-%j.err --wrap="python train_contrastiveVI.py --experiment CRiCD4_Run1_Illumina --testing"
'''

import os
import scvi
import torch
import anndata
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import rapids_singlecell as rsc 
import scipy.sparse
import mudata
import genomic_features as gf
from typing import Optional, Union, List

device = torch.device("cuda")


def _convert_oak_path(path):
        """Helper function to convert oak paths between different mount points"""
        if not os.path.exists(path):
            return path.replace('/oak/stanford/groups/pritch/', '/mnt/oak/')
        return path

# ---- Feature selection ---- # 
def get_ribo_mito_genes():
    '''Get ribosomal and mitochondrial genes to exclude.'''
    ensdb = gf.ensembl.annotation(species="Hsapiens", version="108")
    genes = ensdb.genes()

    ribo_gene_ids = genes[genes.description.str.startswith('ribosomal protein') & (genes.gene_biotype == 'protein_coding')].gene_id.values
    mito_gene_ids = genes[genes.seq_name == 'MT'].gene_id.values
    mito_ribo_gene_ids = genes[genes.description.str.startswith('mitochondrial ribosomal protein') & (genes.gene_biotype == 'protein_coding')].gene_id.values
    return ribo_gene_ids, mito_gene_ids, mito_ribo_gene_ids

def feature_selection(
    adata: anndata.AnnData, 
    n_hvgs: int = 5000, 
    filter_ribo_mito: bool = True, 
    subset_adata: bool = True, 
    subset_obs: str = 'integration_sample_id',
    ):
    '''Save table of selected features to use for integration.'''
    filter_genes = []
    if filter_ribo_mito:
        ribo_gene_ids, mito_gene_ids, mito_ribo_gene_ids = get_ribo_mito_genes()
        filter_genes.extend(ribo_gene_ids.tolist())
        filter_genes.extend(mito_gene_ids.tolist())
        filter_genes.extend(mito_ribo_gene_ids.tolist())

    if subset_adata:
        adata = sc.pp.sample(adata, fraction=0.1, copy=True)
    
    rsc.get.anndata_to_GPU(adata)
    rsc.pp.calculate_qc_metrics(adata)

    # Filter out highly and lowly expressed outlier genes
    highly_expressed_outiers = adata.var_names[(adata.var['mean_counts'] > 15) & (adata.var['pct_dropout_by_counts'] < 5)].tolist()
    lowly_expressed_outliers = adata.var_names[adata.var['pct_dropout_by_counts'] > 99.5].tolist()
    filter_genes.extend(highly_expressed_outiers)
    filter_genes.extend(lowly_expressed_outliers)

    # Subset
    adata = adata[:, ~adata.var_names.isin(filter_genes)].copy()

    # Get HVGs
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    rsc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, flavor='seurat_v3')
    hvg_df = adata.var[adata.var['highly_variable']][['gene_name', 'gene_ids']]
    return(hvg_df)

# --- Train model --- #

def train_contrastiveVI(
    adata, 
    background_indices,
    target_indices,
    batch_key='integration_technical_batch_concat',
    max_epochs = 500,
    use_observed_lib_size = False,
    wasserstein_penalty = 0
    ):
    scvi.external.ContrastiveVI.setup_anndata(adata, batch_key=batch_key)
    contrastive_vi_model = scvi.external.ContrastiveVI(
        adata, n_salient_latent=100, n_background_latent=20, use_observed_lib_size=use_observed_lib_size, wasserstein_penalty = wasserstein_penalty
    )

    contrastive_vi_model.train(
        background_indices=background_indices,
        target_indices=target_indices,
        batch_size=512,
        early_stopping=True,
        max_epochs=max_epochs,
    )
    
    adata.obsm["contrastiveVI_salient"] = contrastive_vi_model.get_latent_representation(representation_kind="salient")
    adata.obsm["contrastiveVI_background"] = contrastive_vi_model.get_latent_representation(representation_kind="background")
    return(contrastive_vi_model)

if __name__ == "__main__":
    import yaml
    parser = argparse.ArgumentParser(description='Embedding with contrastiveVI')
    parser.add_argument('--config', type=str,
                       default='/mnt/oak/users/emma/bin/GWT_perturbseq_analysis/metadata/experiments_config.yaml',
                       help='Path to experiment config YAML file')
    parser.add_argument('--experiment', type=str,
                       help='Experiment ID to process')
    parser.add_argument('--testing', action='store_true', help='test run with few epochs')
    args = parser.parse_args()
    
    if args.testing:
        n_epochs = 10
    else:
        n_epochs = 500
    
    config_file = args.config
    experiment_name = args.experiment

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    config = config[experiment_name]
    datadir = _convert_oak_path(config['datadir'])
    sample_metadata_csv = _convert_oak_path(config['sample_metadata'])
    sample_metadata = pd.read_csv(sample_metadata_csv, index_col=0)
    sgrna_library_metadata = pd.read_csv('../../metadata/sgRNA_library_curated.csv', index_col=0)

    adata = sc.read_h5ad(f'{datadir}/tmp/{experiment_name}_merged.gex.h5ad', backed=True)
    adata.obs = pd.read_csv(f'{datadir}/{experiment_name}_merged.gex.lognorm.postQC_obs.csv', compression='gzip', index_col=0)
    adata.var_names = adata.var['gene_ids'].values
    
    # Exclude cells with QC mask 
    adata = adata[adata.obs['QC_mask'] & (adata.obs['guide_id'] != 'multi_sgRNA')].to_memory()
    print('Subsetting done')

    # Feature selection
    hvgs = feature_selection(adata, subset_adata=True)['gene_ids']
    adata = adata[:, np.intersect1d(hvgs, adata.var['gene_ids'])].copy()
    print('Feature selection done')

    model = train_contrastiveVI(
        adata, 
        batch_key='donor_id',
        background_indices = np.where(adata.obs["guide_type"] == "non-targeting")[0],
        target_indices = np.where(adata.obs["guide_type"] == "targeting")[0],
        max_epochs = n_epochs
    )

    # Save outputs
    # Save highly variable genes
    hvgs_df = pd.DataFrame({'gene_ids': hvgs})
    hvgs_df.to_csv(f'{experiment_name}.contrastiveVI_hvgs.csv', index=False)
    model.save(f'{experiment_name}.model_contrastiveVI', overwrite=True, save_anndata=False)
    np.savez(f'{experiment_name}.X_contrastiveVI_salient.npz', array=adata.obsm["contrastiveVI_salient"])
    np.savez(f'{experiment_name}.X_contrastiveVI_background.npz', array=adata.obsm["contrastiveVI_background"])

    # Make embeddings and save
    rsc.pp.neighbors(adata, use_rep='contrastiveVI_background', n_neighbors=50)
    rsc.tl.umap(adata)

    perturb_adata = adata[adata.obs['guide_type'] == 'targeting'].copy()
    rsc.pp.neighbors(perturb_adata, use_rep='contrastiveVI_salient', n_neighbors=50)
    rsc.tl.umap(perturb_adata)

    mudata.MuData({'background':adata, 'salient':perturb_adata}).write_h5mu(f'{experiment_name}.model_contrastiveVI.h5mu')
    

