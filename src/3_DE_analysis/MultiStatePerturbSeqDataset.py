import pandas as pd
import numpy as np
import anndata
from anndata import AnnData
import scanpy as sc
import pertpy
import logging
from tqdm import tqdm
from pathlib import Path
from typing import List, Dict, Optional, Union, Set
from pertpy.tools._differential_gene_expression._base import LinearModelBase
from pertpy.tools._differential_gene_expression._edger import EdgeR
from pertpy.tools._differential_gene_expression._pydeseq2 import PyDESeq2
import scipy.sparse


class MultistatePerturbSeqDataset:
    """Multi-state Perturb-seq dataset class."""

    def __init__(
        self,
        data: AnnData,
        sample_cols: list[str] = None,
        perturbation_type: str = None,
        target_col: str = "target",
        sgrna_col: str = "sgrna",
        state_col: str = "state",
        control_level: str = "NO-TARGET",
        **pseudobulk_kwargs,
    ):
        """Initialize the dataset object."""
        if isinstance(data, AnnData):
            if "MultistatePerturbSeqDataset" in data.uns:
                # When loading saved objects
                self.adata = data
                # Load all attributes in uns
                self.cell_state_obs = self.adata.uns["MultistatePerturbSeqDataset"]["state_col"]
                self.control_level = self.adata.uns["MultistatePerturbSeqDataset"]["control_level"]
                self.is_pseudobulked = self.adata.uns["MultistatePerturbSeqDataset"]["is_pseudobulked"]
                self.filtered_genes = self.adata.uns["MultistatePerturbSeqDataset"]["filtered_genes"]
                # self.model_fit_DE_targets = self.adata.uns["MultistatePerturbSeqDataset"]["model_fit_DE_targets"]

                del self.adata.uns["MultistatePerturbSeqDataset"]
            else:
                adata = data.copy()
                # Rename columns target and guide
                adata.obs.rename(columns={target_col: "target"}, inplace=True)
                adata.obs.rename(columns={sgrna_col: "sgrna"}, inplace=True)

                # Convert target labels to use hyphens instead of underscores
                adata.obs["target"] = adata.obs["target"].str.replace("_", "-")
                self.control_level = control_level.replace("_", "-")

                self.cell_state_obs = state_col
                if self.control_level not in adata.obs["target"].unique():
                    raise ValueError("Control level not found in target column.")
                
                # Store unique cell states
                self.cell_states = None
                if self.cell_state_obs is not None:
                    self.cell_states = sorted(adata.obs[self.cell_state_obs].unique().tolist())

                # Define sample ID
                if isinstance(sample_cols, str):
                    sample_cols = [sample_cols]
                if state_col not in sample_cols:
                    sample_cols.append(state_col)
                if 'sgrna' not in sample_cols:
                    sample_cols.append('sgrna')
                adata.obs["sample_id"] = adata.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)
                adata.obs["perturbation_type"] = perturbation_type
                self.design_params = sample_cols + ['target']

                # Track pseudobulk state
                self.is_pseudobulked = False

                self.adata = adata
                # Build design matrix
                self.design_matrix = self._build_design_matrix()

                self.filtered_genes = False
                 # Dictionary to store fitted models by design formula
                self.fitted_models = {}


        else:
            raise ValueError("data must be an AnnData object.")

    # def __init__(
    #     self,
    #     adata: anndata.AnnData,
    #     dataset_config: dict,
    #     cell_state_obs: str = 'author_cell_state'
    # ):
    #     self.cell_state_obs = cell_state_obs
    #     self.config = dataset_config
    #     self.design_params = self.config['design_params']
        
    #     # Load dataset
    #     self.adata = adata
        
        
        
    #     # Store unique sample IDs
    #     self.sample_ids = self.design_matrix['sample_id'].tolist()

    #     # Store unique cell types
    #     self.cell_states = None
    #     if self.cell_state_obs is not None:
    #         self.cell_states = sorted(self.adata.obs[self.cell_state_obs].unique().tolist())
        
    #     # Dictionary to store fitted models by design formula
    #     self.fitted_models = {}

    #     # Track pseudobulk state
    #     self.is_pseudobulked = False

    def pseudobulk(
        self,
        sample_col: str = "sample_id",
        min_cells: int = 3
    ) -> None:
        """
        Create a pseudobulk version of the dataset by aggregating cells
        by sample ID and cell type.
        
        Args:
            sample_col: Column name for sample identifiers
            **pseudobulk_kwargs: Arguments passed to decoupler.get_pseudobulk
            
        Note:
            This modifies the AnnData object in place and sets is_pseudobulked flag
        """
        if self.is_pseudobulked:
            logging.warning("Data is already pseudobulked. Skipping.")
            return
            
        logging.info("Performing pseudobulk aggregation...")
        
        try:    
            adata_bulk = sc.get.aggregate(self.adata, by=sample_col, func=['sum', 'count_nonzero'])
            adata_bulk.obs = pd.concat([adata_bulk.obs, self.design_matrix], axis=1)
            adata_bulk.X = adata_bulk.layers['sum']
            del adata_bulk.layers['sum']
            
            # Convert to sparse matrix to save memory
            adata_bulk.X = scipy.sparse.csr_matrix(adata_bulk.X)

            total_cells_sample = self.calculate_cell_state_fractions()['total_cells']
            keep_samples = total_cells_sample[total_cells_sample >= min_cells].index.tolist()
            adata_bulk = adata_bulk[adata_bulk.obs_names.isin(keep_samples)].copy()
            self.design_matrix = self.design_matrix.loc[adata_bulk.obs_names].copy()
            
            # Update the object
            self.adata = adata_bulk
            self.is_pseudobulked = True
            
            # Clear any cached models since data has changed
            self.fitted_models = {}
            
            logging.info(f"Pseudobulk complete. Shape: {self.adata.shape}")
            
        except Exception as e:
            logging.error(f"Error during pseudobulk: {str(e)}")
            raise
    
    def _get_condition_mask(self, conditions: dict) -> pd.Series:
        """
        Create boolean mask for cells matching all conditions.
        
        Args:
            conditions: Dict of column:value pairs, where value can be single value or list
            
        Returns:
            Boolean mask for matching cells
        """
        masks = []
        for col, values in conditions.items():
            if isinstance(values, (list, tuple)):
                masks.append(self.adata.obs[col].isin(values))
            else:
                masks.append(self.adata.obs[col] == values)
        
        return pd.concat(masks, axis=1).all(axis=1)
    
    def run_target_DE(
        self, 
        design_formula: str = '~ donor_id + target', 
        test_targets: List[str] = None,
        test_state=None, 
        min_counts_per_gene = 3,
        return_model = False,
        n_cpus = 10
        ):
        '''
        Run differential expression analysis for each target compared to control.
        
        Args:
            design_formula: Formula string for the design matrix (default: '~ donor_id + target')
            test_targets: List of target names to test. If None, all targets will be tested.
            test_state: Cell state(s) to test. Can be a string or list of strings. 
                        If None, all states will be tested.
            min_counts_per_gene: Minimum number of counts required for a gene to be included in analysis.
            return_model: If True, returns the fitted model along with results.
            
        Returns:
            pandas.DataFrame with differential expression results for each tested target, or
            tuple of (model, pandas.DataFrame) if return_model is True.
        '''
        adata = self.adata

        all_res_df = pd.DataFrame()
        # Convert test_state to a list if it's a string
        if test_state is None:
            test_state = adata.obs[self.cell_state_obs].unique().tolist()
        if isinstance(test_state, str):
            test_state = [test_state]
            
        for st in test_state:
            if test_targets is None:
                adata_state = adata[adata.obs[self.cell_state_obs] == st]
            else:
                adata_state = adata[(adata.obs[self.cell_state_obs] == st) & (adata.obs['target'].isin(test_targets + [self.control_level]))]
                
            adata_state = adata_state[:, adata_state.X.sum(0) >= min_counts_per_gene].copy()

            model = pertpy.tl.PyDESeq2(adata_state, design=design_formula)
            model.fit(n_cpus = n_cpus, quiet=True)
    
            all_targets = adata_state.obs['target'].unique().tolist()
            all_targets.remove(self.control_level)

            for t in tqdm(all_targets, desc="Testing targets"):
                t_contrast = (model.cond(target = t) - model.cond(target = self.control_level)) 
                res_df = model.test_contrasts(t_contrast)
                res_df[self.cell_state_obs] = st
                res_df['contrast'] = t
                all_res_df = pd.concat([all_res_df, res_df])
        
        all_res_df = all_res_df.reset_index().drop('index', axis=1)
        
        if return_model:
            return model, all_res_df    
        else:
            return all_res_df

    # def run_target_DE_comparison(
    #     self,
    #     comparison: Dict,
    #     force_fit_model: bool = False
    #     ) -> pd.DataFrame:
    #     """
    #     Run differential expression analysis for a given comparison.
        
    #     Args:
    #         comparison: Dictionary containing comparison configuration with keys:
    #             - name: Name of the comparison
    #             - experimental: Dict of experimental condition values
    #             - control: Dict of control condition values 
    #             - covariates: Optional list of covariates
                
    #     Returns:
    #         DataFrame containing differential expression results
    #     """
    #     test_covs = list(set(list(comparison['experimental'].keys()) + list(comparison['control'].keys())))
    #     covariates = comparison.get('covariates', []) or []

    #     temp_col = '_temp_comparison'
    #     self.adata.obs[temp_col] = self.adata.obs[test_covs].astype(str).agg('_'.join, axis=1)
    #     experimental_level = "_".join([comparison['experimental'][c] for c in test_covs])
    #     control_level = "_".join([comparison['control'][c] for c in test_covs])

    #     # Check other models with the same covariates
    #     existing_model = None
    #     if not force_fit_model:
    #         for model_name, model_config in self.fitted_models.items():
    #             if (set(model_config['test_covs']) == set(test_covs) and 
    #                 set(model_config['covariates']) == set(covariates)):
    #                 existing_model = model_config['model']
    #                 break

    #     # Fit new model if needed
    #     if existing_model is not None:
    #         model = existing_model
    #     else:
    #         design_formula = f"~{' + '.join(covariates)} + {self.cell_state_obs} + {self.cell_state_obs}*{temp_col}"
    #         model = pertpy.tl.PyDESeq2(self.adata, design=design_formula)
    #         model.fit(quiet=True)
            
    #         self.fitted_models[comparison['name']] = {
    #             'model': model,
    #             'test_covs': test_covs,
    #             'covariates': covariates
    #         }

    #     # Test contrasts
    #     all_res_df = pd.DataFrame()
    #     for ct in self.cell_states:
    #         contrast = model.cond(**{self.cell_state_obs:ct, temp_col:experimental_level}) - model.cond(**{self.cell_state_obs:ct, temp_col:control_level})
    #         res_df = model.test_contrasts(contrast)
    #         res_df = res_df.set_index("variable")
    #         res_df = pd.concat([res_df, self.adata.var], axis=1)
    #         res_df['cell_state'] = ct
    #         res_df['comparison'] = comparison['name']
    #         res_df['contrast'] = ct
    #         all_res_df = pd.concat([all_res_df, res_df])
        
    #     return all_res_df
    
    def _build_design_matrix(self) -> pd.DataFrame:
        """
        Build design matrix from the dataset.
        
        Returns:
            DataFrame with rows for each unique sample and columns for design parameters
        """
        # Get unique samples and their metadata
        sample_meta = (self.adata.obs
                      .groupby('sample_id')
                      .first()
                      .reset_index())
        
        # Select relevant columns
        meta_cols = ['sample_id'] + [
            col for col in self.design_params 
            if col in sample_meta.columns
        ]
        design_matrix = sample_meta[meta_cols]
        
        # Check for missing required columns
        missing_cols = set(self.design_params) - set(design_matrix.columns)
        if missing_cols:
            print(f"Warning: Missing design parameters: {missing_cols}")
            
        return design_matrix.set_index('sample_id')
        
    def calculate_cell_state_fractions(self) -> pd.DataFrame:
        """
        Calculate the fraction of each cell type within each sample.
        
        Returns:
            DataFrame with rows for each sample and columns for:
            - Cell type fractions
            - Total cell count
            - Design parameters from self.design_params
        """
        # Validate required columns
        if self.cell_state_obs not in self.adata.obs.columns:
            raise ValueError(f"Missing cell type column: {self.cell_state_obs}")
        
        # Get total cells per sample
        sample_totals = self.adata.obs.groupby('sample_id').size()
        
        # Count cells of each type per sample
        type_counts = (self.adata.obs
                      .groupby(['sample_id', self.cell_state_obs])
                      .size()
                      .unstack(fill_value=0))
        
        # Calculate fractions
        fractions = type_counts.div(sample_totals, axis=0)
        
        # Add total cell counts
        fractions['total_cells'] = sample_totals
        
        # Add design parameters
        for col in self.design_params:
            if col in self.design_matrix.columns:
                fractions[col] = self.design_matrix[col]
        
        fractions = fractions.reset_index()
        fractions.columns.name = None
        
        return fractions.set_index('sample_id')
    
    def save(self, h5ad_file: str):
        """
        Save the dataset to an h5ad file and a companion pickle file for the models.

        Parameters
        ----------
        h5ad_file : str
            Path to save the h5ad file. The models will be saved to a companion
            pickle file with '_models.pkl' appended to the base name.
        """
        # Convert h5ad_file to Path object for easier manipulation
        h5ad_path = Path(h5ad_file)

        # Save all other attributes in uns
        self.adata.uns["MultistatePerturbSeqDataset"] = {
            "state_col": self.cell_state_obs,
            "control_level": self.control_level,
            "is_pseudobulked": self.is_pseudobulked,
            "filtered_genes": self.filtered_genes,
            # "model_fit_DE_targets": self.model_fit_DE_targets, #TODO: save trained models
        }

        # Save the AnnData object
        self.adata.write(h5ad_file)

        # Clean up uns
        del self.adata.uns["MultistatePerturbSeqDataset"]