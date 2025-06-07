import numpy as np
import pandas as pd
from typing import Dict, Union
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import RepeatedKFold
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


from scipy.stats import spearmanr, pearsonr, kendalltau
from sklearn.metrics import r2_score
import scipy


import matplotlib.pyplot as plt
import seaborn as sns

def score_predictions(y: np.ndarray, y_pred: np.ndarray, decimals: int = 2) -> Dict:
    """
    Calculate prediction metrics including correlations and rank agreements.
    """
    # Calculate correlations with fallback for constant arrays
    constant = np.all(y == y[0])
    constant_pred = np.all(y_pred == y_pred[0])
    if not constant and not constant_pred:
        spearman_corr, spearman_p = spearmanr(y, y_pred)
        pearson_corr, pearson_p = pearsonr(y, y_pred) 
        kendall_corr, kendall_p = kendalltau(y, y_pred)
    else:
        spearman_corr, spearman_p = 0, 1
        pearson_corr, pearson_p = 0, 1
        kendall_corr, kendall_p = 0, 1
    
    metrics = {
        'r2': round(r2_score(y, y_pred), decimals),
        'spearman': round(spearman_corr, decimals),
        'pearson': round(pearson_corr, decimals),
        'kendall': round(kendall_corr, decimals),
        'spearman_p': round(spearman_p, decimals),
        'pearson_p': round(pearson_p, decimals),
        'kendall_p': round(kendall_p, decimals)
    }
   
    return metrics

class Perturb2TissueModel:
    def __init__(
        self, 
        n_splits: int = 5,
        n_repeats: int = 10,
        random_state: int = 42,
        pca_transform: bool = False,
        n_pcs: int = 50,
        positive: bool = False
    ):
        """
        Initialize the model with cross-validation parameters.
        
        Args:
            n_splits: Number of folds in k-fold cross-validation (controls size of train/test split)
            n_repeats: Number of times to repeat the cross-validation
            random_state: Random seed for reproducibility
        """
        self.n_splits = n_splits
        self.n_repeats = n_repeats
        self.random_state = random_state
        self.pca_transform = pca_transform
        self.n_pcs = n_pcs
        self.positive = positive
    
    def evaluate_single_split(
        self, 
        X: pd.DataFrame, 
        y: pd.Series,
        train_idx: np.ndarray,
        test_idx: np.ndarray,
        return_model: bool = False
    ) -> Dict:
        """Evaluate a single train-test split."""
        # Split data
        X_train, X_test = X.iloc[train_idx].values, X.iloc[test_idx].values
        y_train, y_test = y.iloc[train_idx].values, y.iloc[test_idx].values

        # Scale the data
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        
        if self.pca_transform:
            # Reduce dimensions with PCA
            pca = TruncatedSVD(n_components=self.n_pcs)
            pca.fit(X_train)
            pca_X_train = pca.transform(X_train)
            pca_X_test = pca.transform(X_test)
            X_train_input = pca_X_train.copy()
            X_test_input = pca_X_test.copy()
        else:
            X_train_input = X_train.copy()
            X_test_input = X_test.copy()
            pca = None

        # Fit model with internal cross-validation for hyperparameter selection
        model = ElasticNetCV(
            l1_ratio=np.logspace(-2, 2, 20),
            alphas=np.logspace(-3, 3, 30),
            tol=1e-3,
            max_iter=3000,
            cv=4,
            random_state=self.random_state,
            positive=self.positive
        )

        # Fit and predict
        model.fit(X_train_input, y_train)
        y_pred_train = model.predict(X_train_input)
        y_pred_test = model.predict(X_test_input)
        
        # Convert numpy arrays to Python floats before rounding
        train_score = score_predictions(y_train.astype(float), y_pred_train.astype(float))
        test_score = score_predictions(y_test.astype(float), y_pred_test.astype(float))

        results = {
            'alpha': float(model.alpha_),
            'l1_ratio': float(model.l1_ratio_),
        }
        for key, value in train_score.items():
            results[f'train_{key}'] = value
        for key, value in test_score.items():
            results[f'test_{key}'] = value
        if return_model:
            results['model'] = model
            results['pca'] = pca
            results['scaler'] = scaler
        return results
    
    def fit(
        self, 
        X: pd.DataFrame, 
        y: pd.Series,
        model_id: str,
        frac_top_hvgs = None,
    ) -> pd.DataFrame:
        """
        Model latent transcriptional effect (DE estimates) from physiological or unknown condition as the weighted sum of effects
        of known perturbations, across multiple splits. Uses elastic net regression with cross-validation to find optimal 
        hyperparameters and evaluate model performance.

        Args:
            X: DataFrame of differential expression estimates from known perturbations.
               Shape (n_genes, n_perturbations) where columns are perturbation names and rows are gene IDs.
            y: Series of differential expression estimates from the physiological/unknown condition.
               Shape (n_genes,) with gene IDs as index matching X.
            model_id: String identifier for the model, used to track results across different conditions/experiments.
            frac_top_hvgs: Optional float between 0 and 1. If provided, uses only this fraction of genes with highest
                          mean absolute expression for training. If None, uses all genes.

        Returns:
            DataFrame containing evaluation metrics for each fold:
                - alpha: Optimal L1 ratio hyperparameter
                - l1_ratio: Optimal elastic net mixing parameter
                - train_*/test_*: Various performance metrics on training and test sets
                - fold: Index of cross-validation fold
                - model_id: Input model identifier
        """
        # Initialize cross-validation
        rkf = RepeatedKFold(
            n_splits=self.n_splits,
            n_repeats=self.n_repeats,
            random_state=self.random_state
        )
        
        # Input validation
        if not isinstance(X, pd.DataFrame):
            raise TypeError("X must be a pandas DataFrame")
        if not isinstance(y, pd.Series):
            raise TypeError("y must be a pandas Series")
        if not isinstance(model_id, str):
            raise TypeError("model_id must be a string")
        if frac_top_hvgs is not None:
            if not isinstance(frac_top_hvgs, (int, float)):
                raise TypeError("frac_top_hvgs must be a number between 0 and 1")
            if frac_top_hvgs <= 0 or frac_top_hvgs > 1:
                raise ValueError("frac_top_hvgs must be between 0 and 1")
        
        # Check that X and y have compatible indices
        if not all(idx in X.index for idx in y.index):
            raise ValueError("y contains indices not present in X")
        results = []
        self.split_ixs = []

        # Select highly variable features in X
        if frac_top_hvgs is not None:
            variances = pd.Series(X.var(axis=1))
            means = pd.Series(np.abs(X).mean(axis=1))
            # top_features = variances.nlargest(int(len(variances) * frac_top_hvgs)).index.tolist()
            top_features = means.nlargest(int(len(means) * frac_top_hvgs)).index.tolist()
            print(f'Training on {len(top_features)}/{X.shape[0]} genes')
        else:
            top_features = X.index.tolist()

        X_filtered = X.loc[top_features]
        y_filtered = y.loc[top_features]

        for fold_idx, (train_idx, test_idx) in enumerate(rkf.split(X_filtered)):
            result = self.evaluate_single_split(X_filtered, y_filtered, train_idx, test_idx, return_model=True)
            result.update({
                'fold': fold_idx,
                'model_id': model_id,
                'train_idx': train_idx,  # Store train indices
                'test_idx': test_idx     # Store test indices
            })
            results.append(result)
            self.split_ixs.append((train_idx, test_idx))  # Store indices in split_ixs
        
        self.eval = pd.DataFrame(results).drop(['model', 'pca', 'scaler'], axis=1)
        self.models = [x['model'] for x in results]
        self.pcas = [x['pca'] for x in results]
        self.scalers = [x['scaler'] for x in results]
        self.perturbation_names = X.columns.tolist()
        self.model_features = top_features
    
    def get_coefs(self) -> pd.DataFrame:
        '''Get mean and SEM of model coefficients for each known perturbation.'''

        mean_coefs_df = pd.DataFrame()
        for mod, pca in zip(self.models, self.pcas):
            if pca is not None:
                coefs = np.matmul(
                    pca.components_.T, 
                    mod.coef_
                )
            else:
                coefs = mod.coef_
            coefs_df = pd.DataFrame(coefs, index=self.perturbation_names)
            mean_coefs_df = pd.concat([mean_coefs_df, coefs_df], axis=1)

        # Calculate mean and sem across methods and store in mean_coefs_df
        mean_coefs_df['coef_mean'] = mean_coefs_df.mean(axis=1)
        mean_coefs_df['coef_sem'] = mean_coefs_df.sem(axis=1)
        return(mean_coefs_df[['coef_mean', 'coef_sem']])
    
    def _mod_residuals(self, mod, pca, scaler, X, y):
        X_scaled = scaler.transform(X)
        X_input = pca.transform(X_scaled) if pca is not None else X_scaled
        y_pred = mod.predict(X_input)
        residuals = y - y_pred
        return(residuals)
    
    def get_residuals(self, X, y):
        all_res = pd.DataFrame()
        for mod, pca, scaler in zip(self.models, self.pcas, self.scalers):
            # Convert X to numpy array to avoid feature names warning
            X_array = X.values if hasattr(X, 'values') else X
            res = self._mod_residuals(mod, pca, scaler, X_array, y)
            all_res = pd.concat([all_res, res], axis=1)
        # Calculate mean and sem across models
        all_res['residual_mean'] = all_res.mean(axis=1)
        all_res['residual_sem'] = all_res.sem(axis=1)
        return(all_res[['residual_mean', 'residual_sem']])
    
    def get_prediction(self, X, return_splits=False):
        """Get mean and SEM of predictions across models.
        
        Args:
            X: Input features
            return_splits: If True, return predictions for each split instead of mean/SEM
        """
        all_preds = pd.DataFrame()
        for mod, pca, scaler in zip(self.models, self.pcas, self.scalers):
            # Convert X to numpy array to avoid feature names warning
            X_array = X.values if hasattr(X, 'values') else X
            # Scale the data
            X_scaled = scaler.transform(X_array)
            # Apply PCA transformation if needed
            X_input = pca.transform(X_scaled) if pca is not None else X_array
            preds = pd.DataFrame(mod.predict(X_input), index=X.index)
            all_preds = pd.concat([all_preds, preds], axis=1)
        
        if return_splits:
            return all_preds
        
        # Calculate mean and sem across models
        all_preds['pred_mean'] = all_preds.mean(axis=1)
        all_preds['pred_sem'] = all_preds.sem(axis=1)
        return all_preds[['pred_mean', 'pred_sem']]

    def get_r2(self, X, y):
        """Calculate R² values for each feature's contribution to the prediction.
        
        Args:
            X: Input features
            y: Target values
            
        Returns:
            DataFrame with R² values for each feature, sorted in descending order
        """
        # Get coefficients from the model
        coefs = self.get_coefs()
        
        # Initialize empty dataframe to store r2 values
        r2_df = pd.DataFrame(index=X.columns, columns=['r2'])
        
        # Calculate r2 for each feature
        for col in X.columns:
            # Get feature values
            X_g = X[col]
            # Mask self-perturbation and fill with mean
            X_g = X_g.mask(X_g.index == col).fillna(X_g.mean())
            # Calculate predictions using coefficients
            pred_from_g = X_g * coefs['coef_mean'].loc[col]
            # Calculate correlation and r2
            corr, _ = scipy.stats.pearsonr(pred_from_g, y)
            r2_df.loc[col, 'r2'] = corr**2
            
        # Sort by r2 values and return
        return r2_df.sort_values('r2', ascending=False)

    def plot_residuals(self, X, y, annotate_top_n=0):
        """Plot residuals vs predicted values using mean residuals across models.
        
        Args:
            X: Input features
            y: Target values
            annotate_top_n: Number of top/bottom points to annotate with their index names
        """
        residuals = self.get_residuals(X, y)
        
        # Get mean predictions across models
        y_pred = self.get_prediction(X)['pred_mean']
        
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.errorbar(y_pred, residuals['residual_mean'], 
                   yerr=residuals['residual_sem'],
                   fmt='o', markersize=4, alpha=0.5, 
                   elinewidth=1, capsize=0, color='black')
        ax.axhline(y=0, color='darkgrey', linestyle='--')
        ax.set_xlabel('Predicted DE estimate', fontsize=12)
        ax.set_ylabel('Residuals', fontsize=12)
        ax.tick_params(axis='both', labelsize=11)
        
        if annotate_top_n > 0:
            # Get indices of top and bottom n residuals
            sorted_idx = residuals['residual_mean'].abs().sort_values(ascending=False)
            top_idx = sorted_idx.index[:annotate_top_n]
            
            # Add annotations
            for idx in top_idx:
                ax.annotate(str(idx), 
                    (y_pred[idx],
                     residuals.loc[idx, 'residual_mean']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=11)
        
        return fig
    
    def plot_prediction(
            self,
            X: pd.DataFrame,
            y: pd.Series,
            annotate_top_n: int = 0,
            plot_metric: str = None,
            ax: plt.Axes = None,
            return_ax: bool = False) -> Union[plt.Figure, plt.Axes]:
        """Plot actual vs predicted values using mean predictions across models.
        
        Args:
            X: Input features as a pandas DataFrame
            y: Target values as a pandas Series
            annotate_top_n: Number of top/bottom points to annotate with their index names
            plot_metric: Optional metric name to display in plot title (e.g. 'r2', 'mse')
            ax: Optional matplotlib Axes object to plot on
            return_ax: If True, return the Axes object instead of Figure
            
        Returns:
            matplotlib Figure object containing the plot or Axes object if return_ax=True
        """
        # Get mean predictions and standard errors across models
        predictions = self.get_prediction(X)
        y_pred = predictions['pred_mean']
        y_pred_sem = predictions['pred_sem']
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        
        ax.errorbar(y, y_pred, yerr=y_pred_sem,
                   fmt='o', markersize=2, alpha=0.5,
                   elinewidth=1, capsize=0, color='black')
        
        # # Add diagonal line representing perfect prediction
        # lims = [
        #     min(min(y), min(y_pred)),
        #     max(max(y), max(y_pred))
        # ]
        # ax.plot(lims, lims, '--', color='darkgrey')
        
        ax.set_xlabel('Actual DE estimate', fontsize=12)
        ax.set_ylabel('Predicted DE estimate', fontsize=12)
        ax.tick_params(axis='both', labelsize=11)
        
        if annotate_top_n > 0:
            # Get indices of top and bottom n values
            sorted_idx = y_pred.abs().sort_values(ascending=False).index
            top_idx = sorted_idx[:annotate_top_n]
            
            # Add annotations
            for idx in top_idx:
                ax.annotate(idx,
                    (y[idx], y_pred[idx]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=11)
        
        if plot_metric is not None:
            # Get evaluation metrics
            eval_summary = self.summarize_eval()
            metric_col = f'test_{plot_metric}'
            if metric_col in eval_summary.columns:
                metric_mean = eval_summary[metric_col].mean()
                metric_sem = eval_summary[f'{metric_col}_se'].mean()
                ax.set_title(f'Test {plot_metric}: {metric_mean:.3f} ± {metric_sem:.3f}',
                           fontsize=12)
        
        if return_ax:
            return ax
        else:
            if 'fig' in locals():
                return fig
            else:
                return ax.figure
    
    def plot_coefs(self, figsize=(10,6), show_labels=True, top_n=None, bottom_n=None, ax=None, return_ax=False) -> Union[plt.Figure, plt.Axes]:
        '''Plot mean coefficients with error bars, sorted by coefficient value.
        
        Args:
            figsize: Figure size as (width, height) tuple
            show_labels: Whether to show feature labels on y-axis
            top_n: If provided, show only the top n coefficients by absolute value
            bottom_n: If provided, show only the bottom n coefficients by absolute value
            ax: Optional matplotlib axes to plot on
            return_ax: If True, return the axes object instead of the figure
        
        Returns:
            Either the figure or axes object depending on return_ax parameter
        '''
        
        # Get coefficients and sort by mean value
        coefs = self.get_coefs()
        coefs = coefs.sort_values('coef_mean', ascending=True)
        
        # Filter coefficients if top_n or bottom_n is specified
        if top_n is not None or bottom_n is not None:
            # Sort by absolute value for top_n selection
            coefs_by_abs = coefs.sort_values('coef_mean', key=abs, ascending=False)
            
            selected_indices = []
            if top_n is not None:
                selected_indices.extend(coefs_by_abs.index[:top_n])
                
            # Filter and re-sort by actual value
            coefs = coefs.loc[selected_indices].sort_values('coef_mean', ascending=True)
        
        # Create plot or use provided axes
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        
        y_pos = np.arange(len(coefs))
        
        # Plot points and error bars using errorbar instead of separate scatter and hlines
        ax.errorbar(coefs['coef_mean'], y_pos, 
                   xerr=coefs['coef_sem'],
                   fmt='o',  # Use circle markers
                   color='black',
                   capsize=3,  # Add small caps to error bars
                   capthick=1,  # Cap thickness
                   elinewidth=1)  # Error bar line width
        # Configure y-axis
        ax.set_yticks(y_pos)
        if show_labels:
            ax.set_yticklabels(coefs.index, fontsize=11)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel('Model coefficient', fontsize=12)
        ax.tick_params(axis='x', labelsize=11)
        
        # Add vertical line at x=0
        ax.axvline(x=0, color='darkgrey', linestyle='--', linewidth=0.5)
        
        if ax is None:
            plt.tight_layout()
        
        if return_ax:
            return ax
        else:
            if ax is None:
                return fig
            else:
                return ax.figure

    def summarize_eval(self) -> pd.DataFrame:
        """Summarize evaluation results across folds with means and standard errors."""
        def sem(x):
            return x.std() / np.sqrt(len(x))
        
        numeric_cols = ['alpha', 'l1_ratio'] + [col for col in self.eval.drop(['train_idx', 'test_idx'], axis=1).columns 
                                            if col.startswith(('train_', 'test_'))]
        
        # Calculate means
        means = (self.eval
            .groupby(['model_id'])[numeric_cols]
            .mean()
            .round(4))
        
        # Calculate SEMs
        sems = (self.eval
            .groupby(['model_id'])[numeric_cols]
            .agg(sem)
            .round(4))
        
        # Rename SEM columns
        sems.columns = [f"{col}_se" for col in sems.columns]
        
        # Combine means and SEMs
        summary = pd.concat([means, sems], axis=1)
        
        # Sort columns to keep related metrics together
        cols = []
        for col in means.columns:
            cols.extend([col, f"{col}_se"])
        
        return summary[cols]