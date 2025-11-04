'''
Adapting R code from Ahlman-Eltze et al. 2024
https://github.com/const-ae/linear_perturbation_prediction-Paper/blob/main/benchmark/src/run_linear_pretrained_model.R
'''
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import solve
import scipy

from sklearn.decomposition import TruncatedSVD
from anndata import AnnData

def solve_y_axb(Y, G=None, P=None, ridge_lambda=0.01):
    """
    Solve matrix equations of the form Y = GXP with ridge regularization.
    
    Parameters:
    -----------
    Y : numpy.ndarray
        The target matrix
    G : numpy.ndarray or None
        Left multiplication matrix
    P : numpy.ndarray or None
        Right multiplication matrix
    ridge_lambda : float
        Ridge regularization parameter for G and P

    Returns:
    --------
    dict
        A dictionary containing the k-dimentional solution W and the center (mean) of Y
    """
    # Input validation
    if not isinstance(Y, (np.ndarray, sparse.spmatrix)):
        raise TypeError("Y must be a numpy array or sparse matrix")
    if G is not None and not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise TypeError("G must be a numpy array or sparse matrix")
    if P is not None and not isinstance(P, (np.ndarray, sparse.spmatrix)):
        raise TypeError("P must be a numpy array or sparse matrix")
    
    # Center Y by row means
    center = np.mean(Y, axis=1, keepdims=True)
    Y_centered = Y - center
    
    if G is not None and P is not None:
        # Check dimensions
        if Y.shape[0] != G.shape[0]:
            raise ValueError("Number of rows in Y must match number of rows in G")
        if Y.shape[1] != P.shape[1]:
            raise ValueError("Number of columns in Y must match number of columns in P")
            
        # Solve Y = GKP
        G_term = G.T @ G + ridge_lambda * np.eye(G.shape[1])
        P_term = P @ P.T + ridge_lambda * np.eye(P.shape[0])
        tmp = solve(G_term, G.T @ Y_centered @ P.T @ solve(P_term, np.eye(P_term.shape[0])))
        
    elif P is None and G is not None:
        # Solve Y = GK
        G_term = G.T @ G + ridge_lambda * np.eye(G.shape[1])
        tmp = solve(G_term, G.T @ Y_centered)
        
    elif G is None and P is not None:
        # Solve Y = KP
        P_term = P @ P.T + ridge_lambda * np.eye(P.shape[0])
        tmp = Y_centered @ P.T @ solve(P_term, np.eye(P_term.shape[0]))
        
    else:
        raise ValueError("Either G or P must be non-null")
    
    # Replace NaN values with 0
    tmp = np.nan_to_num(tmp, nan=0.0)
    
    return {"K": tmp, "center": center.flatten()}



class PerturbLinearEmbedding:
    """
    A linear model for perturbation prediction using PCA embeddings of genes.
    """
    
    def __init__(self, pca_dim = 30, ridge_lambda=0.1):
        """
        Initialize the PerturbLinearEmbedding model.
        
        Parameters:
        -----------
        ridge_lambda : float, default=1e-4
            Regularization parameter for ridge regression.
        """
        self.ridge_lambda = ridge_lambda
        self.pca_dim = pca_dim
        
        self.baseline = None
        self.perturbations = None
        self.genes = None
        self.K = None
        self.center = None
    
    def fit(self, train_adata, use_layer = 'mean', control = 'NTC'):
        """
        Fit the linear model to the data.
        
        Parameters:
        -----------
        train_adata : AnnData
            AnnData object containing the training data.
        use_layer : str, default='mean'
            The layer in the AnnData object to use for fitting the model.
            
        Returns:
        --------
        self : PerturbLinearEmbedding
            The fitted model.
        """
        # Store mean expression baseline
        self.baseline = train_adata[control].layers[use_layer].flatten()
        train_adata.layers['change'] = train_adata.layers[use_layer] - self.baseline
        keep_targets = train_adata.obs_names.tolist()
        keep_targets.remove(control)
        train_adata = train_adata[keep_targets].copy()

        Y = train_adata.layers['change'].T.copy()
        self.train_perturbations = train_adata.obs_names
        self.genes = train_adata.var_names

        # Compute PCA on the genes
        pca = TruncatedSVD(n_components=self.pca_dim)
        Y = Y.toarray() if hasattr(Y, 'toarray') else Y
        gene_emb = pca.fit_transform(Y)
        gene_emb_df = pd.DataFrame(gene_emb, index=self.genes)
        pert_emb_df = gene_emb_df.loc[self.train_perturbations].copy()

        # Fit LM 
        coefs = solve_y_axb(Y=Y, G=gene_emb_df.values, P=pert_emb_df.values.T, ridge_lambda=self.ridge_lambda)
        
        self.K = coefs['K']
        self.center = coefs['center']
        self.gene_embeddings = gene_emb_df.copy()
    
    def predict(self, perturbations = None, add_baseline = True):
        """
        Predict the effects of new perturbations using the fitted model.
        
        Parameters:
        -----------
        G_new : numpy.ndarray or scipy.sparse matrix, optional
            New gene embedding matrix.
        P_new : numpy.ndarray or scipy.sparse matrix, optional
            New perturbation embedding matrix.
            
        Returns:
        --------
        Y_pred : numpy.ndarray
            The predicted effects matrix.
        """
        if perturbations is None:
            perturbations = self.train_perturbations
        
        # Get pre-trained embedding for perturbations of interest
        gene_emb_df = self.gene_embeddings
        pert_emb_df = gene_emb_df.loc[perturbations].copy()

        # Predict effect on genes
        GWPt = np.matmul(np.matmul(gene_emb_df.values, self.K), pert_emb_df.values.T)
        prediction = GWPt + self.center[:, np.newaxis] 
        if add_baseline:
            prediction = prediction + self.baseline[:, np.newaxis]
        
        prediction_df = pd.DataFrame(prediction, index=self.genes, columns=perturbations)
        return(prediction_df)
        
def eval_prediction(groundtruth, predictions):
    '''Evaluate perturbation prediction using vectorized operations
    
    Parameters:
    -----------
    groundtruth : pd.DataFrame
        Ground truth values for each perturbation
    predictions : pd.DataFrame
        Predicted values for each perturbation
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing perturbation names, their L2 distances, and cosine distances
    '''
    # Ensure inputs are aligned
    if not groundtruth.columns.equals(predictions.columns):
        raise ValueError("Ground truth and predictions must have the same columns")
    
    # Compute L2 distances for all perturbations at once
    l2_distances = np.sqrt(np.sum((groundtruth - predictions) ** 2, axis=0))
    
    # Compute MSE for all perturbations at once
    rmse_values = np.mean((groundtruth - predictions) ** 2, axis=0)
    
    # Compute cosine distances and pearson correlations for all perturbations
    cosine_distances = []
    pearson_correlations = []
    pearson_pvalues = []
    for col in groundtruth.columns:
        gt_vec = groundtruth[col].values
        pred_vec = predictions[col].values
        # Compute cosine distance
        cosine_dist = 1 - np.dot(gt_vec, pred_vec) / (np.linalg.norm(gt_vec) * np.linalg.norm(pred_vec))
        cosine_distances.append(cosine_dist)
        # Compute pearson correlation
        pearson_corr, pearson_pval = scipy.stats.pearsonr(gt_vec, pred_vec)
        pearson_correlations.append(pearson_corr)
        pearson_pvalues.append(pearson_pval)
    
    # Create results DataFrame directly
    results_df = pd.DataFrame({
        'perturbation': groundtruth.columns,
        'l2_distance': l2_distances,
        'rmse': rmse_values,
        'cosine_distance': cosine_distances,
        'pearson_r': pearson_correlations,
        'pearson_pval': pearson_pvalues
    })
    
    return results_df
