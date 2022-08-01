import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
def normalize(A, axis=0):
    '''

    Parameters
    ----------
    A
    axis

    Returns
    -------

    '''
    return normalize_column(A, axis=axis)

def normalize_column(A, axis=0):

    """ perform l2 normalization column-wize of given matrix

    Parameters:
        A : the matrix that works on
        axis : switch of column-wize and row-wize.
            axis=0: column-wize
            axis=1: row-wize
    """

    if (axis == 0):
        return np.divide(A, np.sqrt(np.sum(A**2, 0)))
    else:
        At = np.transpose(A)
        return np.transpose(np.divide(At, np.sqrt(np.sum(At**2, 0))))

def clr(X):
    '''
    Center log ratio transforms input matrix
    Parameters
    ----------
    X : DataFrame or array_like
        input matrix

    Returns
    -------
    x_norm : DataFrame or array_like
        normalized matrix
    '''
    if type(X) is pd.DataFrame:
        x=X.to_numpy()
    else:
        x=X.copy()
    x_norm=np.log1p(x/np.exp(np.mean(np.log1p(x), axis=1)).reshape((-1,1)))
    if type(X) is pd.core.frame.DataFrame:
        x_norm=pd.DataFrame(x_norm, columns=X.columns, index=X.index)
    return x_norm


def subsample_celltype(adata, obs_name="cell_types", n_cells=700, random_state=0):
    '''
    Subsets data set to contain no more than a given amount of any one cell type.

    Parameters
    ----------
    adata : AnnData
        input dataset
    obs_name : str, optional
        field in `adata.obs` corresponding to cell types.  Default: "cell_types"
    n_cells : int, optional
        maximum number of cells for any cell type.
    random_state : int, optional
        seed used for reproducibility.

    Returns
    -------
    adata : AnnData
        subsampled dataset
    '''

    adatas=dict()
    for ct in np.unique(adata.obs[obs_name]):
        adatas[ct]=adata[(adata.obs[obs_name] == ct)]
        if adatas[ct].n_obs>n_cells:
            adatas[ct]=sc.pp.subsample(adatas[ct], n_obs=n_cells, random_state=random_state, copy=True)
    return ad.concat(adatas, label="cell_types")

