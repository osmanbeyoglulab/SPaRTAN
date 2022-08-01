import pandas as pd
from sklearn.preprocessing import normalize
from anndata import AnnData
import numpy as np


def load_dorthea():
    '''

    Returns
    -------
    TF-gene interaction matrix used in the tutorial.rst and Nucleic Acids Research publication.
    '''

    #return pd.read_csv("https://raw.githubusercontent.com/osmanbeyoglulab/SPaRTAN/main/data/GenevsTF_SPaRTAN_sample.csv", index_col=0)
    return pd.read_csv(pd.read_csv("https://sites.pitt.edu/~xim33/data/SPaRTAN_input/D_pbmc.csv", index_col=0))

def load_tf_gene(adata=None, database="DoRothEA", expression_cutoff=0.001, min_genes_per_tf=10, max_jac_sim=0.8):
    '''

    Parameters
    ----------
    adata : AnnData, optional
        AnnData object containing gene expression.  Used to filter TFs that are not expressed in mRNA
    database : {"DoRothEA", "hTFtarget"}
        Name of database to use.  Default is "DoRothEA"
    expression_cutoff : float, optional
        Threshold to use for filtering TFs not expressed in adata.  If `adata` is specified, TFs with mean mRNA expression
        in the adata object less than `expression_cutoff` will be filtered out.
    min_genes_per_tf : int, optional
        TFs with less than this many target genes are filtered out.
    max_jac_sim : float, optional
        Pairs or groups of TFs that all have jacard similarity higher than `max_jac_sim` are merged together.


    Returns
    -------
    TF-gene interaction matrix
    '''



def pbmc(ct):
    '''
    Retreive AnnData object containing pre-processed data for PBMC dataset.

    Parameters
    ----------
    ct : {'B', 'CD4mem', 'CD4nav', 'CD8', 'CD14+MONO', 'CD16+MONO', 'DC', 'NK'}
        Cell type for which data is obtained

    Returns
    -------
    adata : AnnData
        PBMC data for cell type `ct`
    '''

    dataset_D="D_pbmc"
    dataset_P="Ppbmc5kNextGEM_"+ct
    dataset_Y="Ypbmc5kNextGEM_"+ct

    D_url="https://sites.pitt.edu/~xim33/data/SPaRTAN_input/"+dataset_D+".csv"
    P_url="https://sites.pitt.edu/~xim33/data/SPaRTAN_input/"+dataset_P+".csv"
    Y_url="https://sites.pitt.edu/~xim33/data/SPaRTAN_input/"+dataset_Y+".csv"

    D_ori = pd.read_csv(D_url, index_col=0)
    P_ori = pd.read_csv(P_url, index_col=0)
    Y_ori = pd.read_csv(Y_url, index_col=0)

    D_mat = D_ori.values
    P_mat = P_ori.values
    Y_mat = Y_ori.values

    D_mat = normalize(D_mat, norm="l2", axis=0)
    Y_mat = normalize(Y_mat, norm="l2",axis=0)
    P_mat = normalize(P_mat, norm="l2", axis=1)

    adata=AnnData(X=Y_ori.T)


    adata.obsm["protein"]=pd.DataFrame(P_mat, index=adata.obs_names, columns=P_ori.columns)
    adata.obsm["protein_raw_counts"]=P_ori
    adata.varm["tf_gene"]=pd.DataFrame(D_mat, index=adata.var_names, columns=D_ori.columns)
    adata.layers["normalized"]=Y_mat.T

    adata.obs["training"]= np.random.rand(adata.n_obs) < 0.8

    return adata
