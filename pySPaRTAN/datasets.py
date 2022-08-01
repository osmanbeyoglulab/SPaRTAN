import pandas as pd
from sklearn.preprocessing import normalize
from anndata import AnnData
import numpy as np


def load_dorthea():
    #return pd.read_csv("https://raw.githubusercontent.com/osmanbeyoglulab/SPaRTAN/main/data/GenevsTF_SPaRTAN_sample.csv", index_col=0)
    return pd.read_csv(pd.read_csv("https://sites.pitt.edu/~xim33/data/SPaRTAN_input/D_pbmc.csv", index_col=0))
def pbmc(ct):
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
