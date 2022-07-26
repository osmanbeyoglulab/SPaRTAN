
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Xiaojun Ma
# Created Date: Mar 18 10:54:00 PDT 2020
# =============================================================================
"""
This script contains the major class of SPaRTAN model and its dependencies.

This script requires numpy, Scipy, matplotlib to be installed within the Python
environment you are running this script in

This script requires Cython modules present in the current directory

This file contains the following classes and functions

    class Timer: a class to convert time period in seconds to the format of h:m:s

    class pySPaRTAN: The major class for SPaRTAN, establishing an interaction matrix between
    surface proteins (P) and TFs (D) that predict target gene expression (Y).

    function normalize_column(): perform l2 normalization column-wize of given matrix

"""
import numpy as np
import pySPaRTAN.cythKronPlus as krnP
import pySPaRTAN.cythLeastR as leastR
import scipy.linalg, scipy.sparse
import gc
from copy import deepcopy
import statsmodels.stats.multitest
from scipy import stats
import warnings
from tqdm.auto import tqdm
import pandas as pd
from anndata import AnnData
from sklearn.preprocessing import normalize
import numpy as np
import anndata

class SPaRTAN:
    """
    The major class for SPaRTAN, establishing an interaction matrix between
    surface proteins (P) and TFs (D) that predicts target gene expression (Y).

    Methods
    -------
    fit(self, D, P, Y, lamda=0.001, rsL2=0.001,
        spectrumP=0.7):
        train a SPaRTAN model

    ar_model2w(self):
        converts a trained model to intermidiat vaiable W

    ar_reconstruction(self, pred_test=None):
        reconstruction function

    predict(self, P_test=None):
        predict target gene expression

    get_corr(self, Y_pred, Y_test, plot=False):
        get the correlation between predicted Y_pred and Y_test

    get_W(self):
        get coefficient matrix

    get_projP(self, Y=None):
        get projected protein expression

    get_projD(self, P=None):
        get projected TF activity
    """
    def __init__(self, lamda=0.01, alpha=0.5, corrtype='pearson',grid_search_steps=10, n_folds=5, verbose=False, solver="SVD", spectrum=0.7):
        ''''''
        self.protein_key=None
        self.gene_key=None

        self.solver=solver
        self.lamda=lamda
        self.alpha=alpha
        self.corrtype=corrtype
        self.n_folds=n_folds
        self.spectrum=spectrum

        if corrtype=='spearman':
            self.corr=scipy.stats.spearmanr
        elif corrtype=='pearson':
            self.corr=scipy.stats.pearsonr
        else:
            warnings.warn("Unrecognized corrtype, using Spearman correlation")
            self.corr=scipy.stats.spearmanr

        self.grid_search_steps=grid_search_steps
        self.verbose=verbose

        self.obs_names=None
        self.gene_names=None
        self.tf_names=None
        self.protein_names=None

    def fit(self, adata=None, D="tf_gene", P="protein", Y="log1p"):
        '''

        IMPORTANT: the protein and gene expression matrices must be normalized prior to calling this function.  See ... for recomendations

        Parameters
        ----------
        adata : AnnData, if specified, the object contains the TF-gene interaction matrix, gene expression, and protein expression
        D : str or array-like, either the TF-gene interaction matrix or the key in adata.varm containing the TF-gene interaction matrix.
        P : str or array-like, either the protein expression  matrix or the key in adata.obsm containing the protein expression interaction matrix.
        Y : str or array-like, either the gene expression  matrix or the key in adata.obsm containing the gene expression interaction matrix.
        Returns
        -------

        '''

        if adata is not None:

            self.protein_key=P
            self.gene_key=Y

            D=adata.varm[D]
            P=adata.obsm[P]
            Y=adata.layers[Y]

            self.obs_names=adata.obs_names
            self.gene_names=adata.var_names

        if type(D) is pd.DataFrame or anndata._core.views.DataFrameView:
            self.tf_names=D.columns
        if type(P) is pd.DataFrame or anndata._core.views.DataFrameView:
            self.protein_names=P.columns
        if type(Y) is pd.DataFrame:
            self.obs_names=Y.index
            self.gene_names=Y.columns

        if scipy.sparse.issparse(Y):
            Y=Y.todense()
        if scipy.sparse.issparse(D):
            D=D.todense()
        if scipy.sparse.issparse(P):
            P=P.todense()


        D=np.asarray(D)
        P=np.asarray(P)
        Y=np.asarray(Y)

        if type(self.lamda) is list or type(self.alpha) is list:
            self.lamda, self.alpha=self._cv_grid_search(D, P, Y)

        self.YTD=Y.dot(D)
        self.P=P
        self.YTY=Y.dot(Y.T)
        self.Y=Y.T
        self.D=D

        self.svdYTD=np.linalg.svd(self.YTD)
        self.svdP=np.linalg.svd(self.P)

        self._fit(self.lamda,self.alpha)

        return self

    def _cv_grid_search(self, D, P, Y):

        n_steps=self.grid_search_steps

        if self.alpha[0]==self.alpha[1]:
            alpha_list=[self.alpha[0]]
        else:
            alpha_list=np.linspace(self.alpha[0], self.alpha[1], num=n_steps)


        if self.lamda[0]==self.lamda[1]:
            lamda_list=[self.lamda_list[0]]
        else:
            lamda_list=np.geomspace(self.lamda[0], self.lamda[1], num=n_steps)


        best_lamda=[None]*self.n_folds
        best_alpha=[None]*self.n_folds
        folds=np.round(np.random.rand(P.shape[0])*self.n_folds)

        for fold in range(self.n_folds):
            #("fold: "+str(fold))
            P_train=P[folds==fold,:]
            P_test=P[~(folds==fold),:]
            Y_train=Y[folds==fold,:]
            Y_test=Y[~(folds==fold),:]

            self.YTD=Y_train.dot(D)
            self.P=P_train
            self.YTY=Y_train.dot(Y_train.T)
            self.Y=Y_train.T
            self.D=D

            self.svdYTD=np.linalg.svd(self.YTD)
            self.svdP=np.linalg.svd(self.P)

            best_score=0
            for l in lamda_list:
                for a in alpha_list:
                    self._fit(l,a)
                    s=self.score(P=P_test, Y=Y_test)
                    #print("\t l: " +str(np.round(l,3)) + "\ta: "+str(np.round(a,3)) +" \t score: "+str(s))
                    if s>best_score:
                        best_lamda[fold]=l
                        best_alpha[fold]=a
                        best_score=s

        alpha=np.mean(best_alpha)
        lamda=np.mean(best_lamda)

        return (lamda, alpha)

    def _fit(self,l,a):
        if self.solver=="SVD":
            self.W=aff_reg_svd(self.YTD, self.P.T, self.YTY,l)

        elif self.solver=="Kron":
            # svdA=self.svdYTD, svdB=self.svdP
            self.W=aff_reg_kron(np.double(self.D), np.double(self.P), np.double(self.Y),l*a,l*(1-a), spectrumP=self.spectrum)

    def ar_reconstruction(self, pred_test=None):
        """ reconstruction function
        Parameters
        ----------
        pred_test: prediction on test data

        """
        A = self.Y.T @ pred_test
        B = scipy.linalg.orth(self.Y)
        cm = scipy.linalg.lstsq(B, self.Y)[0]
        ct = scipy.linalg.lstsq(cm.T, A)[0]
        pred = B @ ct
        return pred

    def predict(self, P_test=None):
        """ predict target gene expression

        Parameters
        ----------
        P_test: Protein expression on test data

        Returns
        -------
        Y_pred: array of shape (N, Mtest)
                The predicted Y matrix on test data set which has N genes and Mtest cells
        """

        pred = self.D @ (self.W @ P_test.T)

        aff_rec = self.ar_reconstruction(pred)

        self.Y_pred = aff_rec
        return self.Y_pred.T

    def score(self,adata=None, P=None, Y=None):
        if adata is not None:
            P=adata.obsm[self.protein_key]
            Y=adata.layers[self.gene_key]
        elif P is None or Y is None:
            P=self.P
            Y=self.Y

        Y_pred=self.predict(P)

        return self.corr(Y_pred.flatten(), Y.flatten())[0]

    def get_W(self):
        # get coefficient matrix
        return self.W

    def get_projP(self, Y=None):
        """ get projected protein expression

        Parameters
        ----------
        Y:  array of shape (optional, default is (N, M))
            input gene expression with N genes and M cells

        Returns
        -------
        projP: array of shape (M, S)
               projected protein expression with M cells and S proteins

        """
        if Y is None:
            Y = self.Y
        W = self.W
        return (Y.T @ self.D @ W).T

    def get_projD(self, P=None, n_trials=0, verbose=False):
        """ get projected TF activity
        Parameters
        ----------
        P: array of shape (optional, default is (M, S) )
           input protein expression with M cells and S proteins

        Returns
        -------
        projD:  array of shape (Q, M)
            projected TF activities with Q TFs and M cells

        """

        if P is None:
            P = self.P
            obs_names=self.obs_names
        elif type(P) is pd.DataFrame or type(P) is anndata._core.views.DataFrameView:
            obs_names=P.index
            P=np.asarray(P)
        else:
            obs_names=None

        tf_activity=self.W @ P.T
        tf_activity=tf_activity.T
        if obs_names is not None and self.tf_names is not None:
            tf_activity=pd.DataFrame(tf_activity, index=obs_names, columns=self.tf_names)

        if n_trials==0:
            return tf_activity
        else:
            permuted_model=deepcopy(self)
            tf_p_val=np.zeros(tf_activity.shape)


            pbar=range(n_trials)
            if verbose:
                pbar=tqdm(range(n_trials))

            for i in pbar:
                permuted_model.Y=permuted_model.Y[np.random.permutation(self.Y.shape[0]), :]
                permuted_model.Y=permuted_model.Y[:, np.random.permutation(self.Y.shape[1])]

                permuted_model.YTD=permuted_model.Y.T.dot(permuted_model.D)
                permuted_model.svdYTD=np.linalg.svd(permuted_model.YTD)

                permuted_model._fit(self.lamda,self.alpha)
                tf_activity_perm=(permuted_model.W @ P.T).T

                tf_p_val+=(np.abs(tf_activity_perm)>np.abs(tf_activity))

            tf_p_val=(tf_p_val+1)/(n_trials+1)

            # tf_p_val=pd.DataFrame(statsmodels.stats.multitest.multipletests(
            #     np.asarray(tf_p_val).flatten(),
            #     alpha=0.15,
            #     method='holm'
            # )[1].reshape(tf_activity.shape),index=tf_activity.index,columns=tf_activity.columns)


            return tf_activity, tf_p_val
    #     def get_tf_activity(self, adata=None, P=None):
    #         pass
    def get_tf_protein_cor(self,P=None):
        if P is None:
            P=self.P
        tf_activity=self.get_projD(P)

        # X=scipy.stats.pearsonr(
        #     tf_activity,
        #     P,
        #     axis=0)[0]



        X=scipy.stats.zscore(tf_activity, axis=0).T.dot(
            scipy.stats.zscore(P, axis=0)
        )/P.shape[0]

        if self.tf_names is not None and self.protein_names is not None:
            X=pd.DataFrame(
                X,
                index=self.tf_names,
                columns=self.protein_names
            )

        return X




def aff_reg_kron(
        D,
        P,
        Y,
        l1,
        l2,
        svdA=None,
        svdB=None,
        spectrumP=0.7,
        spectrumA=1
):

    """ trains a SPaRTAN model

    Parameters
    ----------
    D : array of shape (N, Q)
        The data matrix with N genes and Q TFs

    P : array of shape (M, S)
        The data matrix with M cells and S proteins

    Y : array of shape (N, M)
        The data matrix with N genes and M cells

    lamda : float > 0, default=0.001
        LASSO regularization for linear regression

    rsL2 : float > 0, default=0.001
        ridge regularization for linear regression

    corrtype: string, default='spearman'
        correlation type used to evaluate the performance

    """


    # transformation
    A = Y.T @ D
    B = P.T
    Y = Y.T @ Y

    # SVD(A) SVD(B)
    if svdA is None:
        svdA= np.linalg.svd(A)
    if svdB is None:
        svdB=np.linalg.svd(B)

    UA, SA, VhA = svdA
    VA = VhA.T
    UB, SB, VhB = svdB
    VB = VhB.T

    a_cum_spectrum = np.cumsum(SA) / sum(SA)
    b_cum_spectrum = np.cumsum(SB) / sum(SB)

    da = np.nonzero(a_cum_spectrum >= spectrumA)[0][0] + 1
    db = np.nonzero(b_cum_spectrum >= spectrumP)[0][0] + 1

    Ua = UA[:, :da]
    Sa = SA[:da]
    Va = VA[:, :da]

    Ub = UB[:, :db]
    Sb = SB[:db]
    Vb = VB[:, :db]

    Yv = (Y.T).flatten()

    Vb = Vb.copy(order='C')
    Ua = Ua.copy(order='C')
    L = krnP.kron(Vb, Ua)

    d = np.eye(Y.shape[0], Y.shape[1])
    cidex = np.where(d.flatten() != 0)
    diag = np.array(cidex, dtype=np.int32).flatten()

    # make it c-like contiguous array
    Yv = Yv.copy(order='C')
    diag = diag.copy(order='C')

    L, Yv = krnP.removeDiagC(L, Yv, diag)

    opts = dict()
    opts['rsL2'] = l2

    # reshape Yv to 2darry
    Yv = Yv.reshape(Yv.shape[0], 1)
    beta, b = leastR.LeastR(L, Yv, l1, opts)

    del L, Yv
    gc.collect()

    Sa = np.diag(Sa)
    Sb = np.diag(Sb)

    m1 = Va
    m2 = np.linalg.pinv(Sa)
    m3 = beta.reshape(Va.shape[1], Ub.shape[1], order="F")
    m4 = np.linalg.pinv(Sb)
    m5 = Ub.T
    W = m1 @ m2 @ m3 @ m4 @ m5

    return W

def aff_reg_svd(a, b, c, l2, svdA=None, svdB=None):

    (d1, d2)=a.shape
    (d3, d4)=b.shape

    if svdA is None:
        [ua,sa,va]=np.linalg.svd(a, full_matrices=False)
    else:
        [ua,sa,va]=svdA
    if svdB is None:
        [ub,sb,vb]=np.linalg.svd(b, full_matrices=False)
    else:
        [ub,sb,vb]=svdB

    scale=np.multiply(sa.reshape((len(sa), 1)), sb.reshape((1,len(sb)))) **2
    scale=np.divide(1, scale+(d1*d4*l2))

    return va.T.dot(np.multiply(scale, va.dot(a.T.dot(c).dot(b.T)).dot(ub))).dot(ub.T)
