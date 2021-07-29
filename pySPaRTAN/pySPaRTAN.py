
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
import cythKronPlus as krnP
import cythLeastR as leastR
import scipy.linalg
import functools
import time
import gc
import matplotlib.pyplot as plt
from scipy import stats

class Timer:
    """ a class to convert time in seconds to the format of h:m:s
    
    Methods
    -------
    def __init__(self):
        initiate a timer

    def restart(self):
        restart a timer

    def get_time_hhmmss(self):
        return the period = end_time - start_time in (h, m, s) format
        
    """

    def __init__(self):
        # initiate a timer
        self.start = time.time()

    def restart(self):
        # restart a timer
        self.start = time.time()

    def get_time_hhmmss(self):
        # return the period = end_time - start_time
        # in (h, m, s) format
        end = time.time()
        m, s = divmod(end - self.start, 60)
        h, m = divmod(m, 60)
        time_str = "%02d:%02d:%02d" % (h, m, s)
        return time_str


def normalize_column(A, T=0):
    """ perform l2 normalization column-wize of given matrix
    
    Parameters:
        A : the matrix that works on
        T : switch of column-wize and row-wize.
            T=0: column-wize
            T=1: row-wize
            
    """
    if (T == 0):
        return np.divide(A, np.sqrt(np.sum(A**2, 0)))
    else:
        At = np.transpose(A)
        return np.transpose(np.divide(At, np.sqrt(np.sum(At**2, 0))))


class pySPaRTAN:
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

    def fit(self, D, P, Y, lamda=0.001, rsL2=0.001, corrtype='spearman'):
 
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
        
        spectrumA = 1
        spectrumP = 0.7

        self.D = D
        self.P = P
        self.Y = Y
        self.corrtype = corrtype

        # transformation
        A = self.Y.T @ self.D
        B = self.P.T
        Y = self.Y.T @ self.Y

        # SVD(A) SVD(B)
        UA, SA, VhA = np.linalg.svd(A)
        VA = VhA.T
        UB, SB, VhB = np.linalg.svd(B)
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
        opts['rsL2'] = 0

        # reshape Yv to 2darry
        Yv = Yv.reshape(Yv.shape[0], 1)
        beta, b = leastR.LeastR(L, Yv, lamda, opts)

        del L, Yv
        gc.collect()

        self.beta = beta
        self.Ua = Ua
        self.Ub = Ub
        self.Sa = np.diag(Sa)
        self.Sb = np.diag(Sb)
        self.Va = Va
        self.Vb = Vb
        self.lamda = lamda

    def ar_model2w(self):
        # converts a trained model to W
        m1 = self.Va
        m2 = np.linalg.pinv(self.Sa)
        m3 = self.beta.reshape(self.Va.shape[1], self.Ub.shape[1], order="F")
        m4 = np.linalg.pinv(self.Sb)
        m5 = self.Ub.T
        ww = m1 @ m2 @ m3 @ m4 @ m5
        return ww

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
        if P_test is not None:
            self.P_test = P_test

        w = self.ar_model2w()
        pred = self.D @ (w @ self.P_test.T)

        aff_rec = self.ar_reconstruction(pred)

        self.Y_pred = aff_rec
        return self.Y_pred

    def get_corr(self, Y_pred, Y_test, plot=False):
        """ get the correlation between predicted Y_pred and Y_test
        
        Parameters
        ----------
        Y_pred: array of shape (N, Mtest)
                predicted gene expression with N genes and Mtest cells
                
        Y_test: array of shape (N, Mtest)
               gene expression test data with N genes and Mtest cells
        plot: whether to plot the correlation between Y_pred and Y_test, default is False
        
        
        Returns
        -------
        corr: float 0 <= value <= 1
              spearman/pearson corrlatioin between flattened Y_pred and Y_test
        
        """
        if self.corrtype == 'spearman':
            corr = stats.spearmanr(Y_test.ravel(order='F'), Y_pred.ravel(order='F'))[0]
        else:
            corr = stats.pearsonr(Y_test.ravel(order='F'), Y_pred.ravel(order='F'))[0]                           

        if plot:
            plt.plot(Y_test.ravel(order='F'), Y_pred.ravel(order='F'),
                     linestyle='none', marker='+')
            plt.title('reconstruction of Y test, corr={:.2f}'.format(corr))

        return corr

    def get_W(self):
        # get coefficient matrix
        
        self.W = self.ar_model2w()
        return self.W

    def get_projP(self, Y=None):
        """ get projected protein expression
        
        Parameters
        ----------
        Y:  array of shape (optional, default is (N, M) )
            input gene expression with N genes and M cells
            
        Returns
        -------
        projP: array of shape (M, S)
               projected protein expression with M cells and S proteins
               
        """
        if Y is None:
            Y = self.Y
        W = self.ar_model2w()
        return (Y.T @ self.D @ W).T

    def get_projD(self, P=None):
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
        W = self.ar_model2w()
        return W @ P.T
