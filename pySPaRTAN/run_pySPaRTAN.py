"""
This script intends to use pySPaRTAN module to generate predicted matrices used in the paper.

Input data
----------
D: dataframe of shape (N, Q)
   The data frame with N genes and Q TFs
   
P: dataframe of shape (M, S)
   The data frame with M cells and S proteins  
   
Y: array of shape (N, M)
   The data frame with N genes and M cells
   
location: default in the directory ../data/inputs
   
Output data
-----------
projD: dataframe of shape (Q, M) 
       projected TF activities with Q TFs and M cells 
       
projP: dataframe of shape (S, M)
      projected protein expression with S proteins and M cells    
     
location: location: default in the directory ../data/outputs
 
Hyperparameters
---------------
      
pySPaRTAN has 2 Hyperparameters that can be adjusted: lamda and rsL2
We can run pySPaRTAN by specifying some values to those parameters or using default ones in the script.
We can also use cross-validation at first to get the optional values for those
hyperparameters, and then run pySPaRTAN to generate the projections.

Command lines
-------------
When running this script from command line, the following parameters can be added to the command:
    
    --input_dir : directory of input files, default="../data/inputs"
    
    --output_dir : directory of output files, default="../data/outputs"
    
    --dataset_D : dataframe of gene X TF 
                  Requires .csv format, only contains file name, not include ".csv" extension
                  
    --dataset_P : dataframe of cell X protein
                  Requires .csv format, only contains file name, not include ".csv" extension
        
    --dataset_Y : dataframe of gene X cell
                  Requires .csv format, only contains file name, not include ".csv" extension
                  

    --lamda :  float > 0.0, default=0.001
            LASSO regularization for linear regression 
            
    --rsL2 : float > 0.0, default=0.001
            ridge regularization for linear regression
            
    --normalization : string, default = "l2"
                     type of normalizion performed on matrices,
                     if set to "", then no normalization
                     
    --fold : int >=0, default = 0
             how many folds to be used when doing cross-validation.
             if set to 0, it means using default/specified hyper-parameters, 
             do not conduct cross-validation
 
    --correlation : string, ("pearson" or "spearman") default = "pearson"
                    type of correlation coefficient
                     
System requirements
------------------                         
This script requires numpy, pandas, sklearn to be installed in the python running environment

"""
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import normalize
import numpy as np
import cythKronPlus as krnP
import cythLeastR as leastR
import scipy.linalg
import gc
# import matplotlib.pyplot as plt
from scipy import stats
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
        opts['rsL2'] = rsL2

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

    def get_corr(self, Y_pred, Y_test):#, plot=False):
        """ get the correlation between predicted Y_pred and Y_test

        Parameters
        ----------
        Y_pred: array of shape (N, Mtest)
                predicted gene expression with N genes and Mtest cells

        Y_test: array of shape (N, Mtest)
               gene expression test data with N genes and Mtest cells
       # plot: whether to plot the correlation between Y_pred and Y_test, default is False


        Returns
        -------
        corr: float 0 <= value <= 1
              spearman/pearson corrlatioin between flattened Y_pred and Y_test

        """
        if self.corrtype == 'spearman':
            corr = stats.spearmanr(Y_test.ravel(order='F'), Y_pred.ravel(order='F'))[0]
        else:
            corr = stats.pearsonr(Y_test.ravel(order='F'), Y_pred.ravel(order='F'))[0]

        #         if plot:
        #             plt.plot(Y_test.ravel(order='F'), Y_pred.ravel(order='F'),
        #                      linestyle='none', marker='+')
        #             plt.title('reconstruction of Y test, corr={:.2f}'.format(corr))

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
parser = ArgumentParser(description='test', formatter_class=RawTextHelpFormatter)

parser.add_argument("--input_dir", help="string, default='../data/inputs'\n\
Directory of input files",
                    type=str, default="../data/inputs")
parser.add_argument("--output_dir", help="string, default: '../data/outputs'\n\
Directory of output files",
                    type=str, default="../data/outputs")
parser.add_argument("--dataset_D", help="string, default='Dpbmc'\n\
File name of (gene X TF) dataframe.\n\
Requires .csv format,\n\
only contains file name, not include '.csv' extension",
                    type=str, default="Dpbmc")
parser.add_argument("--dataset_P", help="string, default='Ppbmc5kn_CD16'\n\
File name of (cell X protein) dataframe.\n\
Requires .csv format,\n\
only contains file name, not include '.csv' extension",
                    type=str, default="Ppbmc5kn_CD16")
parser.add_argument("--dataset_Y", help="string, default='Ypbmc5kn_CD16'\n\
File name of (gene X cell) dataframe.\n\
Requires .csv format,\n\
only contains file name, not include '.csv' extension",
                    type=str, default="Ypbmc5kn_CD16")
parser.add_argument("--lamda", help="float, value>0, default=0.001\n\
LASSO regularization for linear regression.",
                    type=float, default=0.001)
parser.add_argument("--rsL2", help="float, value>0, default=0.001\n\
Ridge regularization for linear regression,",
                    type=float, default=0.001)
parser.add_argument("--normalization", help="string, default='l2'\n\
Type of normalizion performed on matrices,\n\
if set to empty(''), then no normalization ",     
                    type=str, default="l2")
parser.add_argument('--fold', help="int, value>=0, default=0\n\
How many folds for the cross_validation.\n\
value=0 means no cross_validation and using default/specified parameters",
                    type=int, default=0)
parser.add_argument("--correlation", help="string, value='pearson' or 'spearman'\n\
default='pearson'\n\
Type of correlation coefficient",     
                    type=str, default="pearson")

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

print("Read in datasets D, P, and Y ...")
D_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_D+'.csv'), index_col=0)
P_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_P+'.csv'), index_col=0)
Y_ori = pd.read_csv(os.path.join(args.input_dir, args.dataset_Y+'.csv'), index_col=0)

TF_name = list(D_ori.columns)
cell_name = list(Y_ori.columns)
gene_name = list(Y_ori.index)
protein_name = list(P_ori.columns)

D_mat = D_ori.values
P_mat = P_ori.values
Y_mat = Y_ori.values

# normalize the dataset
if args.normalization != "":
    D_mat = normalize(D_mat, norm=args.normalization, axis=0)
    Y_mat = normalize(Y_mat, norm=args.normalization, axis=0)
    P_mat = normalize(P_mat, norm=args.normalization, axis=1)

# create the object of SPaRTAN
reg = pySPaRTAN()

# cross-validate to determine optimal parameters
fold = args.fold
if fold != 0:  # using cross validation to determine the optimal parameters
    
    lamdas = [0.001, 0.01, 0.1, 0.2, 0.3]
    rsL2s = [0.001, 0.01, 0.1]

    lenlamdas = len(lamdas)
    lenrsL2s = len(rsL2s)
    
    corr_all = np.zeros((lenlamdas, lenrsL2s))

    for l in range(0, lenlamdas):
        for r in range(0, lenrsL2s):
            print("cross validating lambda={}, rsL2={}"
                  .format(lamdas[l], rsL2s[r]))


            Y_pred_all = np.zeros(Y_mat.shape)
            kf = KFold(n_splits=fold)
            for train_index, test_index in kf.split(P_mat):

                # split the data into train and test set
                P_train, P_test = P_mat[train_index, :], P_mat[test_index, :]
                Y_train, Y_test = Y_mat[:, train_index], Y_mat[:, test_index]

                # normalize the train and test set
                if args.normalization != "":
                    Y_train = normalize(Y_train, axis=0)
                    Y_test = normalize(Y_test, axis=0)

                    P_train = normalize(P_train, axis=1)
                    P_test = normalize(P_test, axis=1)

                    
                # train the model
                reg.fit(D_mat, P_train, Y_train, lamdas[l], rsL2s[r], args.correlation)

                # get predicted value Y_pred  on P_test
                Y_pred = reg.predict(P_test)

                # save Y_pred to the whole matrix
                Y_pred_all[:,test_index] = Y_pred


            corr = reg.get_corr(Y_pred_all, Y_mat)
            corr_all[l, r] = corr

    # retrive the best parameters
    max_l, max_r = np.unravel_index(
        corr_all.argmax(), corr_all.shape
    )
    
    lamda_best = lamdas[max_l]
    rsL2_best = rsL2s[max_r]
       
    print("lamda_best={}, rsL2_best={}"
          .format(lamda_best, rsL2_best))

else:  # fold ==0: using default/specified paramters

    lamda_best = args.lamda
    rsL2_best = args.rsL2

print("Processing ...")
# re-train the model
reg.fit(D_mat, P_mat, Y_mat, lamda_best, rsL2_best, args.correlation)

# retrieve projD, projP
projD = reg.get_projD()
projP = reg.get_projP()

df_projP = pd.DataFrame(data=projP, index=protein_name, columns=cell_name)
df_projD = pd.DataFrame(data=projD, index=TF_name, columns=cell_name)

set_name = args.dataset_P[1:]
outfile_projP = os.path.join(args.output_dir, "projP_{}_{}_{}_{}.csv".format
                (set_name, lamda_best, rsL2_best,args.correlation))
outfile_projD = os.path.join(args.output_dir, "projD_{}_{}_{}_{}.csv".format
                (set_name, lamda_best, rsL2_best, args.correlation))

df_projP.to_csv(outfile_projP)
df_projD.to_csv(outfile_projD)

print("Process finished successfully!")
