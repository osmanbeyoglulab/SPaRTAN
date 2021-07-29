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
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import normalize
from pySPaRTAN import pySPaRTAN

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
    
    corr_all_spearman = np.zeros((lenlamdas, lenrsL2s))

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


            corr_spearman = reg.get_corr(Y_pred_all, Y_mat)
            corr_all_spearman[l, r] = corr_spearman

    # retrive the best parameters
    max_l, max_r = np.unravel_index(
        corr_all_spearman.argmax(), corr_all_spearman.shape
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
