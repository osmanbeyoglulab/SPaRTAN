# pySPaRTAN

pySPaRTAN is a  python implementation of the SPaRTAN pipeline (Single-cell Proteomic and RNA based Transcription factor Activity Network) which enables users to infer transcription factor (TF) activities and link cell-surface receptors to TFs by exploiting cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq) datasets with cis-regulatory information.

## Data
In this demo, we use a subset of the CITE-seq data for 5k (Nextgen) PBMC obtained from 10X Genomics website data to train SPaRTAN.
To construct the TF – target gene prior matrix, we downloaded a gene set resource containing TF target-gene interactions from DoRothEA. 
Sample datasets used in our study can be downloaded from http://www.pitt.edu/~xim33/SPaRTAN 

## System requirements
The code was generated in python 3.7. Besides python installation, some other packages from PyPI repository, such as pandas, numpy, sklearn, scipy, and matplotlib need to be installed additionaly. In order to improve the running time performance, we converted some computationally intensive python functions into two Cython modules which are platform dependent. We have built those Cython extensions for Window(.pyd files), Mac(.so files), and Linux system(.so files) respectivly. Cython extensions can be built onsite as well which is explained in the section "Cython extension compilation". 

## Installation
Download the reporsitory from https://github.com/osmanbeyoglulab/SPaRTAN

Install python3.7 or later. pySpaRTAN used the following dependencies as well: pandas, numpy, scipy, sklearn, matplotlib. 

You can install python dependencies by the following commands:
```sh
pip install pandas
pip install numpy
pip install scipy
pip install -U scikit-learn
pip install matplotlib
```
Cython is not reqired to be installed unless the pre-built Cython extensions do not work on your system. 

Our codes have been tested on Linux, Mac, and Windows systems. Please see Prerequisites.xlsx for the version of packages we tested on each operating system.

### Cython extension compilation

If the Cython extensions, which are platform-dependent binary modules, are not compatible with your operating system, additional built of those Cython extensions are needed on your own machine. 

First, install Cython by running the command
```sh
pip install "cython>0.21"    
```
Cython requires a C compiler to be present on the system. Please check [here](https://cython.readthedocs.io/en/latest/src/quickstart/install.html) for C compiler installation on various operating systems

After installing Cython and C compiler, navigate to the directory pySPaRTAN and type in the command:
```sh
python setup.py build_ext --inplace
```
This will generate new Cython extension .so files (or .pyd files on Windows). The previously downloaded .so and .pyd files are renamed to "*_old.so" and "*_old.pyd" 

## Usage

pySPaRTAN module has 3 input datasets, and 2 generated ouputs

**Input datasets**
```sh
D: dataframe of shape (N, Q)
   The dataframe with N genes and Q TFs
   
P: dataframe of shape (M, S)
   The dataframe with M cells and S proteins  
   
Y: dataframe of shape (N, M)
   The dataframe with N genes and M cells    
```
**Outputs**
```sh
projD: dataframe of shape (Q, M) 
       projected TF activities with Q TFs and M cells 
       
projP: dataframe of shape (S, M)
      projected protein expression with S proteins and M cells  
```

**Hyperparameters**

pySPaRTAN model has 2 hyperparameters: lamda and rsL2

    lamda : float > 0, default=0.001
            LASSO regularization for linear regression 
   
    rsL2 : float > 0, default=0.001
           ridge regularization for linear regression

We can run pySPaRTAN by giving some specific values to those hyperparameters or using default setting in the script.
We can also use cross-validation to find the optional values of those hyperparameters for the input datasets, and then run pySPaRTAN to generate the projected outputs.

**Command line**

To run pySPaRTAN module with default parameters, simply execute the following command to generate outputs 
```sh
python run_pySPaRTAN.py
```

To check all optional arguments, type in the command line
```sh
python run_pySPaRTAN.py -h
```
It shows each argument and its description:

    optional arguments:
      
      --input_dir INPUT_DIR
                             string, default='../data/inputs'
                             Directory of input files
      --output_dir OUTPUT_DIR
                             string, default: '../data/outputs'
                             Directory of output files
      --dataset_D DATASET_D
                             string, default='Dpbmc'
                             File name of (gene X TF) dataframe.
                             Requires .csv format,
                             only contains file name, not include '.csv' extension
      --dataset_P DATASET_P
                             string, default='Ppbmc5kn_CD16'
                             File name of (cell X protein) dataframe.
                             Requires .csv format,
                             only contains file name, not include '.csv' extension
      --dataset_Y DATASET_Y
                             string, default='Ypbmc5kn_CD16'
                             File name of (gene X cell) dataframe.
                             Requires .csv format,
                             only contains file name, not include '.csv' extension
      --lamda LAMDA         
                             float, value>0, default=0.001
                             LASSO regularization for linear regression.
      --rsL2 RSL2           
                             float, value>0, default=0.001
                             Ridge regularization for linear regression,
      --normalization NORMALIZATION
                             string, default='l2'
                             Type of normalization performed on matrices,
                             if set to empty(''), then no normalization
      --fold FOLD           
                             int, value>=0, default=0
                             How many folds for the cross_validation.
                             value=0 means no cross_validation and using default/specified parameters
      --correlation
                             string, ['pearson', 'spearman'], default='pearson'
                             Type of correlation coefficient
                             
For example, we can use the following command to run the pySPaRTAN model with 5-fold cross-validation and using spearman correlation:

    python run_pySPaRTAN.py --fold 5 --correlation spearman

