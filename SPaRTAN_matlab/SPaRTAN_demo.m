%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SPaRTAN (Single-cell Proteomic and RNA based Transcription factor
%  Activity Network) trains on parallel gene mRNA and surface protein
%  expression data from a set of cells (e.g. based on CITE-seq),
%  together with curated transcription factor (TF)-target gene interactions,
%  to learn a model that links surface proteins and TFs (D).
%  Specifically, the algorithm learns a weight matrix W between
%  surface proteins and TFs that predicts gene expression (Y)
%  from cell-specific surface protein levels (P) and
%  TF-target gene interactions (D) by solving a bilinear regression framework
%  called affinity regression.
%  Using the trained interaction matrix (W), one can predict TF activities
%  from the cell surface protein expression profile (WP'),
%  or protein expression from the cellular mRNA expression data
%  and the TF-target gene interaction matrix (Y'DW).
%  This demo illustrates
%  -how to build and assess SPaRTAN models
%  -how to predict cell-specific TF activities
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close;
addpath(genpath('SLEP_package_4.1/'));
addpath(genpath('affreg/'));

stream = RandStream.getGlobalStream;
reset(stream,100);

%%%Read in the sample inputs from different files
% Transcription factor (TF)-target gene interactions
D = readtable("../data/inputs/Dpbmc.csv", 'ReadRowNames',true);
% CLR normalized surface protein expression matrix
P = readtable("../data/inputs/Ppbmc5kn_CD16.csv", 'ReadRowNames',true);
% log-normalized Gene expression matrix
Y = readtable("../data/inputs/Ypbmc5kn_CD16.csv", 'ReadRowNames',true);

D = table2array(D); %gene by TF
P = table2array(P); %protein by cell
Y = table2array(Y); %gene by cell

%%Normalize data: D,Y, P were unit-normalized for computational stability
D_norm = normalize_column(D);
P_norm = normalize_column(P')';
Y_norm = normalize_column(Y);

%learning setting
lambda = 0.001; % L1 regularization parameter
rsL2 = 0.001;   % L2 regularization parameter
spectrumD = 1;
spectrumP = 0.7;

% affinity regression train and predict
model = ar_train(D_norm,P_norm,Y_norm,lambda,rsL2,spectrumD,spectrumP);
pred = ar_predict(D_norm,P_norm,Y_norm,model);

Y_pred = [];
corr_type = 'spearman';
for i = 1:size(Y_norm,2)
    Y_pred(i) = corr(pred.rec(:, i), Y_norm(:,i),'type',corr_type);
end

% W not uniquely defined "projections" onto phosphoprotein space (P) and TF space (D) are more meaningful
W = ar_model2w(model); %interaction matrix TF by surface protein

% Projections
% inferred surface protein expression
projP = Y_norm' * D_norm * W;
% inferred transcription factor activity
projD = W * P_norm';

% correlation between TF activity and surface protein expression
pcc_TF_P = corr(projD', P_norm,'type','Pearson');
