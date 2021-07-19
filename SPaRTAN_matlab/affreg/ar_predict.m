%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       predicts on new data with the ar_train model
%

function predictions= ar_predict(D, P_test, Y_train, model)
% D - left matrix
% P_test - test right vector or matrix
% Y_train - Y used in training
% model - the learned model from ar_train
%
%  output: 
%   .rec - reconstruction through the span(Y_train)
%   .pred - direct prediction of Y
%   .train_aff - predictions of affinities Y^T*Y

n_probes = size(Y_train, 1);
n_test = size(P_test, 1);

predictions.rec = zeros([n_probes, n_test]);

w = ar_model2w(model);
pred = D * (w * P_test'); 
train_aff = Y_train' * pred;

aff_rec = ar_reconstruction(Y_train, pred);

predictions.rec(:,:) = aff_rec;
predictions.pred(:,:) = pred;
predictions.train_aff(:,:) = train_aff;    



