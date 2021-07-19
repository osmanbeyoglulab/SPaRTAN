%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function [predictions, idx] = predictBN(Y_test, Y_train)
% find the nearest neighbor in output space Y
idx = knnsearch(Y_train', Y_test');
predictions = Y_train(:,idx);
