%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function [predictions, idx] = predictNN(P_test, P_train, Y_train)
% find the nearest neighbor in input space P
idx = knnsearch(P_train, P_test);
predictions = Y_train(:,idx);