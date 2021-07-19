%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function [predictions, idx] = predictBL(seqs_test, seqs_train, Y_train)
% find the nearest blosum neighbor in input space P
similarity = nwalign(seqs_test, seqs_train);
[~, idx] = max(similarity,[],2);
predictions = Y_train(:,idx);
