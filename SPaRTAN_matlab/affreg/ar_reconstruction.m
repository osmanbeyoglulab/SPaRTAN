%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       reconstruction function as described in som 2.1
%

function pred = ar_reconstruction(y_train, pred_test)
%    pred_mult = D * model.Va * w * P_test';
%    A = y_train' * pred_mult;
A = y_train' * pred_test;

O = orth(y_train);
c = O\y_train;
ct = c'\A; %equivalent to: inv(c')*A;

pred = O*ct;

