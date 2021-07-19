%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       train the affinityregression model
%

function [model] = ar_train(A,B,Y,lambda, rsL2, spectrumA, spectrumB, old_version)
% Rafi Pelossof 2015, MSKCC
%   ar_train trains an affinity regresssion model
%   A - left matrix: measurements x features
%   B - right matrix: samples x features
%   Y - output matrix: measurements x samples
%   lambda - L2 regularization
%   rsL2 - L1 regularization, recommended rsL2 = 0
%   spectrumA - percent spectrum to keep (0,1] for left matrix
%   spectrumB - percent spectrum to keep (0,1] for right matrix
%   affinity_transform - 
%       0: solve AWB'=Y (for small problems)
%       1: solve Y'AWB'=Y'Y (for big problems)
%
%   old_version
%       0 - remove diagonal equations (from Y'Y)
%       1 - keep diagonal equations (from Y'Y)

if nargin < 8
    old_version = 0;
end

%ar transformation
A = Y'*A;
B = B';
Y = Y'*Y;

% solution through lower dimensional SVD(A), SVD(B)
[UA SA VA] = svd(A);
[UB SB VB] = svd(B);
%da = min(size(SA));

a_spectrum = diag(SA);
b_spectrum = diag(SB);
a_cum_spectrum = cumsum(a_spectrum)/sum(a_spectrum);
b_cum_spectrum = cumsum(b_spectrum)/sum(b_spectrum);

da = find(a_cum_spectrum >= spectrumA, 1 , 'first');
db = find(b_cum_spectrum >= spectrumB, 1 , 'first');
%da = length(a_cum_spectrum);

Ua = UA(:,1:da); Sa = SA(1:da,1:da); Va = VA(:,1:da);
Ub = UB(:,1:db); Sb = SB(1:db,1:db); Vb = VB(:,1:db);

Yv = Y(:);
L = kron(Vb,Ua);

% remove elements from the diagonal
if ~old_version
    d = eye(size(Y));
    Yv(find(d)) = [];
    L(find(d),:) = [];
end

opts.rsL2=rsL2;
fprintf('ar_train.m: running LeastR L=(%d x %d), Yv=%d, lambda=%.3f, rsL2=%.3f, specA=%.3f, specB=%.3f\n', size(L,1), size(L,2), length(Yv), lambda, opts.rsL2, spectrumA, spectrumB);
start = tic();
[beta, b] = LeastR(L, Yv, lambda, opts);
end_time = toc(start);

model.Ua = Ua;
model.Ub = Ub;
model.Sa = Sa;
model.Sb = Sb;
model.Va = Va;
model.Vb = Vb;
model.fit_svd.beta = beta;
model.fit_svd.lambda = lambda;
model.fit_svd.pred = b;
model.end_time = end_time;


