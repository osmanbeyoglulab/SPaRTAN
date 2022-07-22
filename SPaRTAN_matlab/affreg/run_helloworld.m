%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       a very simple simulation
%

function run_helloworld
% run affinity regression on a simple simulation

n = 500;
m = 160;
p = 120;
q = 90;
n_test = 10;

D = randn(n, p);
P = randn(m, q);
P_test = P(1:n_test, :);
P_train = P(n_test+1:end, :);
W = randn(p, q);
Y_test  = D*W*P_test';
Y_train = D*W*P_train';

lambda = 0.1;
rsL2 = 0;
spectrumA = 0.95;
spectrumB = 0.9;

% train a model, predict on test data P_test, and reconstruct W
tic
model  = ar_train(D, P_train, Y_train, lambda, rsL2, spectrumA, spectrumB);
toc
Y_pred = ar_predict(D, P_test, Y_train, model);
W_pred = ar_model2w(model);

w_small_pred = model.fit_svd.beta;
w_small = model.Sa * model.Va' * W * model.Ub * model.Sb;

figure, plot(W(:), W_pred(:), 'o'), title(sprintf('reconstruction of W, corr=%.2f', corr(W(:), W_pred(:))))
figure, plot(w_small(:), w_small_pred(:), 'o'), title(sprintf('prediction of truncated \tilde w, corr=%.2f', corr(w_small(:), w_small_pred(:))))
figure, plot(Y_test(:), Y_pred.rec(:), 'o'), title(sprintf('reconstruction of Y test, corr=%.2f', corr(Y_test(:), Y_pred.rec(:))))

fprintf('for more extensive SOM simulation run (run_simulation.m) or can be found at http://bitbucket.org/leslielan/affreg\n')



