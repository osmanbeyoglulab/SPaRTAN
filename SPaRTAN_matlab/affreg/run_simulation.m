%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate simulation figure S23
%

function run_simulation

% this function runs a simulation of affinity regression
% you need to install SLEP to have the code run:
% http://www.public.asu.edu/~jye02/Software/SLEP/
%
% The code creates multiple models for different parameter configurations
% and plots performance for the best achievable solution.


version = 'sim';

n = 50; m = 50;
p = 30; q = 30;

dn = randn(n,p);
pn = randn(m,q);
w = randn(p, q);
y_sim = dn*w*pn';
y_sim = zscore(y_sim);
y_sim = bsxfun(@rdivide, y_sim, sqrt(sum(y_sim.^2))); %unit norm


lambdas = [0.001 0.01 0.05 0.1];
rsL2 =    [0.001 0.01 0.05 0.1];
spectrumA = [1 0.9 0.8 0.7 0.6];
spectrumB = [1 0.9 0.8 0.7 0.6];

[L1 L2 SA SB] = ndgrid(lambdas, rsL2, spectrumA, spectrumB);
params = [L1(:) L2(:) SA(:) SB(:)];

for param_ix = 1:size(params,1)
    lambda = params(param_ix, 1);
    rsL2 = params(param_ix, 2);
    spectrumA = params(param_ix, 3);
    spectrumB = params(param_ix, 4);
    fprintf('%d: starting %d/%d\n',param_ix, param_ix, size(params,1));
    
    models_par{param_ix} = ar_train(dn, pn,y_sim,lambda, rsL2, spectrumA, spectrumB);
    new_pred = ar_predict(dn, pn, y_sim, models_par{param_ix});
    w_pred = ar_model2w(models_par{param_ix});
    cw(param_ix) = corr(w(:), w_pred(:));
    cy(param_ix) = mean(diag(corr(y_sim, new_pred.rec)));
    fprintf('%d: finished running cw=%f, cy=%f, elapsed time: %f secs\n', param_ix, cw(param_ix), cy(param_ix),models_par{param_ix}.end_time);
end

% pick optimal model
[m, max_ix] = max(cy);
model = models_par{max_ix};
pred = ar_predict(dn, pn, y_sim, model);

w_pred = ar_model2w(model);


w_v = w(:);
w_pred_v = w_pred(:);
w_small = model.Sa * model.Va' * w * model.Ub * model.Sb;
w_small_pred = model.fit_svd.beta;

%
% Plots
%

% create a nearest neighbor
pc = corr(pn');
pc = pc - eye(size(pc));
[m mix] = max(pc);
ynn = y_sim(:, mix);

cc_rec = diag(corr(y_sim, pred.rec));
cc_nn = diag(corr(y_sim, ynn));

clrred = [228,26,28]/255;
clrblue = [55,126,184]/255;
markersize = 5;

% plot w corr vs y performance
h = figure;
plot(cw, cy,'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
xlabel('W correlation')
ylabel('Y correlation')
hold on;
[lx, ly, rsq] = plotfit(cw, cy, clrred);
axis square
title(sprintf('Relationship between model accuracy and prediction accuracy, R^2=%f', rsq))
print(h, sprintf('results/sim_%s_cw_cy.pdf',version), '-dpdf');

% plot corr in pred y against nn
h = figure;
plot(cc_nn, cc_rec, 'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
line([0 1],[0 1],'color', clrred)
axis square
axis([0 1 0 1])
xlabel('Nearest neighbor')
ylabel('Affinity regression')
title('Probe intensity prediction correlation')
print(h, sprintf('results/sim_%s_cc_nn_cc_rec.pdf',version), '-dpdf');

% plot w vs w_pred
ix = randsample(numel(w), min(numel(w), 10000));
h = figure;
plot( w_v(ix), w_pred_v(ix),'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
hold on;
[lx, ly, rsq] = plotfit(w_v(ix), w_pred_v(ix), clrred);
axis square
xlabel('W simulated')
ylabel('W predicted')
title(sprintf('Prediction of full model W, R^2=%f', rsq))
print(h, sprintf('results/sim_%s_w.pdf',version), '-dpdf');

% plot w_small vs pred
h = figure;
plot( w_small(:), w_small_pred(:), 'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
hold on;
[lx, ly, rsq] = plotfit(w_small(:), w_small_pred(:), clrred);
xlabel('W tilde simulated')
ylabel('W tilde predicted')
axis square
title(sprintf('Prediction of compressed model W tilde, R^2=%f', rsq))
print(h, sprintf('results/sim_%s_w_small.pdf',version), '-dpdf');

% plot kmer enrichment vs pred
proj = y_sim' * dn;
proj_pred = pred.rec' * dn;
proj_v = proj(:);
proj_pred_v = proj_pred(:);
ix = randsample(numel(proj), min(numel(proj), 10000));
h = figure;
plot(proj_v(ix), proj_pred_v(ix), 'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
hold on;
[lx, ly, rsq] = plotfit(proj_v(ix), proj_pred_v(ix), clrred);
axis square
xlabel('simulated Kmer enrichment')
ylabel('predicted Kmer enrichment')
title(sprintf('kmer enrichment Y''D, R^2 = %f', rsq))
print(h, sprintf('results/sim_%s_kmerenrich.pdf',version), '-dpdf');

% plot left_proj vs pred
proj = y_sim' * dn * w;
proj_pred = pred.rec' * dn * w_pred;
proj_v = proj(:);
proj_pred_v = proj_pred(:);
ix = 1:length(proj_v);%randsample(numel(proj), 10000);
h = figure;
plot(proj_v(ix), proj_pred_v(ix), 'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
hold on;
[lx, ly, rsq] = plotfit(proj_v(ix), proj_pred_v(ix), clrred);
axis square
xlabel('simulated AA projection')
ylabel('predicted AA projection')
title(sprintf('AA kmer projection Y''DW, R^2 = %f', rsq))
print(h, sprintf('results/sim_%s_aa_proj.pdf',version), '-dpdf');

% plot right_proj vs pred
proj = w * pn';
proj_pred = w_pred * pn';
proj_v = proj(:);
proj_pred_v = proj_pred(:);
ix = randsample(numel(proj), min(numel(proj), 10000));
h = figure;
plot(proj_v(ix), proj_pred_v(ix), 'o', 'color', clrblue, 'markersize',markersize, 'markerface', clrblue);
hold on;
[lx, ly, rsq] = plotfit(proj_v(ix), proj_pred_v(ix), clrred);
axis square
xlabel('simulated NT projection')
ylabel('predicted NT projection')
title(sprintf('DNA kmer projection WP'', R^2 = %f', rsq))
print(h, sprintf('results/sim_%s_nn_proj.pdf',version), '-dpdf');



1;

function [lx, ly, rsq] = plotfit(x, y, color)
[lx, ly, rsq] = linfit(x, y);
line (lx, ly, 'color', color) % Plot best fit line 


function [lx, ly, rsq, pPoly] = linfit(x, y)
pPoly = polyfit(x, y, 1); % Linear fit of xdata vs ydata 
lx = [min(x) max(x)]; % find left and right x values 
ly = polyval(pPoly,[min(x), max(x)]); % find y values 
%Compute the residual values as a vector signed numbers
yfit = polyval(pPoly,x);
yresid = y - yfit;
%Square the residuals and total them obtain the residual sum of squares:
SSresid = sum(yresid.^2);
%Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);
%Compute R2 using the formula given in the introduction of this topic:
rsq = 1 - SSresid/SStotal;
