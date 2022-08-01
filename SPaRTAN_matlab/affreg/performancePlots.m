%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function folder =  performancePlots(ar_rec, nn_rec, bl_rec, or_rec, Y_gr, source, scriptname ,dumppred, corr_type)

if nargin < 8
    dumppred = 1;
end

if ~exist('results','dir')
    mkdir('results')
end

current_time = datestr(now,'yyyymmdd_HHMMSS');
folder = sprintf('results/%s_%s', current_time, source);
mkdir (folder)
fprintf('saving in %s\n', folder)

[pathstr,name,ext] = fileparts(scriptname);
copyfile(scriptname, sprintf('%s/%s%s',folder, name, ext));

cc_rec = zeros(size(Y_gr,2),1);
cc_nn  = zeros(size(Y_gr,2),1);
cc_bl  = zeros(size(Y_gr,2),1);
cc_or  = zeros(size(Y_gr,2),1);

Top_rec = zeros(size(Y_gr,2),1);
Top_nn  = zeros(size(Y_gr,2),1);
Top_bl  = zeros(size(Y_gr,2),1);
Top_or  = zeros(size(Y_gr,2),1);

for i = 1:size(Y_gr,2);
    cc_rec(i) = corr(ar_rec(:, i), Y_gr(:,i),'type',corr_type);
    cc_nn(i) =  corr(nn_rec(:, i), Y_gr(:,i),'type',corr_type);
    cc_bl(i) =  corr(bl_rec(:, i), Y_gr(:,i),'type',corr_type);
    cc_or(i) =  corr(or_rec(:, i), Y_gr(:,i),'type',corr_type);
    [s, ix]=sort(Y_gr(:,i));
    gr_top(:,i) = ix((end-100+1):end);
    [s, ix]=sort(ar_rec(:, i));
    ar_top(:,i) = ix((end-100+1):end);
    [s, ix]=sort(nn_rec(:, i));
    nn_top(:,i) = ix((end-100+1):end);
    [s, ix]=sort(bl_rec(:, i));
    bl_top(:,i) = ix((end-100+1):end);
    [s, ix]=sort(or_rec(:, i));
    or_top(:,i) = ix((end-100+1):end);
    Top_rec(i)  = length(intersect(ar_top(:,i), gr_top(:,i)));
    Top_nn(i)   = length(intersect(nn_top(:,i), gr_top(:,i)));
    Top_bl(i)   = length(intersect(bl_top(:,i), gr_top(:,i)));
    Top_or(i)   = length(intersect(or_top(:,i), gr_top(:,i)));
end

Y_norm = quantilenorm(zscore(Y_gr));
Test_top = (Y_norm > prctile(Y_norm(:,1), 99,1)) + 0;
parfor ix = 1:size(Test_top, 2)   % use parfor for faster calculation
    fprintf('%d/%d\n',ix, size(Test_top, 2)) 
    [~,~,~, pr_rec(ix)] = perfcurve(Test_top(:,ix), ar_rec(:, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, pr_nn(ix) ] = perfcurve(Test_top(:,ix), nn_rec(:, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, pr_bl(ix) ] = perfcurve(Test_top(:,ix), bl_rec(:, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, pr_or(ix) ] = perfcurve(Test_top(:,ix), or_rec(:, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
%    [~,~,~,~, pr_rec1(ix)] = prec_rec(ar_rec(:, ix), Test_top(:,ix), 'numThresh',length(Test_top(:,1))/30);
%    [~,~,~,~, pr_nn1(ix) ] = prec_rec(nn_rec(:, ix), Test_top(:,ix), 'numThresh',length(Test_top(:,1))/30);
%    [~,~,~,~, pr_bl1(ix) ] = prec_rec(bl_rec(:, ix), Test_top(:,ix), 'numThresh',length(Test_top(:,1))/30);
%    [~,~,~,~, pr_or1(ix) ] = prec_rec(or_rec(:, ix), Test_top(:,ix), 'numThresh',length(Test_top(:,1))/30);
end


Test_bot = double(Y_norm < prctile(Y_norm(:,1), 50,1));
parfor ix = 1:size(Test_top, 2)   % use parfor for faster calculation
    fprintf('%d/%d\n',ix, size(Test_top, 2)) 
    probes = Test_bot(:, ix) | Test_top(:,ix);
    [~,~,~, prtb_rec(ix)] = perfcurve(Test_top(probes,ix), ar_rec(probes, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, prtb_nn(ix) ] = perfcurve(Test_top(probes,ix), nn_rec(probes, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, prtb_bl(ix) ] = perfcurve(Test_top(probes,ix), bl_rec(probes, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
    [~,~,~, prtb_or(ix) ] = perfcurve(Test_top(probes,ix), or_rec(probes, ix), 1, 'xCrit', 'reca', 'yCrit', 'prec');
%    [~,~,~,~, prtb_rec(ix)] = prec_rec(ar_rec(probes, ix, max_ix), Test_top(probes,ix), 'numThresh',length(Test_top(probes,1))/30);
%    [~,~,~,~, prtb_nn(ix) ] = prec_rec(nn_rec(probes, ix, max_ix), Test_top(probes,ix), 'numThresh',length(Test_top(probes,1))/30);
%    [~,~,~,~, prtb_bl(ix) ] = prec_rec(bl_rec(probes, ix, max_ix), Test_top(probes,ix), 'numThresh',length(Test_top(probes,1))/30);
%    [~,~,~,~, prtb_or(ix) ] = prec_rec(or_rec(probes, ix, max_ix), Test_top(probes,ix), 'numThresh',length(Test_top(probes,1))/30);
end



clrred = [228,26,28]/255;
clrblue = [55,126,184]/255;
clrgreen = [77,175,74]/255;
clrpurple = [152,78,163]/255;
clrblack = [0 0 0];
myBarColors = [clrred; clrblue; clrgreen; clrpurple];

[h_pr_bl_ar p_pr_bl_ar] = kstest2(pr_bl(:), pr_rec(:), 0.05, 'larger');
[h_cc_bl_ar p_cc_bl_ar] = kstest2(cc_bl(:), cc_rec(:), 0.05, 'larger');
[h_pr_nn_ar p_pr_nn_ar] = kstest2(pr_nn(:), pr_rec(:), 0.05, 'larger');
[h_cc_nn_ar p_cc_nn_ar] = kstest2(cc_nn(:), cc_rec(:), 0.05, 'larger');
[~, p_top_nn_ar] = kstest2(Top_nn, Top_rec, 0.05, 'larger');
[~, p_top_bl_ar] = kstest2(Top_bl, Top_rec, 0.05, 'larger');
fid = fopen(sprintf('%s/alleyene_%s_pval.txt', folder, source), 'w');
fprintf(fid, 'performing 2 sided ks test, comparing pr_rec > pr_bl, model\n');
fprintf(fid, 'pvalue = %f\n', p_pr_bl_ar);
fprintf(fid, 'performing 2 sided ks test, comparing pr_rec > pr_nn\n');
fprintf(fid, 'pvalue = %f\n', p_pr_nn_ar);
fprintf(fid, 'performing 2 sided ks test, comparing cc_rec > cc_bl\n');
fprintf(fid, 'pvalue = %f\n', p_cc_bl_ar);
fprintf(fid, 'performing 2 sided ks test, comparing cc_rec > cc_nn\n');
fprintf(fid, 'pvalue = %f\n', p_cc_nn_ar);
fprintf(fid, 'performing 2 sided ks test, comparing top_rec > top_nn\n');
fprintf(fid, 'pvalue = %f\n', p_top_nn_ar);
fprintf(fid, 'performing 2 sided ks test, comparing top_rec > top_bl\n');
fprintf(fid, 'pvalue = %f\n', p_top_bl_ar);
fclose(fid);

h = figure;
plotScatter(pr_bl(:), pr_rec(:), 0.008, 1, [55,126,184]/255, ...
     'Blosum Nearest Neighbor', 'Affinity Regression', 'Test AUPR for top 1% luminescent probes');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/aupr_bl_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(pr_nn(:), pr_rec(:), 0.008, 1, [55,126,184]/255, ...
    'Nearest Neighbor', 'Affinity Regression', 'Test AUPR for top 1% luminescent probes');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/aupr_nn_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
pval = plotScatter(prtb_bl(:), prtb_rec(:), 0.008, 1, [55,126,184]/255, ...
     'Blosum Nearest Neighbor', 'Affinity Regression', 'Test AUPR for top 1% luminescent probes vs. low 50');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/auprtb_bl_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(prtb_nn(:), prtb_rec(:), 0.008, 1, [55,126,184]/255, ...
    'Nearest Neighbor', 'Affinity Regression', 'Test AUPR for top 1% luminescent probes vs. low 50%');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/auprtb_nn_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(cc_bl(:), cc_rec(:), 0.008, 1, [55,126,184]/255, ...
    'Blosum Nearest Neighbor', 'Affinity Regression', 'Test probe correlation');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/cc_bl_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(cc_nn(:), cc_rec(:), 0.008, 1, [55,126,184]/255, ...
    'Nearest Neighbor', 'Affinity Regression', 'Test probe correlation');
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/cc_nn_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(Top_nn(:)/100, Top_rec(:)/100, 0.008, 1, [55,126,184]/255, ...
    'Nearest Neighbor', 'Affinity Regression', sprintf('Prediction rate of top 100 probes p<%f', p_top_nn_ar));
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/top_nn_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');

h = figure;
plotScatter(Top_bl(:)/100, Top_rec(:)/100, 0.008, 1, [55,126,184]/255, ...
    'Blosum Nearest Neighbor', 'Affinity Regression', sprintf('Prediction rate of top 100 probes p<%f', p_top_bl_ar));
axis([-0.01 1.01 -0.01 1.01])
print(h, sprintf('%s/top_bl_vs_ar_%d_%s.pdf', folder, 0, source),'-dpdf');


cc = [cc_or(:) cc_rec(:) cc_bl(:) cc_nn(:)];
pr = [pr_or(:) pr_rec(:) pr_bl(:) pr_nn(:)];
prtb = [prtb_or(:), prtb_rec(:), prtb_bl(:), prtb_nn(:)];
top = [Top_or(:)/100, Top_rec(:)/100, Top_bl(:)/100, Top_nn(:)/100];
pr(isnan(pr)) = 0;
cc_sem = std(cc)/sqrt(size(cc, 1));
pr_sem = std(pr)/sqrt(size(pr, 1));
top_sem = std(top)/sqrt(size(top, 1));
prtb_sem = std(prtb)/sqrt(size(prtb, 1));
S = [cc_sem; prtb_sem; pr_sem;  top_sem];
M = [mean(cc); mean(prtb); mean(pr); mean(top)];
h = figure; 
hbar = barwitherr(S, M, 0.5); 
set(gca,'xticklabel',{'All Probe Correlation', 'Top 1% vs 50% precision', 'Top 1% Precision', 'Top 100 precision'})
axis([0.5 4.5 0.0 min(1, max(M(:))+0.2)])
legend('Oracle Neighbor','Affinity Regression','Blosum Nearest Neighbor','Nearest Neighbor')
ylabel('Performance')

% Set the bar colors
hBarChildren = get(hbar, 'Children');
for i = 1:size(myBarColors, 1)
    set(hBarChildren{i}, 'CDataMapping', 'direct', 'CData', i);
end
colormap(myBarColors);

print(h, sprintf('%s/performance_%d_%s.pdf', folder, 0, source),'-dpdf');
writetable(array2table([M, S], ...
    'rownames',{'All Probe Correlation', 'Top 1% vs 50% precision', 'Top 1% Precision', 'Top 100 precision'},...
    'variablenames',{'oracle_mean','ar_mean', 'blosumnn_mean', 'NearestNeighbor_mean','oracle_stderr','ar_stderr', 'blosumnn_stderr', 'NearestNeighbor_stderr'}),...
    sprintf('%s/performance_%d_%s.csv', folder, 0, source),...
    'writerownames', 1);

if dumppred
    csvwrite(sprintf('%s/pred_ar_rec_%d_%s.csv', folder, 0, source), ar_rec(:,:));
    csvwrite(sprintf('%s/pred_ar_nn_%d_%s.csv', folder, 0, source), nn_rec(:,:));
    csvwrite(sprintf('%s/pred_ar_test_%d_%s.csv', folder, 0, source), Y_gr);
end






