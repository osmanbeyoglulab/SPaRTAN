function pvals = plotScatter(x, y, size, trans, color, xlabel, ylabel, title, test, dashed)
%x is base, y is alternative
if nargin < 9
    test = 1;
    dashed = 1;
end

hold on
fontName = 'Helvetica';
fontSizeTitle = 16;
fontSizeLabel = 14;

transparentScatter(x, y, size, trans, color);
if dashed
    dashline([0 1],[0 1], 2,3,2,3,'color',[0.5 0.5 0.5],'linewidth', 2);
end
axis equal
axis([0 1 0 1]);
set(gca, 'xtick',linspace(0,1,5));
set(gca, 'ytick',linspace(0,1,5));
set(gca,'box','on');
set(gca,'fontsize',10,'fontname','fontName');
set(get(gca,'Title'), 'string', title,'fontsize', fontSizeTitle,'fontname', fontName);
set(get(gca,'xlabel'), 'string', xlabel,'fontsize', fontSizeLabel,'fontname', fontName);
set(get(gca,'ylabel'), 'string', ylabel,'fontsize', fontSizeLabel,'fontname', fontName);


if test
    [h pvals.ks_pval_over] = kstest2(x, y, 0.05, 'larger');
    [h pvals.ks_pval_under] = kstest2(y, x, 0.05, 'larger');
    [h pvals.sr_pval_over] = signrank(x, y, 0.05, 'tail', 'right');
    [h pvals.sr_pval_under] = signrank(y, x, 0.05, 'tail', 'right');
    text(0.05, 0.9, sprintf('%s<%s, p<%.4e,ks', xlabel, ylabel, pvals.ks_pval_over))
    text(0.05, 0.85, sprintf('%s>%s, p<%.4e,ks', xlabel, ylabel, pvals.ks_pval_under))
    text(0.05, 0.8, sprintf('%s<%s, p<%.4ef,sr', xlabel, ylabel, pvals.sr_pval_over))
    text(0.05, 0.75, sprintf('%s>%s, p<%.4ef,sr', xlabel, ylabel, pvals.sr_pval_under))
end

text(0.05, 0.7, sprintf('R^2 = %.3f', corr(x, y)^2))
