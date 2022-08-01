%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function plotGrid(w, ch, names, fontsize)
%w is a matrix of weights
%c is a same sized matrix of charachters
%names are the y axis values

assert(all(size(w) == size(ch)))

if nargin < 4
    fontsize = 9;
end



%# text location and labels
[r c] = size(w);
[xloc yloc] = meshgrid(1:c,1:r);
xloc = xloc(:); yloc = yloc(:);
str = cellstr(ch(:));

%# plot colored cells
h = imagesc(1:c, 1:r, w);
set(gca, 'Box','on', 'XAxisLocation','bottom', 'YDir','reverse', ...
    'XLim',[0 c]+0.5, 'YLim',[0 r]+0.5, 'TickLength',[0 0], ...
    'LineWidth',2, 'Color','none', ...
    'FontWeight','bold', 'FontSize', fontsize, 'DataAspectRatio',[1 1 1]);

%# plot grid
xv1 = repmat((2:c)-0.5, [2 1]); xv1(end+1,:) = NaN;
xv2 = repmat([0.5;c+0.5;NaN], [1 r-1]);
yv1 = repmat([0.5;r+0.5;NaN], [1 c-1]);
yv2 = repmat((2:r)-0.5, [2 1]); yv2(end+1,:) = NaN;

%# plot text
text(xloc, yloc, str, 'FontSize',fontsize, 'HorizontalAlignment','center');
set(gca,'YTick',1:size(w,1))
set(gca,'YTicklabel',names)
