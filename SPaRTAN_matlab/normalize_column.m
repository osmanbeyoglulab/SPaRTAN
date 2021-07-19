function [y] = normalize_column(x)
% 
% Author : Kamlesh Pawar
% Date : April 2014
% Usage : Takes 2D-matrix x and returns y such that each columns of y have unit lenght
% [y] = normalize_column(x)
% Input  : 
%      x : 2D real/complex valued matrix
% Output : 
%      y : Normalized matrix such that each column of y is of unit length
    if (ndims(x)~=2)
        error('Input argument must be a 2D matrix');
    end
    scale = diag(1./sqrt(sum(abs(x).^2)));
    y = x*scale;
end