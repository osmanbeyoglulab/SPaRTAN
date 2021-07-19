%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       converts a trained model from ar_train to W
%

function w = ar_model2w(model)
w = model.Va * pinv(model.Sa) * reshape(model.fit_svd.beta(:,1),size(model.Va,2),size(model.Ub,2)) * pinv(model.Sb) * model.Ub';
