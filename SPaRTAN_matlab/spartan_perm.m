%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Hatice Osmanbeyoglu, UPITT, 2019
%fix the random seed for figure 
clear;
clc;
close;  
addpath(genpath('/Users/osmanbeyogluhu2/OneDrive - University of Pittsburgh/MATLAB/')); % load c-files
addpath(genpath('/MALSAR1.1/')); % load c-files
addpath(genpath('/SLEP_package_4.1/'));
addpath(genpath('/Users/osmanbeyogluhu2/OneDrive - University of Pittsburgh/MATLAB/MALSAR1.1/affinity/pmpack_gamma_1/'));
 
stream = RandStream.getGlobalStream;
reset(stream,100);
%%%load and setup data
%load('/Users/osmanbeyogluhu2/OneDrive - University of Pittsburgh/MATLAB/MALSAR1.1/affinity/10X_pbmc5k_102820.mat');
load('/Users/osmanbeyogluhu2/OneDrive - University of Pittsburgh/MATLAB/MALSAR1.1/affinity/10X_pbmc5k_cibersort_5kb_gene.mat');
 
task=1;
D = double(D);
%D = D(randperm(size(D,1)),:);
p = cell(task ,1);
y = cell(task, 1);
%%
% p{1} = normalize_column(Ppbmc5kB')';
% p{2} = normalize_column(Ppbmc5kcd4')';
% p{3} = normalize_column(Ppbmc5kcd8')';
% p{4} = normalize_column(Ppbmc5kmye')';
% p{5} = normalize_column(Ppbmc5kNK')';
% p{6} = normalize_column(Ppbmc5kgdelta')';
% 
% y{1} = normalize_column(Ypbmc5kB);
% y{2} = normalize_column(Ypbmc5kcd4);
% y{3} = normalize_column(Ypbmc5kcd8);
% y{4} = normalize_column(Ypbmc5kmye);
% y{5} = normalize_column(Ypbmc5kNK);
% y{6} = normalize_column(Ypbmc5kgdelta);
 
% %% center across genes and proteins
p{1} = normalize_column(center(Ppbmc5kB)')';
p{2} = normalize_column(center(Ppbmc5kcd4)')';
p{3} = normalize_column(center(Ppbmc5kcd8)')';
p{4} = normalize_column(center(Ppbmc5kmye)')';
p{5} = normalize_column(center(Ppbmc5kNK)')';
p{6} = normalize_column(center(Ppbmc5kgdelta)')';
 
y{1} = normalize_column(center(Ypbmc5kB')');
y{2} = normalize_column(center(Ypbmc5kcd4')');
y{3} = normalize_column(center(Ypbmc5kcd8')');
y{4} = normalize_column(center(Ypbmc5kmye')');
y{5} = normalize_column(center(Ypbmc5kNK')');
y{6} = normalize_column(center(Ypbmc5kgdelta')');
%% center each cell
% p{1} = normalize_column(center(Ppbmc5kB'))';
% p{2} = normalize_column(center(Ppbmc5kcd4'))';
% p{3} = normalize_column(center(Ppbmc5kcd8'))';
% p{4} = normalize_column(center(Ppbmc5kmye'))';
% p{5} = normalize_column(center(Ppbmc5kNK'))';
% p{6} = normalize_column(center(Ppbmc5kgdelta'))';
% 
% y{1} = normalize_column(center(Ypbmc5kB));
% y{2} = normalize_column(center(Ypbmc5kcd4));
% y{3} = normalize_column(center(Ypbmc5kcd8));
% y{4} = normalize_column(center(Ypbmc5kmye));
% y{5} = normalize_column(center(Ypbmc5kNK));
% y{6} = normalize_column(center(Ypbmc5kgdelta));
%%
D_norm = bsxfun(@times, double(D), sqrt(1./sum(D.^2)));
clear Ypbmc5k Ypbmc5kB Ypbmc5kcd4 Ypbmc5kcd8 Ypbmc5kgdelta Ypbmc5kmye Ypbmc5kNK ...
    Ppbmc5k Ppbmc5kB Ppbmc5kcd4 Ppbmc5kcd8 Ppbmc5kgdelta Ppbmc5kmye Ppbmc5kNK Dpbmc5k
 
clear Ypbmc5kn Ypbmc5knB Ypbmc5kncd4 Ypbmc5kncd8 Ypbmc5kngdelta Ypbmc5knmye Ypbmc5knNK ...
    Ppbmc5kn Ppbmc5knB Ppbmc5kncd4 Ppbmc5kncd8 Ppbmc5kngdelta Ppbmc5kmye Ppbmc5knNK Dpbmc5kn
 
clear Ypbmc1k Ypbmc1kB Ypbmc1kcd4 Ypbmc1kcd8 Ypbmc1kgdelta Ypbmc1kmye Ypbmc1kNK ...
    Ppbmc1k Ppbmc1kB Ppbmc1kcd4 Ppbmc1kcd8 Ppbmc1kgdelta Ppbmc1kmye Ppbmc1kNK Dpbmc1k
%%
 
% optimal learning setting
lambdast = 0;
rsL2st = 0.001;
spectrumA = 1;
spectrumB = 0.65;
 
% affinity regression train and predict
Y_pred = cell(task, 1);
W_st = cell(task ,1);
projP = cell(task, 1);
projD = cell(task, 1);
model = cell(task, 1);
for t = 1:task
    model{t} = ar_train(D_norm,p{t},y{t},lambdast, rsL2st, spectrumA, spectrumB);
    Y_pred{t} = ar_predict(D_norm,p{t},y{t}, model{t});
    W_st{t} = ar_model2w(model{t});
    projP{t} = (y{t}'*D_norm*W_st{t})';
    projD{t} = W_st{t}*p{t}';
end
filename='/Users/osmanbeyogluhu2/OneDrive - University of Pittsburgh/MATLAB/MALSAR1.1/affinity/output/pbmc/pbmc5kmodel/pbmc5k-STL-lambda-0-rsl2-0.001-65-B-cibersort-5kb-gene-center-by-feature.mat';
save(filename,'projP','projD','W_st', 'p','D_norm','y', 'model');
%%
N = 1000;
TFa_perm = cell(N ,1);
Pa_perm = cell(N ,1);
for n = 1:N
    rng('shuffle');
    y_perm = y{1}(randperm(size(y{1},1)),randperm(size(y{1},2)));
    model = ar_train(D_norm, p{1},y_perm,lambdast, rsL2st, spectrumA, spectrumB);
    new_pred = ar_predict(D_norm, p{1}, y{1}, model);
    w_perm = ar_model2w(model);
    TFa_perm{n} = w_perm*p{1}';
    Pa_perm{n} = y{1}'*D_norm*w_perm;
    n
end
%%
% # Determine p-value
% p = sum(abs(Tcalc) >= abs(Tobs)) / I
 
my_tfa = zeros(size(D_norm,2),size(p{1},1));
for t = 1:N
    for tf = 1:size(D_norm,2)
        for patient = 1:size(p{1},1)
            if(abs(TFa_perm{t}(tf,patient)) >= abs(projD{1}(tf,patient)))
                    my_tfa(tf,patient) = my_tfa(tf,patient) + 1;
            end
        end
    end
end
fname = sprintf('my_tfa_pbmc5k_cibersort_5kb_gene_1000_B_cf.csv', datestr(now, 30));
xlswrite(fname, (my_tfa));
 
% my_tfa_left = zeros(size(D_norm,2),size(p{1},1));
% my_tfa_right = zeros(size(D_norm,2),size(p{1},1));
% for t = 1:N
%     for tf = 1:size(D_norm,2)
%         for patient = 1:size(p{1},1)
%             if(TFa_perm{t}(tf,patient)>projD{1}(tf,patient))
%                 my_tfa_left(tf,patient) = my_tfa_left(tf,patient) + 1;
%             end;
%             if(TFa_perm{t}(tf,patient)<projD{1}(tf,patient))
%                 my_tfa_right(tf,patient) = my_tfa_right(tf,patient) + 1;
%             end;
%         end
%     end
% end
% 
% fname_left = sprintf('my_tfa_left_pbmc5k_cibersort_5kb_gene_10000_B_cf.csv', datestr(now, 30));
% fname_right = sprintf('my_tfa_right_pbmc5k_cibersort_5kb_gene_10000_B_cf.csv', datestr(now, 30));
% xlswrite(fname_left, (my_tfa_left+1)/(N+1));
% xlswrite(fname_right, (my_tfa_right+1)/(N+1));
%
% my_pa_left = zeros(size(p{1},2),size(p{1},1));
% my_pa_right = zeros(size(p{1},2),size(p{1},1));
% for n = 1:N
%     for tf = 1:size(p{1},2)
%         for patient = 1:size(p{1},1)
%             if(Pa_perm{n}(patient,tf)>projP{1}(patient,tf))
%                 my_pa_left(tf,patient) = my_pa_left(tf,patient) + 1;
%             end;
%             if(Pa_perm{n}(patient,tf)<projP{1}(patient,tf))
%                 my_pa_right(tf,patient) = my_pa_right(tf,patient) + 1;
%             end;
%         end
%     end
% end
% 
% fname_left = sprintf('my_pa_left_pbmc5k_cibersort_2kb_gene_1000Bwc%s.csv', datestr(now, 30));
% fname_right = sprintf('my_pa_right_pbmc5k_cibersort_2kb_gene_1000Bwc%s.csv', datestr(now, 30));
% xlswrite(fname_left, (my_pa_left+1)/N);
% xlswrite(fname_right,( my_pa_right+1)/N);
%%%
