%%% Experiment 3: Classical vs. Bayesian ALS - signal %%%
clear; close all; clc; close all;
%addpath('C:\Program Files\MATLAB\TikZ\src')

preamb

%% User input 
sigma   = 200;      % standard deviation for prior covariance 
a       = 10;       % distance from ground truth for prior mean
SNR     = 0;        % signal-to-noise ratio
I       = 10;       % number of noisy samples
maxiter = 5;        % maximum iterations in BayALS_MPT
orthog  = false;    % with/without orthogonalization step

%% truth and noisy samples y and experiment parameters
load truth
t               = randi([1 1000]);
[Y,sigmae]      = computenoisymeas(truth{t},100,SNR);
pmpt0           = computeprior(sigma,a,Truthtt{t}.cores);
expparams.truth = truth{t};
expparams.expnr = 3.1;
expparams.eps_m = 1e-15;
expparams.eps_P = 1e-15;
D               = 5; 

%% Computation of low-rank tensor estimate with classical ALS and Bayesian ALS
for i = 1:I
    Yt                 = reshape(Y(:,i),[D D D]);
    [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
    ttCla              = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,expparams);
    for n = 1 : pmpt0.N
        m0 = ttBay.cores{n}.mean;
        P0 = ttBay.cores{n}.cova;
        pmpt0.cores{n} = probCore(ttBay.cores{n}.data,pmpt0.cores{n}.ranks,m0,P0);
    end
    pmpt0 = probMPT(pmpt0.cores);
    if orthog == true
        pmpt0 = changeCanonical(pmpt0,frobnorm(pmpt0.cores{1}),1);
    end
    YhatBay(:,i) = reshape(approxTensor(expresults.mUT).data,[numel(truth{t}) 1]);
    PhatBay(:,i) = sqrt(diag(approxMatrix(expresults.PUT)));
    YhatCla(:,i) = reshape(approxTensor(ttCla).data,[numel(truth{t}) 1]);
    clc; disp(strcat(num2str(100/I*i),'% completed'))
end
clc; disp('done')

%% Plots
figure; hold on;
y2 = [1:numel(truth{t}), fliplr(1:numel(truth{t}))];
inBetween = [(YhatBay(:,1)+2*PhatBay(:,1))', fliplr((YhatBay(:,1)-2*PhatBay(:,1))')];
fill(y2, inBetween,[155/255 191/255 188/255],'LineStyle','none');
h1 = plot(truth{t},'Color',red,'Marker','o','LineWidth',1);
h2 = plot(YhatCla(:,1),'Color',yellow,'LineWidth',1);
h3 = plot(YhatBay(:,1),'Color',blue,'LineWidth',1);
legend([h1,h2,h3],'Truth','ALS','BayALS','Location','southwest')
%matlab2tikz('.\figures\fig_exp3_signal_1.tex');

figure; hold on;
y2 = [1:numel(truth{t}), fliplr(1:numel(truth{t}))];
inBetween = [(YhatBay(:,10)+2*PhatBay(:,10))', fliplr((YhatBay(:,10)-2*PhatBay(:,10))')];
fill(y2, inBetween,[155/255 191/255 188/255],'LineStyle','none');
h1 = plot(truth{t},'Color',red,'Marker','o','LineWidth',1);
h2 = plot(YhatCla(:,10),'Color',yellow,'LineWidth',1);
h3 = plot(YhatBay(:,10),'Color',blue,'LineWidth',1);
legend([h1,h2,h3],'Truth','ALS','BayALS','Location','southwest')
%matlab2tikz('.\figures\fig_exp3_signal_100.tex');
