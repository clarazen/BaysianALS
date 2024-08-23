%%% Experiment on covariance evolution %%%
clear; close all; clc; close all; format compact;
%addpath('C:\Program Files\MATLAB\TikZ\src')

preamb

%% User input 
sigma   = 200;      % standard deviation for prior covariance 
a       = 10;       % distance from ground truth for prior mean
SNR     = 0;        % signal-to-noise ratio
K       = 5;        % number of priors
maxiter = 10;       % maximum iterations in BayALS_MPT
orthog  = true;     % with/without orthogonalization step

%% truth and noisy sample y and experiment parameters
load truth
t               = randi([1 1000]);
[y,sigmae]      = computenoisymeas(truth{t},1,SNR);
D               = 5; 
Yt              = reshape(y,[D D D]);
expparams.expnr = 2;
expparams.truth = truth{t};
expparams.eps_m = 1e-15;
expparams.eps_P = 1e-15;

%% Computation of covariance matrices' Frobenius norm and trace
for k = 1:K
    pmpt0              = computeprior(sigma,a,Truthtt{t}.cores);
    [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
    normPn{k}          = expresults.normPn;
    tracePn{k}         = expresults.tracePn;
    normP(k,:)       = expresults.normP;
    traceP(k,:)      = expresults.traceP;
    clc; disp(strcat(num2str(ceil(100/K*k)),'% completed'))
end
clc; disp('done')

%% Plots
for n = 1:pmpt0.N
    for k = 1:K
        TracePn(:,k) = tracePn{k}(:,n);
    end
    figure; hold on;
    x = 1:maxiter;
    uncertaintyplot(x,TracePn(x,:)','-',blue,lightblue,false,2);
    xlabel('Number of iterations','interpreter','latex')
    ylab = strcat('$\mathrm{tr}(\mathbf{P}_',num2str(n),')$');
    ylabel(ylab,'interpreter','latex')
    xlim([1 maxiter])
    ylim([TracePn(maxiter,1)-500 TracePn(maxiter,1)+1000])
cleanfigure;
%matlab2tikz(strcat('.\figures\fig_exp2_trP_',num2str(n),'.tex'),'width', '\fwidth')
end

for n = 1:pmpt0.N
    for k = 1:K
        NormPn(:,k) = normPn{k}(:,n);
    end
    figure; hold on;
    x = 1:maxiter;
    uncertaintyplot(x,NormPn(x,:)','-',blue,lightblue,false,2);
    xlabel('Number of iterations','interpreter','latex')
    ylab = strcat('$||\mathbf{P}_',num2str(n),'||_\mathrm{F}$');
    ylabel(ylab,'interpreter','latex')
    xlim([1 maxiter])
    ylim([NormPn(maxiter,1)-500 NormPn(maxiter,1)+1000])
cleanfigure;
%matlab2tikz(strcat('.\figures\fig_exp2_normP_',num2str(n),'.tex'),'width', '\fwidth')
end

figure(3); hold on;
set(gca,'XTickLabelMode','manual')
subplot(1,2,1)
uncertaintyplot(x,traceP,'-',blue,lightblue,true,2);
xlabel('Number of iterations','interpreter','latex')
ylabel('$\mathrm{tr}(\mathbf{P}_\mathrm{UT})$','interpreter','latex')
xlim([1 5])

subplot(1,2,2)
uncertaintyplot(x,normP,'-',blue,lightblue,true,2);
xlabel('Number of iterations','interpreter','latex')
ylabel('$||\mathbf{P}_\mathrm{UT}||_\mathrm{F}$','interpreter','latex')
xlim([1 5])
%matlab2tikz('.\figures\fig_exp2_norm_traceP.tex');

