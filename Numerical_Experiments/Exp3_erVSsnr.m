%%% Experiment 3: classical vs. Bayesian ALS - relative error vs. SNR %%%
clear; close all; clc; close all;
%addpath('C:\Program Files\MATLAB\TikZ\src')

preamb

%% User input
sigma   = 0.1;      % standard deviation for prior covariance 
a       = 0.1;      % distance from ground truth for prior mean
K       = 10;       % number of priors
I       = 10;       % number of noisy samples
maxiter = 10;       % maximum iterations in BayALS_MPT
orthog  = false;    % with/without orthogonalization step

%% truth and and experiment parameters
load truth
t               = randi([1 1000]);
expparams.truth = truth{t};
expparams.expnr = 3;
expparams.eps_m = 1e-15;
expparams.eps_P = 1e-15;
D               = 5;

%% Computation of error for different SNRs
SNR = 0:2:25;
count = 1;
for s = 1:numel(SNR)
    [Y,sigmae] = computenoisymeas(truth{t},1,SNR(s));
    Yt         = reshape(Y,[D D D]);
    for k = 1:K
        pmpt0              = computeprior(sigma,a,Truthtt{t}.cores); 
        [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
        erBay(s,k)         = expresults.er;
        [ttCla,expresults] = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,expparams);
        erCla(s,k)         = expresults.er;
        clc; disp(strcat(num2str(ceil(100/(K*numel(SNR))*count)),'% completed'))
        count = count + 1;
    end
end
disp('done')

%% Plot
figure; hold on;
h1 = uncertaintyplot(SNR,erCla','-',yellow,lightyellow,true,2);
h2 = uncertaintyplot(SNR,erBay','-',blue,lightblue,true,2);
xlabel('SNR [dB]')
ylabel('$\varepsilon_\mathrm{truth}$ [-]','interpreter','latex')
xlim([0 25])
legend([h1,h2],'ALS','Algorithm 3.1')
%matlab2tikz('.\figures\fig_exp3_erVSsnr.tex');
