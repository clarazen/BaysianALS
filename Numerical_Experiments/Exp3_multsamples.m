%%% Experiment 3: Classical vs. Bayesian ALS - rel. error vs. number of noisy samples %%%
clear; close all; clc; close all;
%addpath('C:\Program Files\MATLAB\TikZ\src')
addpath('.\2include')

preamb

%% User input
sigma   = 200;      % standard deviation for prior covariance 
a       = 10;       % distance from ground truth for prior mean
SNR     = 0;        % signal-to-noise ratio
K       = 10;       % number of priors
I       = 10;       % number of noisy samples
maxiter = 5;        % maximum iterations in BayALS_MPT
orthog  = false;    % with/without orthogonalization step

%% truth and noisy sample y and experiment parameters
load truth
t               = randi([1 1000]);
[Y,sigmae]      = computenoisymeas(truth{t},I,SNR);
expparams.expnr = 3;
expparams.truth = truth{t};
expparams.eps_m = 1e-15;
expparams.eps_P = 1e-15;
D               = 5;

%% Computation of error for increasing number of noisy samples
count = 1;
for k = 1:K
    pmpt0 = computeprior(sigma,a,Truthtt{t}.cores);
    ttCla = pmpt0;
    for i = 1:I  
        Yt                 = reshape(Y(:,i),[D D D]);
        [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
        erBay(k,i)         = expresults.er(end);
        for n = 1 : pmpt0.N
            m0 = ttBay.cores{n}.mean;
            P0 = ttBay.cores{n}.cova;
            pmpt0.cores{n} = probCore(ttBay.cores{n}.data,pmpt0.cores{n}.ranks,m0,P0);
        end
    pmpt0 = probMPT(pmpt0.cores);
    if orthog == true
        pmpt0 = changeCanonical(pmpt0,frobnorm(pmpt0.cores{1}),1);
    end
    [ttCla,expresults] = MPT_ALS(Tensor(Yt),ttCla,maxiter,orthog,expparams);
    erCla(k,i) = expresults.er;
    clc; disp(strcat(num2str(ceil(100/(K*I)*count)),'% completed'))
    count = count +1;
    end
end
disp('done')

%% Plot
figure; hold on;
h1 = uncertaintyplot(1:I,erCla,'-',yellow,lightyellow,false,2);
h2 = uncertaintyplot(1:I,erBay,'-',blue,lightblue,false,2);
xlabel('Number of noisy samples [-]')
ylabel('$\varepsilon_\mathrm{truth}$ [-]','interpreter','latex')
legend([h1,h2],'ALS','Algorithm 3.1')
%matlab2tikz('.\figures\fig_exp3_multsamples.tex');