%%% Convergence study of likelihood and relative error %%%
clear; close all; clc; close all; format compact
%addpath('C:\Program Files\MATLAB\TikZ\src')

preamb

%% User input
sigma   = 200;      % standard deviation for prior covariance 
a       = 10;       % distance from ground truth for prior mean
SNR     = 0;        % signal-to-noise ratio
K       = 10;       % number of priors
maxiter = 5;        % maximum iterations in BayALS_MPT
orthog  = false;    % with/without orthogonalization step

%% truth and noisy sample y and experiment parameters
load truth
t               = randi([1 1000]);
[y,sigmae]      = computenoisymeas(truth{t},1,SNR);
D               = 5; 
Yt              = reshape(y,[D D D]);
expparams.expnr = 1;
expparams.truth = truth{t};
expparams.eps_m = 1e-15;
expparams.eps_P = 1e-15; 

%% Computation of errors and log likelihood times prior
for k = 1:K 
    pmpt0              = computeprior(sigma,a,Truthtt{t}.cores);
    [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
    er_truth(k,:)      = expresults.er_truth;
    er_meas(k,:)       = expresults.er_meas;
    Likprior(k,:)      = expresults.Likprior;
    clc; disp(strcat(num2str(ceil(100/K*k)),'% completed'))
end
clc; disp('done')

%% plots
figure; hold on;
subplot(1,2,1)
x = 1:maxiter;
htruth = uncertaintyplot(x,er_truth(:,x),'-',blue,lightblue,true,2);
hmeas  = uncertaintyplot(x,er_meas(:,x),'--',blue,lightgreen,true,2);

xlabel('Number of iterations','interpreter','latex')
ylabel('Relative error','interpreter','latex')
legend([htruth,hmeas],'$\varepsilon_\mathrm{truth}$',...
       '$\varepsilon_\mathrm{meas}$',...
       'interpreter','latex')
xlim([1 maxiter])
ylim([0.3 1.2])

subplot(1,2,2)
x = 1:maxiter;
uncertaintyplot(x,Likprior(:,x),'-',blue,lightblue,false,2);

xlabel('Number of iterations','interpreter','latex')
ylabel('$\log \mathrm{Likelihood } \cdot \mathrm{ Prior}$','interpreter','latex')
xlim([1 maxiter])
%matlab2tikz('.\figures\fig_exp1.tex');

