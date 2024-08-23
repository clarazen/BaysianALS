%%% Experiment 3: different priors %%%
clear; close all; clc; close all; format compact;
%addpath('C:\Program Files\MATLAB\TikZ\src')

preamb

%% user input
SNR     = 0;        % signal-to-noise ratio
K       = 5;      % number of priors
maxiter = 5;        % maximum iterations in BayALS_MPT
orthog  = false;    % with/without orthogonalization step


%% truth and noisy sample y and experiment parameters
load truth
t                =604; %= randi([1 1000]);
[Y,sigmae]       = computenoisymeas(truth{t},1,SNR);
D                = 5; 
Yt               = reshape(Y,[D D D]);
expparams.expnr  = 3;
expparams.truth  = truth{t};
expparams.eps_m  = 1e-15;
expparams.eps_P  = 1e-15;

%% Computation of error for different priors
a = [.01 .03 .07 .1 .2 .3 .5 1 2 3 5];
b = a;
count = 1;
for k = 1:K
    for i = 1:numel(b)
        for j = 1:numel(a)
            pmpt0              = computeprior(b(i),a(j),Truthtt{t}.cores);
            [ttBay,expresults] = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,expparams);
            erBay(i,j,k)       = expresults.er;
            [ttCla,expresults] = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,expparams);
            erCla(i,j,k)       = expresults.er;
            clc; disp(strcat(num2str(ceil(100/(K*numel(a)^2)*count)),'% completed'))
            count = count + 1;
        end
    end
end
clc; disp('done')

%% Plots
[X,Y]      = meshgrid(a,b);
mean_erBay = mean(erBay,3);
for i = 1:numel(a)
    for j = 1:numel(b)
        if mean_erBay(i,j) > 1
            mean_erBay(i,j) = 1;
        end
    end
end
figure; hold on;
mean_erCla = mean(mean(mean(erCla,3)));
contour(X,Y,mean_erBay',[ 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9],'LineWidth',2)
xlabel('$b$','interpreter','latex')
ylabel('$a$','interpreter','latex')
C = contour(X,Y,mean_erBay',[mean_erCla mean_erCla],'LineWidth',1,'Linecolor','black','Linestyle','--');
cb = colorbar();
colormap(turbo(10))
set(cb, 'Ticks', [.1,.2,.3,.4,.5,.6,.7,.8,.9,1],'TickLabels',{'10','20','30','40','50','60','70','80','90','>100%'},'FontName','Times','FontSize',28)
caxis([0 1])
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times';
cb.Label.String = 'Relative Error';

