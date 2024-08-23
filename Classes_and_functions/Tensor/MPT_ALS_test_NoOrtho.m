clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

r = 3; 
n = 5; 
SNR = 0;
maxiter = 20;

seed1 = 1;
[G,truth]  = computetruth(r,n,seed1);
[y,sigmae] = computenoisymeas(truth,1,SNR);

seed2 = 10;
pmpt0 = computeprior(0,0,G,seed2);

orthog = false;
Yt     = reshape(y,[n n n]);
ttCla  = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,[]);

%% test 
X2   = pmpt0.cores{2}.data;
X3   = pmpt0.cores{3}.data;

for i = 1:maxiter
    %% update first core
    M23 = reshape(reshape(X2,[r*n r])*X3,[r n*n]);
    U1  = kron(M23',eye(n));  
    m1  = U1\y;
    X1  = reshape(m1,[n r]); 
    %% update second core
    U2  = kron(kron(X3',eye(n)),X1); 
    m2  = U2\y;
    X2  = reshape(m2,[r n r]);
    %% update third core
    M12 = reshape(X1*reshape(X2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    m3  = U3\y;
    X3  = reshape(m3,[r n]);
    %% update second core
    U2  = kron(kron(X3',eye(n)),X1); 
    m2  = U2\y;
    X2  = reshape(m2,[r n r]);
    %% reconstructed signal
    yhat  = reshape(reshape(X1*reshape(X2,[r,n*r]),[n*n r])*X3,[n*n*n 1]);
end

%% comparison
norm(yhat-reshape(approxTensor(ttCla).data,size(yhat)))/norm(yhat)