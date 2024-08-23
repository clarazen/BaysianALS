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
pmpt0 = computeprior(10,10,G,seed2);

orthog    = false;
Yt        = reshape(y,[n n n]);
ttBay     = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,[]);
[mUT,PUT] = UT_pmpt(ttBay,1e-15,1e-15,true);

%% test
X1_0 = pmpt0.cores{1}.data;
X2_0 = pmpt0.cores{2}.data;
X3_0 = pmpt0.cores{3}.data;
P1_0 = pmpt0.cores{1}.cova;
P2_0 = pmpt0.cores{2}.cova;
P3_0 = pmpt0.cores{3}.cova;
    
% prior mean vectors
m1_0 = reshape(X1_0,[numel(X1_0) 1]);
m2_0 = reshape(X2_0,[numel(X2_0) 1]);
m3_0 = reshape(X3_0,[numel(X3_0) 1]);

% initial configuration
X2 = X2_0;
X3 = X3_0;

for iter = 1: maxiter
    %% update first core
    M23 = reshape(reshape(X2,[r*n,r])*X3,[r n*n]);
    U1  = kron(M23',eye(n));  
    P1  = inv( inv(P1_0) + (U1'*U1)/sigmae^2 );
    m1  = P1 * ((U1'*y)/sigmae^2 + inv(P1_0)*m1_0);  
    X1  = reshape(m1,[n r]); 
    %% update second core
    U2  = kron(kron(X3',eye(n)),X1); 
    P2  = inv( inv(P2_0) + (U2'*U2)/sigmae^2 );
    m2  = P2 * ((U2'*y)/sigmae^2 + inv(P2_0)*m2_0);
    X2  = reshape(m2,[r n r]);
    %% update third core
    M12 = reshape(X1*reshape(X2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    P3  = inv( inv(P3_0) + (U3'*U3)/sigmae^2 );
    m3  = P3 * ((U3'*y)/sigmae^2 + inv(P3_0)*m3_0);
    X3  = reshape(m3,[r n]); 
end
%% Reconstructed signal
[yhat,covhat] = UT(X1,X2,X3,P1,P2,P3,r,n);

%% comparison
norm(yhat-reshape(approxTensor(mUT).data,size(yhat)))/norm(yhat)
norm(approxMatrix(PUT)-covhat)/norm(covhat)
