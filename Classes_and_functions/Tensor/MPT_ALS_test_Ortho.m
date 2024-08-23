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
maxiter = 100;

seed1 = 1;
[G,truth]  = computetruth(r,n,seed1);
[y,sigmae] = computenoisymeas(truth,1,SNR);

seed2 = 10;
pmpt0 = computeprior(0,0,G,seed2);

orthog = true;
Yt     = reshape(y,[n n n]);
ttCla  = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,[]);

%% test
X2_0   = pmpt0.cores{2}.data;
X3_0   = pmpt0.cores{3}.data;

%% Bring prior into site-k-mixed canonical form with norm in first core
[Q,R] = qr(X3_0',0);
Q3_0  = Q';
%
X2_0  = reshape(X2_0,[r*n r])*R';
X2_0  = reshape(X2_0,[r n r]);
%
X2_0  = reshape(X2_0,[r n*r]);
[Q,~] = qr(X2_0',0);
Q2_0  = reshape(Q',[r n r]);

% initial configuration
Q2 = Q2_0; Q3 = Q3_0; 

for i = 1:maxiter
    %% update first core
    M23 = reshape(reshape(Q2,[r*n r])*Q3,[r n*n]);
    U1  = kron(M23',eye(n));  
    m1  = U1'*y;
    X1  = reshape(m1,[n r]); 
    % orthogonalize first core
    [Q1,~] = qr(X1,0);
    %% update second core
    U2  = kron(kron(Q3',eye(n)),Q1); 
    m2  = U2'*y;
    X2  = reshape(m2,[r n r]);
    % orthogonalize second core
    [Q,~] = qr(reshape(X2,[r*n r]),0);
    Q2 = reshape(Q,[r n r]);
    %% update third core
    M12 = reshape(Q1*reshape(Q2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    m3  = U3'*y;
    X3  = reshape(m3,[r n]); 
    % orthogonalize third core
    [Q,~] = qr(X3',0);
    Q3 = reshape(Q',[r n]);
    %% update second core
    U2  = kron(kron(Q3',eye(n)),Q1); 
    m2  = U2'*y;
    X2  = reshape(m2,[r n r]);
    % orthogonalize second core
    [Q,R] = qr((reshape(X2,[r n*r]))',0);
    Q2 = reshape(Q',[r n r]);
    % move norm to first core
    X1 = Q1*R';
    %% Reconstructed signal
    yhat  = reshape(reshape(X1*reshape(Q2,[r,n*r]),[n*n r])*Q3,[n*n*n 1]);
end

%% Comparison
norm(yhat-reshape(approxTensor(ttCla).data,size(yhat)))/norm(yhat)