clear; close all; clc; format compact;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

r = 3; 
n = 5; 
SNR = 0;
maxiter = 10;

seed1        = randi([1,100]);
[G,Truthtt]  = computetruth(r,n,seed1);
Truth        = approxTensor(Truthtt);
truth        = reshape(Truth.data, [numel(Truth.data),1]);
[y,~]        = computenoisymeas(truth,1,SNR);

pmpt0        = Truthtt;

%% test
G1   = reshape(pmpt0.cores{1}.data,[n r]);
G2   = pmpt0.cores{2}.data;
G3   = pmpt0.cores{3}.data;
y0n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
norm(y0n)
%% Bring prior into site-k-mixed canonical form with norm in first core
[Q,R] = qr(G3',0);
Q3_0  = Q';
%
X2_0  = reshape(G2,[r*n r])*R';
%
X2_0  = reshape(X2_0,[r n*r]);
[Q,R] = qr(X2_0',0);
Q2_0  = reshape(Q',[r n r]);
X1_0 = G1*R';

% initial configuration
X1 = X1_0; Q2 = Q2_0; Q3 = Q3_0; 
y0o = reshape(reshape(X1*reshape(Q2,[r,n*r]),[n*n r])*Q3,[n*n*n 1]);
norm(y0n)-norm(y0o);

for i = 1:maxiter
    %% update first core
    M23 = reshape(reshape(Q2,[r*n r])*Q3,[r n*n]);
    U1  = kron(M23',eye(n));  
    m1  = U1'*y;
    norm1o = norm(m1)
    X1  = reshape(m1,[n r]); 
    % orthogonalize first core
    [Q1,~] = qr(X1,0);
    %% update second core
    U2  = kron(kron(Q3',eye(n)),Q1); 
    m2  = U2'*y;
    norm2o = norm(m2)
    G2  = reshape(m2,[r n r]);
    % orthogonalize second core
    [Q,~] = qr(reshape(G2,[r*n r]),0);
    Q2 = reshape(Q,[r n r]);
    %% update third core
    M12 = reshape(Q1*reshape(Q2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    m3  = U3'*y;
    norm3o = norm(m3)
    G3  = reshape(m3,[r n]); 
    % orthogonalize third core
    [Q,~] = qr(G3',0);
    Q3 = reshape(Q',[r n]);
    %% update second core
    U2  = kron(kron(Q3',eye(n)),Q1); 
    m2  = U2'*y;
    norm2o = norm(m2)
    G2  = reshape(m2,[r n r]);
    % orthogonalize second core
    [Q,R] = qr((reshape(G2,[r n*r]))',0);
    Q2 = reshape(Q',[r n r]);
    % move norm to first core
    X1 = Q1*R';
    yhat_o  = reshape(reshape(X1*reshape(Q2,[r,n*r]),[n*n r])*Q3,[n*n*n 1]);
    er_o = norm(truth-yhat_o)/norm(truth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% update first core
    M23 = reshape(reshape(G2,[r*n r])*G3,[r n*n]);
    U1  = kron(M23',eye(n));  
    m1  = U1\y;
    G1  = reshape(m1,[n r]); 
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    norm1n = norm(yhat_n)
    %% update second core
    U2  = kron(kron(G3',eye(n)),X1); 
    m2  = U2\y;
    G2  = reshape(m2,[r n r]);
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    norm2n = norm(yhat_n)
    %% update third core
    M12 = reshape(G1*reshape(G2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    m3  = U3\y;
    G3  = reshape(m3,[r n]);
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    norm3n = norm(yhat_n)
    %% update second core
    U2  = kron(kron(G3',eye(n)),G1); 
    m2  = U2\y;
    G2  = reshape(m2,[r n r]);
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    norm2n = norm(yhat_n)
    %% reconstructed signal
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    er_n = norm(truth-yhat_n)/norm(truth)
    %% Comparison
    d_no = norm(yhat_o-yhat_n)
end
er_o-er_n
