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

orthog    = true;
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

        
%% Bring prior into site-1-mixed-canonical form
X1_0 = reshape(X1_0,[n r]);
[Q,R] = qr(X3_0',0);
Q3_0  = Q';
P3_0  = kron(eye(n),pinv(R)') *P3_0* kron(eye(n),pinv(R));
%
X2_0  = reshape(X2_0,[r*n r])*R';
X2_0  = reshape(X2_0,[r n r]);
P2_0  = kron(R, eye(r*n)) *P2_0* kron(R',eye(r*n)); 
%
X2_0  = reshape(X2_0,[r n*r]);
[Q,R] = qr(X2_0',0);
Q2_0  = reshape(Q',[r n r]);
P2_0  = kron(eye(r*n),pinv(R)') *P2_0* kron(eye(r*n),pinv(R));
%
X1_0  = X1_0*R';
P1_0  = kron( R,eye(n) ) *P1_0* kron( R',eye(n) );

% initial configuration
m1_0 = reshape(X1_0, [numel(X1_0) 1]);
Q2 = Q2_0; Q3 = Q3_0; 

for i = 1:maxiter
    %% update first core
    M23      = reshape(reshape(Q2,[r*n,r])*Q3,[r n*n]);
    U1       = kron(M23',eye(n));  
    P1       = pinv( pinv(P1_0) + eye(size(P1_0))/sigmae^2 );
    m1       = P1 * ((U1'*y)/sigmae^2 + pinv(P1_0)*m1_0);
    X1       = reshape(m1,[n r]); 
    % orthogonalize first core
    [Q1,R]   = qr(X1,0);
    P1       = kron( pinv(R)',eye(n) ) *P1* kron( pinv(R),eye(n) );
    % orthogonalize first core prior
    [Q1_0,R] = qr(X1_0,0);
    P1_0     = kron( pinv(R)',eye(n) ) *P1_0* kron( pinv(R),eye(n) );
    % move norm to second core prior
    X2_0     = R*reshape(Q2_0,[r n*r]);
    m2_0     = reshape(X2_0,[numel(X2_0) 1]);
    P2_0     = kron( eye(r*n),R ) *P2_0* kron( eye(r*n),R' );
    %% update second core
    U2       = kron(kron(Q3',eye(n)),Q1); 
    P2       = pinv( pinv(P2_0) + eye(size(P2_0))/sigmae^2 );
    m2       = P2 * ((U2'*y)/sigmae^2 + pinv(P2_0)*m2_0);
    X2       = reshape(m2,[r n r]); 
    % orthogonalize second core
    [Q,R]    = qr(reshape(X2,[r*n r]),0);
    Q2       = reshape(Q,[r n r]);
    P2       = kron( pinv(R)',eye(r*n) ) *P2* kron( pinv(R),eye(r*n) );
    % orthogonalize second core prior
    [Q,R]    = qr(reshape(X2_0,[r*n r]),0);
    Q2_0     = reshape(Q,[r n r]);
    P2_0     = kron( pinv(R)',eye(r*n) ) *P2_0* kron( pinv(R),eye(r*n) );
    % move norm to third core prior
    X3_0     = R*Q3_0;
    m3_0     = reshape(X3_0,[numel(X3_0) 1]);
    P3_0     = kron( eye(n),R) *P3_0* kron(eye(n),R' );
    %% update third core
    M12      = reshape(Q1*reshape(Q2,[r n*r]),[n*n r]);
    U3       = kron(eye(n),M12);
    P3       = pinv( pinv(P3_0) + eye(size(P3_0))/sigmae^2 );
    m3       = P3 * ((U3'*y)/sigmae^2 + pinv(P3_0)*m3_0);
    X3       = reshape(m3,[r n]); 
    % orthogonalize third core
    [Q,R]    = qr(X3',0);
    Q3       = reshape(Q',[r n]);
    P3       = kron( eye(n),pinv(R)') *P3* kron(eye(n),pinv(R) );
    % orthogonalize prior of third core
    X3_0     = reshape(m3_0,[r n]);
    [Q,R]    = qr(X3_0',0);
    Q3_0     = reshape(Q',[r n]);
    P3_0     = kron( eye(n),pinv(R)') *P3_0* kron(eye(n),pinv(R) );
    % move norm to prior of second core
    X2_0     = reshape(Q2_0,[r*n r])*R';
    m2_0     = reshape(X2_0,[numel(X2_0) 1]);
    P2_0     = kron( R, eye(r*n) ) *P2_0* kron( R',eye(r*n) );      
    %% update second core
    U2       = kron(kron(Q3',eye(n)),Q1); 
    P2       = pinv( pinv(P2_0) + eye(size(P2_0))/sigmae^2 );
    m2       = P2 * ((U2'*y)/sigmae^2 + pinv(P2_0)*m2_0);
    X2       = reshape(m2,[r n r]); 
    % orthogonalize second core
    [Q,R]    = qr((reshape(X2,[r n*r]))',0);
    Q2       = reshape(Q',[r n r]);
    P2       = kron( eye(r*n),pinv(R)' ) *P2* kron( eye(r*n),pinv(R) );
    % move norm to first core
    X1       = Q1*R';
    m1       = reshape(X1,[numel(X1) 1]);
    P1       = kron( R,eye(n) ) *P1* kron( R',eye(n) );
    % orthogonalize prior of second core
    [Q,R]    = qr((reshape(X2_0,[r n*r]))',0);
    Q2_0     = reshape(Q',[r n r]);
    P2_0     = kron( eye(r*n),pinv(R)' ) *P2_0* kron( eye(r*n),pinv(R) );
    % move norm to prior of first core
    X1_0     = Q1_0*R';
    m1_0     = reshape(X1_0,[numel(X1_0) 1]);
    P1_0     = kron( R,eye(n) ) *P1_0* kron( R',eye(n) );    
end 
%% Reconstructed signal
   [yhat,covhat] = UT(X1,Q2,Q3,P1,P2,P3,r,n);
    norm(yhat-reshape(approxTensor(mUT).data,size(yhat)))/norm(yhat)
    norm(approxMatrix(PUT)-covhat)/norm(covhat)