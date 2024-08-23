clear; close all; clc;

addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

r = 3; 
n = 5; 
SNR = 200;
maxiter = 1;

seed1 = 1;
[G,truth]  = computetruth(r,n,seed1);
[y,sigmae] = computenoisymeas(truth,1,SNR);

seed2 = 10;
pmpt0 = computeprior(10,10,G,seed2);

%% test
G1_0 = pmpt0.cores{1}.data;
G2_0 = pmpt0.cores{2}.data;
G3_0 = pmpt0.cores{3}.data;
P1_0 = pmpt0.cores{1}.cova;
P2_0 = pmpt0.cores{2}.cova;
P3_0 = pmpt0.cores{3}.cova;

%% Bring prior into site-1-mixed-canonical form
G1_0 = reshape(G1_0,[n r]);
[Q,R] = qr(G3_0',0);
Q3_0  = Q';
P3o_0  = kron(eye(n),pinv(R)') *P3_0* kron(eye(n),pinv(R));
%
X2_0  = reshape(G2_0,[r*n r])*R';
X2_0  = reshape(X2_0,[r n r]);
P2o_0  = kron(R, eye(r*n)) *P2_0* kron(R',eye(r*n)); 
%
X2_0  = reshape(X2_0,[r n*r]);
[Q,R] = qr(X2_0',0);
Q2_0  = reshape(Q',[r n r]);
P2o_0  = kron(eye(r*n),pinv(R)') *P2o_0* kron(eye(r*n),pinv(R));
%
X1_0  = G1_0*R';
P1o_0  = kron( R,eye(n) ) *P1_0* kron( R',eye(n) );

X1 = X1_0; Q2 = Q2_0; Q3 = Q3_0; 
P1o = P1o_0;
P2o = P2o_0;
P3o = P3o_0;

% prior mean vectors
m1o_0 = reshape(X1_0,[numel(X1_0) 1]);
m2o_0 = reshape(Q2_0,[numel(Q2_0) 1]);
m3o_0 = reshape(Q3_0,[numel(Q3_0) 1]);


G1 = G1_0; G2 = G2_0; G3 = G3_0;
m1_0 = reshape(G1_0,[numel(G1_0) 1]);
m2_0 = reshape(G2_0,[numel(G2_0) 1]);
m3_0 = reshape(G3_0,[numel(G3_0) 1]);


for i = 1:maxiter
    %% update first core
    M23      = reshape(reshape(Q2,[r*n,r])*Q3,[r n*n]);
    U1       = kron(M23',eye(n));  
    P1o      = pinv( pinv(P1o_0) + eye(size(P1o_0))/sigmae^2 );
    m1       = P1o * ((U1'*y)/sigmae^2 + pinv(P1_0)*m1_0);
    X1       = reshape(m1,[n r]); 
    % orthogonalize first core
    [Q1,R]   = qr(X1,0);
    P1o      = kron( pinv(R)',eye(n) ) *P1o* kron( pinv(R),eye(n) );
    % orthogonalize first core prior
    [Q1_0,R] = qr(G1_0,0);
    P1o_0     = kron( pinv(R)',eye(n) ) *P1o_0* kron( pinv(R),eye(n) );
    % move norm to second core prior
    G2_0     = R*reshape(Q2_0,[r n*r]);
    m2_0     = reshape(G2_0,[numel(G2_0) 1]);
    P2o_0     = kron( eye(r*n),R ) *P2o_0* kron( eye(r*n),R' );
    %% update second core
    U2       = kron(kron(Q3',eye(n)),Q1); 
    P2       = pinv( pinv(P2o_0) + eye(size(P2o_0))/sigmae^2 );
    m2       = P2 * ((U2'*y)/sigmae^2 + pinv(P2o_0)*m2_0);
    X2       = reshape(m2,[r n r]); 
    % orthogonalize second core
    [Q,R]    = qr(reshape(X2,[r*n r]),0);
    Q2       = reshape(Q,[r n r]);
    P2o       = kron( pinv(R)',eye(r*n) ) *P2o* kron( pinv(R),eye(r*n) );
    % orthogonalize second core prior
    [Q,R]    = qr(reshape(G2_0,[r*n r]),0);
    Q2_0     = reshape(Q,[r n r]);
    P2o_0     = kron( pinv(R)',eye(r*n) ) *P2o_0* kron( pinv(R),eye(r*n) );
    % move norm to third core prior
    G3_0     = R*Q3_0;
    m3_0     = reshape(G3_0,[numel(G3_0) 1]);
    P3o_0    = kron( eye(n),R) *P3o_0* kron(eye(n),R' );
    %% update third core
    M12      = reshape(Q1*reshape(Q2,[r n*r]),[n*n r]);
    U3       = kron(eye(n),M12);
    P3o      = pinv( pinv(P3o_0) + eye(size(P3o_0))/sigmae^2 );
    m3       = P3o * ((U3'*y)/sigmae^2 + pinv(P3o_0)*m3_0);
    X3       = reshape(m3,[r n]); 
    % orthogonalize third core
    [Q,R]    = qr(X3',0);
    Q3       = reshape(Q',[r n]);
    P3o       = kron( eye(n),pinv(R)') *P3o* kron(eye(n),pinv(R) );
    % orthogonalize prior of third core
    G3_0     = reshape(m3_0,[r n]);
    [Q,R]    = qr(G3_0',0);
    Q3_0     = reshape(Q',[r n]);
    P3o_0     = kron( eye(n),pinv(R)') *P3o_0* kron(eye(n),pinv(R) );
    % move norm to prior of second core
    G2_0     = reshape(Q2_0,[r*n r])*R';
    m2_0     = reshape(G2_0,[numel(G2_0) 1]);
    P2o_0     = kron( R, eye(r*n) ) *P2o_0* kron( R',eye(r*n) );      
    %% update second core
    U2       = kron(kron(Q3',eye(n)),Q1); 
    P2o       = pinv( pinv(P2o_0) + eye(size(P2o_0))/sigmae^2 );
    m2       = P2o * ((U2'*y)/sigmae^2 + pinv(P2o_0)*m2_0);
    X2       = reshape(m2,[r n r]); 
    % orthogonalize second core
    [Q,R]    = qr((reshape(X2,[r n*r]))',0);
    Q2       = reshape(Q',[r n r]);
    P2o       = kron( eye(r*n),pinv(R)' ) *P2o* kron( eye(r*n),pinv(R) );
    % move norm to first core
    X1       = Q1*R';
    m1       = reshape(X1,[numel(X1) 1]);
    P1o      = kron( R,eye(n) ) *P1o* kron( R',eye(n) );
    % orthogonalize prior of second core
    [Q,R]    = qr((reshape(G2_0,[r n*r]))',0);
    Q2_0     = reshape(Q',[r n r]);
    P2o_0     = kron( eye(r*n),pinv(R)' ) *P2o_0* kron( eye(r*n),pinv(R) );
    % move norm to prior of first core
    X1_0     = Q1_0*R';
    m1_0     = reshape(G1_0,[numel(G1_0) 1]);
    P1o_0    = kron( R,eye(n) ) *P1o_0* kron( R',eye(n) );   
    yhat_o  = reshape(reshape(X1*reshape(Q2,[r,n*r]),[n*n r])*Q3,[n*n*n 1]);
    norm(truth-yhat_o)/norm(truth)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% update first core
    M23 = reshape(reshape(G2,[r*n,r])*G3,[r n*n]);
    U1  = kron(M23',eye(n));  
    P1  = inv( inv(P1_0) + (U1'*U1)/sigmae^2 );
    m1  = P1 * ((U1'*y)/sigmae^2 + inv(P1_0)*m1_0);  
    G1  = reshape(m1,[n r]); 
    %% update second core
    U2  = kron(kron(G3',eye(n)),G1); 
    P2  = inv( inv(P2_0) + (U2'*U2)/sigmae^2 );
    m2  = P2 * ((U2'*y)/sigmae^2 + inv(P2_0)*m2_0);
    G2  = reshape(m2,[r n r]);
    %% update third core
    M12 = reshape(X1*reshape(G2,[r n*r]),[n*n r]);
    U3  = kron(eye(n),M12);
    P3  = inv( inv(P3_0) + (U3'*U3)/sigmae^2 );
    m3  = P3 * ((U3'*y)/sigmae^2 + inv(P3_0)*m3_0);
    G3  = reshape(m3,[r n]); 
    %% update second core
    U2  = kron(kron(G3',eye(n)),G1); 
    P2  = inv( inv(P2_0) + (U2'*U2)/sigmae^2 );
    m2  = P2 * ((U2'*y)/sigmae^2 + inv(P2_0)*m2_0);
    G2  = reshape(m2,[r n r]);
    %% reconstructed signal
    yhat_n  = reshape(reshape(G1*reshape(G2,[r,n*r]),[n*n r])*G3,[n*n*n 1]);
    norm(truth-yhat_n)/norm(truth)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Comparison
    norm(yhat_o-yhat_n)
end 