clear; close all; clc; format compact;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

%% Truth: low rank tensor train contracted to the true tensor
r = 5; n = 10;
rng(1)
X1 = randn(1,n,r);
rng(5)
X2 = randn(r,n,r);
rng(7)
X3 = randn(r,n,1);
x1 = reshape(X1,[n,r]);

P1 = eye(numel(X1));
P2 = eye(numel(X2));
P3 = eye(numel(X3));
m1 = reshape(X1,[numel(X1) 1]);
m2 = reshape(X2,[numel(X2) 1]);
m3 = reshape(X3,[numel(X3) 1]);
C{1} = probCore(X1,[1,r],m1,P1);
C{2} = probCore(X2,[r,r],m2,P2);
C{3} = probCore(X3,[r,1],m3,P3);

mpt       = probMPT(C);

epsm      = 1e-20;
epsP      = 1e-20;
[mUT,PUT] = UT_pmpt(mpt,epsm,epsP,false);
mUT_ten   = approxTensor(mUT).data;
mUT_vec   = reshape(mUT_ten,[numel(mUT_ten) 1]);
%sd        = sqrt(diag(approxMatrix(PUT)));

X1 = reshape(mpt.cores{1}.data,[n r]);
X2 = mpt.cores{2}.data;
X3 = mpt.cores{3}.data;
P1 = mpt.cores{1}.cova;
P2 = mpt.cores{2}.cova;
P3 = mpt.cores{3}.cova;
[mUT_,PUT_]= UT(X1,X2,X3,P1,P2,P3,r,n);
sd_ = sqrt(diag(PUT_));

norm(mUT_-mUT_vec)/norm(mUT_)
%norm(sd_-sd)/norm(sd_)


function [mUT,PUT] = UT(X1,X2,X3,P1,P2,P3,r,n)

    n1 = numel(X1);
    n2 = numel(X2);
    n3 = numel(X3);
    M = n1 + n2 + n3; 
    
    % parameters
    kappa  = 3-M; beta = 2; alpha = .001;
    lambda = alpha^2*(M+kappa) - M;
    
    % joint mean
    m  = [reshape(X1,[n1 1]); reshape(X2,[n2 1]); reshape(X3,[n3 1])];
    
    % joint covariance
    cholP                                    = zeros(M);  
    cholP(1:n1,1:n1)                         = chol(P1)';
    cholP(n1+1:n1+n2,n1+1:n1+n2)             = chol(P2)';
    cholP(n1+n2+1:n1+n2+n3,n1+n2+1:n1+n2+n3) = chol(P3)';
    
    % weights
    w = 1/(2*(M+lambda));
    wm0 = lambda/(M+lambda);
    wc0 = lambda/(M+lambda) + (1-alpha^2+beta);
    
    
    %% 1.) Computation of sigma points and 2.) Propagation thru non linear function
    for i = 1:M
        sp = m + sqrt(M+lambda)*cholP(:,i);
        c21 = reshape(sp(1:n1),[n r]);
        c22 = reshape(sp(n1+1:n1+n2),[r n r]);
        c23 = reshape(sp(n1+n2+1:n1+n2+n3),[r n]);
        Stp(:,i) = reshape(reshape(c21*reshape(c22,[r,n*r]),[n*n r])*c23,[n*n*n 1]);
        
        sm = m - sqrt(M+lambda)*cholP(:,i);
        c31 = reshape(sm(1:n1),[n r]);
        c32 = reshape(sm(n1+1:n1+n2),[r n r]);
        c33 = reshape(sm(n1+n2+1:n1+n2+n3),[r n]);
        Stm(:,i) = reshape(reshape(c31*reshape(c32,[r,n*r]),[n*n r])*c33,[n*n*n 1]);
    end
    
    St1 = reshape(reshape(X1*reshape(X2,[r,n*r]),[n*n r])*X3,[n*n*n 1]);
    St  = [Stp, Stm];
    
    %% 3.) Computation of mean and covariance 
    mUT = wm0*St1;
    for i = 1:2*M
        mUT = mUT + w*St(:,i);
    end
    
    PUT = wc0*(St1-mUT)*(St1-mUT)';
    for i = 1:2*M
        D = (St(:,i) - mUT);
        PUT = PUT + w*(D*D');
    end
    
end

