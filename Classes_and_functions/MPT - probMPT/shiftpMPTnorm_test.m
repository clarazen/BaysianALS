close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

r = 3; n = 5;
X1   = randn(1,n,r);
X2   = randn(r,n,r);
X3   = randn(r,n,1);
x1   = reshape(X1,[n,r]);

P1   = eye(numel(X1));
P2   = eye(numel(X2));
P3   = eye(numel(X3));
m1   = reshape(X1,[numel(X1) 1]);
m2   = reshape(X2,[numel(X2) 1]);
m3   = reshape(X3,[numel(X3) 1]);
C{1} = probCore(X1,[1,r],m1,P1);
C{2} = probCore(X2,[r,r],m2,P2);
C{3} = probCore(X3,[r,1],m3,P3);
tt   = probMPT(C);
y    = reshape(approxTensor(tt).data,[125 1]);
tt   = shiftpMPTnorm(tt,1,1);
tt   = shiftpMPTnorm(tt,2,1);
tt   = changeCanonical(tt,frobnorm(tt.cores{3}),3);



