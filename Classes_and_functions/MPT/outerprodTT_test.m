close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A = Tensor(1e7*randn(2, 2, 2));
ttA = MPT_SVD(A,[2 2 2],1e-15);
ttB = ttA;
ttm = outerprodTT(ttA,ttB);
M   = approxMatrix(ttm);
A_  = approxTensor(ttA).data;
M_  = reshape(A_,[8 1])*reshape(A_,[1 8]);
norm(M-M_)/norm(M_)