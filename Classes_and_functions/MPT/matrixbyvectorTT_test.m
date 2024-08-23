close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A    = Tensor(randn(64,4));
ttm  = MPT_SVD(A,[4 4; 4 1; 4 1],1e-15);
b    = randn(4,1);

tt = matrixbyvectorTT(ttm,b);
Ab = approxTensor(tt).data;

norm(reshape(Ab,[64 1])-A.data*b)/norm(A.data*b)
