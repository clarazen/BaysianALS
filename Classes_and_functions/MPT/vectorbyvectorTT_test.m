close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A   = Tensor(randn(5,5,5));
ttA = MPT_SVD(A,[5 5 5],1e-15);
a   = reshape(A.data,[125 1]);

b   = randn(10,1);

ttm = vectorbyvectorTT(ttA,b);
tt  = matrixbyvectorTT(ttm,b);

norm(a*b'*b-reshape(approxTensor(tt).data,[125 1]))/norm(a*b'*b)