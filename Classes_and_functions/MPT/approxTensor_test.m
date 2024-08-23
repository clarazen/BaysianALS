clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A   = Tensor(randn(5,5,5));
a   = reshape(A.data,[125 1]);
tt  = MPT_SVD(A,[5,5,5],1e-15);
D   = approxTensor(tt).data-A.data;
d   = reshape(D,[125 1]);
norm(d)/norm(a)
