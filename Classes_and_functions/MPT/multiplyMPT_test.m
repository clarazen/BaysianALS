close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

At   = Tensor(randn(5,5,5));
a    = reshape(At.data,[125 1]);
ttA  = MPT_SVD(At,[5 5 5],1e-15);
b    = 5;
mpt  = multiplyMPT(ttA,b);

aA   = approxTensor(mpt).data;

norm(reshape(aA,[125 1])-a*b)/norm(a*b)