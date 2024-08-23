clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A     = randn(64,64);
diagA = diag(A);
A     = Tensor(A);
ttm   = MPT_SVD(A,[4 4; 4 4; 4 4],1e-15);
test  = diagTTm(ttm);

max(abs(diagA-test))