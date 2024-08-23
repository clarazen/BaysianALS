clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A     = Tensor(randn(64,64));
ttm   = MPT_SVD(A,[4 4; 4 4; 4 4],1e-15);
norm(approxMatrix(ttm)-A.data,'fro')/norm(A.data,'fro')
