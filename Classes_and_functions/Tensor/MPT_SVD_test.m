clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

A = randn(24,24);

mpo = MPT_SVD(Tensor(A),[2 2; 3 3; 4 4],1e-15);

norm(approxMatrix(mpo)-A)/norm(A)