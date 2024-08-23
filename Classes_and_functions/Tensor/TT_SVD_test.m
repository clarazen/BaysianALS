clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

Yt = randn(5,5,5);

tt = TT_SVD(Tensor(Yt),1e-15);

frobnorm(Tensor(approxTensor(tt).data-Yt))/frobnorm(Tensor(Yt))