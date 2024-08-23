close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

A = Tensor(randn(5,5,5));
a = reshape(A.data,[125 1]);
ttA = MPT_SVD(A,[5 5 5],1e-15);

B = Tensor(randn(5,5,5));
b = reshape(B.data,[125 1]);
ttB = MPT_SVD(B,[5 5 5],1e-15);

dot(a,b)-innerprodTT(ttA,ttB)