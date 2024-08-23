clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

%   TEST:
%   A+A-A 

A    = Tensor(randn(3,4,3));
ttA  = MPT_SVD(A,[3 4 3],1e-15);
mttA = multiplyMPT(ttA,-1);
tts  = {ttA,ttA,mttA};
test = addTTs(tts);
Y    = approxTensor(test);
frobnorm(Tensor(A.data-Y.data))