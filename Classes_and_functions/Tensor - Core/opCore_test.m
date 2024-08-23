clear; close all; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')
addpath('..\..\Numerical_Experiments\functions')

Yt = randn(5,5,5);
tt = TT_SVD(Tensor(Yt),1e-15);

OP = 1;
for n = 1:tt.N
    G{n} = opCore(tt.cores{n});
end

test1 = reshape(approxTensor(MPT(G)).data,[25^3 1]);

test2 = reshape(approxTensor(outerprodTT(tt,tt)).data,[25^3 1]);

norm(test1-test2)/norm(test2)