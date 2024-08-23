close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

Y  = Tensor(randn(5,5,5));
y  = reshape(Y.data,[125 1]);

r = 3; n = 5;
P1 = eye(numel(X1));
P2 = eye(numel(X2));
P3 = eye(numel(X3));
m1 = reshape(X1,[numel(X1) 1]);
m2 = reshape(X2,[numel(X2) 1]);
m3 = reshape(X3,[numel(X3) 1]);
C{1} = probCore(X1,[1,r],m1,P1);
C{2} = probCore(X2,[r,r],m2,P2);
C{3} = probCore(X3,[r,1],m3,P3);

tt   = probMPT(C);

tt.canonical.norm - norm(y)
Y1 = approxTensor(tt).data;
tt = shiftMPTnorm(tt,3,-1);
tt.canonical.norm - norm(y)
Y2 = approxTensor(tt).data;
tt = shiftMPTnorm(tt,2,-1);
tt.canonical.norm - norm(y)
Y3 = approxTensor(tt).data;
tt = shiftMPTnorm(tt,1,1);
tt.canonical.norm - norm(y)
Y4 = approxTensor(tt).data;
tt = shiftMPTnorm(tt,2,1);
tt.canonical.norm - norm(y)
Y5 = approxTensor(tt).data;

norm(reshape(Y1,[125 1])-y)/norm(y)
norm(reshape(Y2,[125 1])-y)/norm(y)
norm(reshape(Y3,[125 1])-y)/norm(y)
norm(reshape(Y4,[125 1])-y)/norm(y)
norm(reshape(Y5,[125 1])-y)/norm(y)