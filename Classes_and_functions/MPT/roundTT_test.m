close all; clear; clc;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

n1 = 10; 
n2 = 5;
n3 = 5;

G1{1} = Core(randn(1,n1,3),[1 3]);
G1{2} = Core(randn(3,n2,3),[3 3]);
G1{3} = Core(randn(3,n3,1),[3 1]);
tt1 = MPT(G1);

G2{1} = Core(randn(1,n1,3),[1 3]);
G2{2} = Core(randn(3,n2,3),[3 3]);
G2{3} = Core(randn(3,n3,1),[3 1]);
tt2 = MPT(G2);

tts  = addTTs({tt1,tt2}); 
ttsr = roundTT(tts,1e-15);
frobnorm(approxTensor(tts))-frobnorm(approxTensor(ttsr))

