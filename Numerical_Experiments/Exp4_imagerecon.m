%%% Image reconstruction %%%
clear; clc; close all; format compact; format long; 

preamb

image = double(rgb2gray(imread('./figures/Cat.png')));
subplot(1,3,1); imshow(uint8(image))
truth = reshape(image,[numel(image) 1]);
Ysize = [4 4 4 4 4 4 4 4];
Truth = Tensor(reshape(truth,Ysize));
tt = MPT_SVD(Truth,Ysize,0.1);
Ylr = reshape(approxTensor(tt).data,[numel(truth) 1]);
subplot(1,3,2); imshow(uint8(reshape(Ylr,size(image))))

expparams.truth = truth;
expparams.expnr = 4;
orthog          = false;
%% Vectorization of corrupted tensor 
SNR = 0;
Y = [];
for i = 1:10
    noise = randn(size(truth));
    noise_norm = norm(truth)/(10^(SNR/20));
    noise = noise/norm(noise)*noise_norm;
    y = Ylr+noise;
    Y = [Y y];
end
subplot(1,3,3); imshow(uint8(reshape(y,size(image))))

%% Data for ALS
maxiter = 3;    
sigmae  = noise_norm^2/(numel(noise)-1);

%% prior 
N = 8;
for j = 1 : N
    m0 = randn(tt.ranks(j)*tt.Ysize(j)*tt.ranks(j+1),1);
    P0 = 1000^2*eye(tt.ranks(j)*tt.Ysize(j)*tt.ranks(j+1));
    G0  = reshape(m0,[tt.ranks(j),tt.Ysize(j),tt.ranks(j+1)]);
    prior.cores{j} = probCore(G0,[tt.ranks(j) tt.ranks(j+1)],m0,P0);
end
pmpt0 = probMPT(prior.cores);

%% ALS 
figure(2); hold on;
for i = 1:size(Y,2)
    tic
    Yt          = reshape(Y(:,i),[4 4 4 4 4 4 4 4]);
    ttCla       = MPT_ALS(Tensor(Yt),pmpt0,maxiter,orthog,[]);
    YhatCla{i}  = reshape(approxTensor(ttCla).data,[numel(truth) 1]);
    pmpt0       = ttCla;
    erCla(i)    = norm(YhatCla{i}-truth)/norm(truth);
    toc
end
subplot(2,2,1)
imshow(uint8(reshape(YhatCla{1},size(image))))
subplot(2,2,2)
imshow(uint8(reshape(YhatCla{10},size(image))))
%% Bayesian ALS
for i = 1:size(Y,2)
    tic
    Yt = reshape(Y(:,i),[4 4 4 4 4 4 4 4]);
    ttBay = MPT_BayALS(Tensor(Yt),pmpt0,sigmae,maxiter,orthog,[]);
    for j = 1 : pmpt0.N
        m0 = ttBay.cores{j}.mean;
        P0 = ttBay.cores{j}.cova;
        pmpt0.cores{j} = probCore(ttBay.cores{j}.data,pmpt0.cores{j}.ranks,m0,P0);
    end
    pmpt0 = probMPT(pmpt0.cores);
    if orthog == true
        pmpt0 = changeCanonical(pmpt0,frobnorm(pmpt0.cores{1}),1);
    end
    toc
    YhatBay{i} = reshape(approxTensor(ttBay).data,[numel(truth) 1]);
    erBay(i)   = norm(YhatBay{i}-truth)/norm(truth);
    toc
end
subplot(2,2,3)
imshow(uint8(reshape(YhatBay{1},size(image))))
subplot(2,2,4)
imshow(uint8(reshape(YhatBay{10},size(image))))

