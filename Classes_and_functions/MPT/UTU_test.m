close all; clear; clc; format compact;

addpath('..\MPT')
addpath('..\MPT - probMPT')
addpath('..\Tensor')
addpath('..\Tensor - Core')
addpath('..\Tensor - Core - probCore')

Y  = Tensor(randn(5,5,5));
tt = MPT_SVD(Y,[5 5 5],1e-15);
    
[G_left, G_right] = supercore(tt,1);
U = kron(kron(G_right',eye(tt.Ysize(1))),G_left); 
norm(UTU(tt,1)-U'*U)/norm(U'*U)

[G_left, G_right] = supercore(tt,2);
U = kron(kron(G_right',eye(tt.Ysize(2))),G_left); 
norm(UTU(tt,2)-U'*U)/norm(U'*U)

[G_left, G_right] = supercore(tt,3);
U = kron(kron(G_right',eye(tt.Ysize(3))),G_left); 
norm(UTU(tt,3)-U'*U)/norm(U'*U)

function [G_left, G_right] = supercore(tt,n)
    G = tt.cores;
    r = tt.ranks;
    if n == 1
        G_left = 1;
    else
        G_left = lr_unfolding(G{1},"left");
        for i = 2:n-1
            G_left = G_left * lr_unfolding(G{i},"right");
            G_left = reshape(G_left,[r(i+1) numel(G_left)/r(i+1)]);
        end
        G_left = reshape(G_left,[numel(G_left)/r(n) r(n)]);
    end

    if n == numel(G)
        G_right = 1;
    else
        G_right = lr_unfolding(G{n+1},"left");
        for i = n+2 : numel(G)
            G_right = G_right * lr_unfolding(G{i},"right");
            G_right = reshape(G_right,[r(i+1) numel(G_right)/r(i+1)]);
        end
        G_right = reshape(G_right,[r(n+1) numel(G_right)/r(n+1)]);
    end 
end
