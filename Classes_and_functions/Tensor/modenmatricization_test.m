clear; close all; clc;

%   Resources:
%   Kolda, Bader: Tensor Decomposition and Applications, 2009., Page 460

Xt(:,:,1) = [1 4 7 10; 2 5 8 11; 3 6 9 12];
Xt(:,:,2) = [1 4 7 10; 2 5 8 11; 3 6 9 12]+12;
tensor = Tensor(Xt);

X1 = modenmatricization(tensor,1);
X2 = modenmatricization(tensor,2);
X3 = modenmatricization(tensor,3);