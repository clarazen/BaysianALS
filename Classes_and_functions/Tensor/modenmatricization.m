function X = modenmatricization(tensor,n)
%%  DESCRIPTION
%   Computes the mode-n matricization matrix for a given tensor
%
%   INPUT
%   tensor   object of Tensor class (I_1 ... I_N)
%   n        mode
%   OUTPUT
%   Xn       mode-n matricization matrix (I_n x I_1...I_n-1 I_n+1...I_N)
%
%   Resources:
%   [1] Kolda, Bader: Tensor Decomposition and Applications, 2009. 
%       Page 460
%
%   November 2019, Clara Menzen
%%   
    Xt = tensor.data;
    D  = tensor.dims;
    N  = tensor.order;

    Xt  = permute(Xt,[n 1:n-1 n+1:N]);
    X   = reshape(Xt,[D(n) numel(Xt)/D(n)]);
end