function norm = frobnorm(tensor)
%%  DESCRIBTION
%   Computes the Frobenius norm of a Tensor object
%   INPUT
%   tensor  Tensor object
%   OUTPUT
%   norm    Frobenius norm of Tensor object
%
%   December 2019, Clara Menzen
%%
    y     = reshape(tensor.data, [numel(tensor.data) 1]);
    norm  = sqrt(y'*y);
end
        