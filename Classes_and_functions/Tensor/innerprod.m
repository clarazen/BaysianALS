function p = innerprod(tensor1,tensor2)
%%  DESCRIPTION   
%   Computes the inner product for two same sized tensors. Output is
%   scalar.
%
%   INPUT
%   tensor1    tensor 1
%   tensor2    tensor 2
%
%   OUTPUT
%   p          inner product of tensor 1 and 2
%
%                 ___        
%               /     \ 
%              |       |            ___
%               \ ___ /           /     \ 
%                 | |      ->    |       |
%                 |_|             \ ___ /
%               /     \
%              |       |
%               \ ___ /  
%
%   Resources:
%   Kolda, Bader: Tensor Decomposition and Applications, 2009. Page 458
%
%   November 2019, Clara Menzen
%%  
        Xt = tensor1.data;
        Yt = tensor2.data;
        y1 = reshape(Xt, [1 numel(tensor1.data)]);
        y2 = reshape(Yt, [1 numel(tensor2.data)]);
        p  = y1*y2';
end 