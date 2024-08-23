function XtnA = nmodeprod(tensor,A,n)
%%  DESCRIPTION
%   Computes the product between a given tensor and matrix in mode n. The
%   dimension which is summed over is the second one of the matrix A.
%
%   INPUT
%   tensor   Tensor object (I_1 x ... x I_n x ... x I_N)
%   A        matrix (J x I_n)
%   n        mode
%   OUTPUT
%   XtnA     n-mode product (I_1 x ... x J x ... x I_N)
%        ___         ___
%    ___/ A \_______/ X \_________________________
%     J \___/  I_n  \___/ I_1...I_n-1 I_n+1...I_N
%
%              ____
%             /XtnA\
%    =====>   \____/
%             /||||\
%           I_1  J  I_N
%
%   Resources:
%   [1] Kolda, Bader: Tensor Decomposition and Applications, 2009. Page 460
%
%   November 2019, Clara Menzen
%% 
    Xt   = tensor.data;
    dims = tensor.dims;
    dimP = [n 1:n-1 n+1:tensor.order];

    X = permute(Xt,dimP);
    X = reshape(X,[size(Xt,n) numel(X)/size(Xt,n)]);

    XtnA      = A*X; 
    newdim    = dims(dimP);
    newdim(1) = size(A,1);
    XtnA      = reshape(XtnA,newdim);
    XtnA      = permute(XtnA,[2:n 1 n+1:tensor.order]);
    
end