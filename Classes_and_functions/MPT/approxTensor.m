function B = approxTensor(mpt)
%%  DESCRIPTION
%   Computes the TT-approximation B to A from an MPT object
%
%   INPUT
%   mpt       MPT object
%   OUTPUT
%   B         Approximated tensor
%    
%   Resources:
%   V. Oseledets: Tensor-Train Decomposition
%
%   July 2020, Clara Menzen
%% 

    B = mpt.cores{1}.data;
    B = reshape(B, [size(B,2) size(B,3)]);
    for i = 2:mpt.N
        Gi = mpt.cores{i};
        B = B * lr_unfolding(Gi,"right");
        B = reshape(B, [numel(B)/Gi.dims(end), Gi.dims(end)]);
    end
    sizesY = [];
    for i = 1:mpt.N
        sizesY = [sizesY mpt.cores{i}.dims(2)];
    end
    B = reshape(B,sizesY);
    B = Tensor(B);
end