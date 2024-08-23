function A = approxMatrix(ttm)
%   Computes the matrix from a TTm, that represents a square matrix. For
%   all long indices it is I_i = J_i;
%
%   INPUT
%   ttm         TT matrix
%
%   OUTPUT
%   A           matrix
%
%   2020, Clara Menzen
%% 
    newdims = [];
    NSU = ttm.N;
    for i = 1:NSU
        newdims = [newdims sqrt(ttm.Ysize(i)) sqrt(ttm.Ysize(i))];
    end

    newperm = [1:2:NSU*2-1 2:2:NSU*2];

    ttm = approxTensor(ttm).data;
    ttm = reshape(ttm,newdims);
    ttm = permute(ttm,newperm);
    A = reshape(ttm,[sqrt(numel(ttm)) sqrt(numel(ttm))]);
    
end