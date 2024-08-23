function ttm = vectorbyvectorTT(y,b)
    %% DESCRIPTION   
    % Computes the product of a vector in TTm format with a vector by
    % multiplying the first core with the vector with a rank-1 connection.
    % Used in UT_pmpt.
    %
    % INPUT
    % y     TT 
    % b     Vector
    % OUTPUT
    %
    % November 2020, Clara Menzen
    %%
    G    = y.cores;
    M    = size(b,1);
    s    = y.cores{1}.dims;
    C1   = y.cores{1}.data;
    C1   = reshape(C1,[numel(C1) 1]);
    C1   = C1*b';
    C1   = reshape(C1,[s(1) s(2) s(3) M]);
    C1   = permute(C1,[1 2 4 3]);
    C1   = reshape(C1,[s(1) s(2)*M s(3)]);
    G{1} = changeCdata(G{1},C1,[size(C1,1) size(C1,3)]);
    ttm  = MPT(G);

end