function tt = matrixbyvectorTT(A,b)
    %% DESCRIPTION   
    % Computes the product of a matrix in TTm format with a vector by
    % multiplying the first core with the vector. Used in UT_pmpt.
    %
    % INPUT
    % A     TTm with first core that has a long index
    % b     Vector of size of the one of the long indices of A
    % OUTPUT
    %
    % November 2020, Clara Menzen
%%
    M     = size(b,1);
    G     = A.cores;
    s     = size(A.cores{1}.data);
    C1    = reshape(A.cores{1}.data,[s(1) s(2)/M M s(3)]);
    C1    = permute(C1,[1 2 4 3]);
    C1    = reshape(C1, [s(1)*s(2)/M*s(3) M] );
    C1    = C1*b;
    C1    = reshape(C1,[s(1) s(2)/M s(3)]); 
    G{1}  = changeCdata(A.cores{1},C1,[size(C1,1) size(C1,3)]);
    tt    = MPT(G);
    
end