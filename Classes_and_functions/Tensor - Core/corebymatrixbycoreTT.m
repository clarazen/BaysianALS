function cmc = corebymatrixbycoreTT(core,A,s)
% Computes a core by matrix by core multiplication
% used in computeUTcova

    M  = size(A,1);
    S1 = core.data;
    s1 = core.dims;
    r1 = core.ranks;
    S1 = reshape(S1,[s1(2)/M M r1(2)]);
    S1 = permute(S1,[1 3 2]);
    S1 = reshape(S1,[s1(2)/M*r1(2) M]);
    S1 = S1*A*S1';
    S1 = reshape(S1,[s(2) r1(2) s(2) r1(2)]);
    S1 = permute(S1,[1 3 2 4]);
    cmc = Core(reshape(S1,[1 s(2)^2 r1(2)^2]),[1 r1(2)^2]); 
    
end