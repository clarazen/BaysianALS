function ttm = outerprodTT(ttA,ttB)

    N = ttA.N;
    GA = ttA.cores;
    GB = ttB.cores;
    rA = ttA.ranks;
    rB = ttB.ranks;
    sA = ttA.Ysize;
    sB = ttB.Ysize;
    
    for i = 1:N
        gA = reshape(GA{i}.data,[numel(GA{i}.data) 1]);
        gB = reshape(GB{i}.data,[1 numel(GB{i}.data)]);
        Gdata = kron(gA,gB);
        Gdata = reshape(Gdata,[rA(i) sA(i) rA(i+1) rB(i) sB(i) rB(i+1)]);
        Gdata = permute(Gdata,[1 4 2 5 3 6]);
        Gdata = reshape(Gdata,[rA(i)*rB(i) sA(i)*sB(i) rA(i+1)*rB(i+1)]);
        G{i}  = Core(Gdata,[rA(i)*rB(i) rA(i+1)*rB(i+1)]);
    end
    ttm = MPT(G);

end