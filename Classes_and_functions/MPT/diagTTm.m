function diagA = diagTTm(ttm)

    for n = 1:ttm.N
        tmp   = permute(ttm.cores{n}.data,[1 3 2]);
        dims  = ttm.cores{n}.dims;
        tmp   = reshape(tmp,[dims(1)*size(ttm.cores{n}.data,3) dims(2)]);
        Cdata = zeros(size(tmp,1), sqrt(size(tmp,2)));
        
        ind = find(reshape(eye(sqrt(ttm.Ysize(n))),[ttm.Ysize(n) 1]));
        for j = 1:sqrt(ttm.Ysize(n))
            Cdata(:,j) = tmp(:,ind(j));
        end
        Cdata = reshape(Cdata,[dims(1) size(ttm.cores{n}.data,3) sqrt(ttm.Ysize(n))]);
        Cdata = permute(Cdata,[1 3 2]);
        G{n}  = Core(Cdata,[ttm.ranks(n) ttm.ranks(n+1)]);
    end
    diagttm = MPT(G);
    diagA   = approxTensor(diagttm).data;
    diagA   = reshape(diagA,[numel(diagA) 1]);
end