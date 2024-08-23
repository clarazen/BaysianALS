function [mpt,frn] = roundTT(mpt, epsilon)
%%  DESCRIBTION
%   Computes the rounding of a TT, meaning that it reduces the ranks 
%   First, the cores are orthogonalized from right to left and then the
%   compression with a truncated SVD is performed to reduce the ranks
%
%   INPUT
%   G        tensor in TT-format -> cores
%   epsilon  required accuracy
%   desranks desired ranks
%   OUTPUT
%   B        tensor in TT-format with TT-ranks rk_hat equal to the
%            delta-ranks of the unfoldings Ak of A where 
%            delta = epsilon/sqrt(d-1)||A||_F. The computed approximation 
%            satisfies ||A-B||_F <= epsilon||A||_F
%
%   Resources:
%   [1] I. V. Oseledets: Tensor-Train Decomposition, 2011, p.2305: Algorithm 2
%
%   December 2019, Clara Menzen
%%   
    G  = mpt.cores;
    N  = mpt.N;
    n  = linspace(N,2,N-1);
    % Right-to-left orthogonalization
    for i = 1:numel(n)
        k      = n(i);
        Gtmp   = lr_unfolding(G{k},'right');
        [Q,R]  = qr(Gtmp',0);
        Gdata  = reshape(Q', [ size(Q,2) mpt.Ysize(k) numel(Q)/(size(Q,2)*mpt.Ysize(k))]);
        G{k}   = changeCdata(G{k},Gdata,[size(Gdata,1) size(Gdata,3)]);
        mpt.cores{k} = G{k};
        Gdata  = nmodeprod(G{k-1},R,3);
        G{k-1} = changeCdata(G{k-1},Gdata,[size(Gdata,1) size(Gdata,3)]);
        mpt.cores{k-1} = G{k-1};
    end
    % Compression of the orthogonalized representation
    for k = 1:N-1
        [U,S,V] = svd(lr_unfolding(G{k},'left'),'econ');
        sv = diag(S);
        delta = epsilon/sqrt(N-1)*sqrt(innerprodTT(mpt,mpt));
        St = S; j = 0;
        frn = sv(end);
        while frn <= delta && j<size(S,1)-1
            St  = St(1:end-1,1:end-1);
            j   = j+1;
            frn = norm(sv(end-j:end));
        end
        Ut = U(:,1:size(U,2)-j);
        Vt = V(:,1:size(V,2)-j);

        Gdata  = reshape(Ut, [numel(Ut)/(mpt.Ysize(k)*size(Vt,2)) mpt.Ysize(k) size(Vt,2)] );
        G{k}   = changeCdata(G{k},Gdata,[size(Gdata,1) size(Gdata,3)]);
        Gdata  = nmodeprod(G{k+1},(Vt*St)',1);
        G{k+1} = changeCdata(G{k+1},Gdata,[size(Gdata,1) size(Gdata,3)]);

        mpt.cores{k}   = G{k};
        mpt.cores{k+1} = G{k+1};
        mpt = changeMPT(mpt,G);
    end
end