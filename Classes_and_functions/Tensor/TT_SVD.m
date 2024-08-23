function tt = TT_SVD(tensor,eps)
%   Computes the cores of a TT for the given tensor and accuracy

%   INPUT
%   tensor   Tensor object
%   esp      accuracy
%   OUTPUT
%   tt       MPT object corrisponding to a TT

%   Resources:
%   V. Oseledets: Tensor-Train Decomposition, 2011, p.2301: Algorithm 1

%   November 2019, Clara Menzen
%%
    N = tensor.order;             
    D = tensor.dims;         
    % Initialization: truncation parameter
    delta = eps / sqrt(N-1) * frobnorm(tensor);
    % Temporary tensor and r0 = rprev = 1
    C = tensor.data; rprev = 1;
    for k = 1 : N-1
        n = D(k);
        C = reshape( C, [ rprev*n, numel(C) / (rprev* n) ]);
        % truncated svd 
        [U,S,V] = svd(C,'econ');
        St = S; j = 0;
        sv = diag(S);
        frn = sv(end);
        while frn <= delta && j<size(S,1)-1
            St  = St(1:end-1,1:end-1);
            j  = j+1;
            frn = norm(sv(end-j:end));
        end
        Ut = U(:,1:size(U,2)-j);
        Vt = V(:,1:size(V,2)-j);
        C = Ut*St*Vt';
        rcurr = rank(C);
        % new core
        G{k} = Core(reshape(Ut,[rprev,n,rcurr]),[rprev rcurr]);
        rprev = rcurr;
        C = St*Vt';
    end
    G{N} = Core(C,[rprev 1]);
    
    tt = MPT(G);
    tt = changeCanonical(tt,frobnorm(G{N}),N);
   
end