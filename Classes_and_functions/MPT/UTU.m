function utu = UTU(mpt,n)
%%  DESCRIPTION
%   computes U'*U for the update equation of mean and covariance without
%   orhtogonalization step
%
%   INPUT
%   mpt         MPT object
%   OUTPUT
%   utu         U'*U
%
%   October 2020, Clara Menzen
%%
    G = mpt.cores;
    N = mpt.N;
    J = [1:n-1 n+1:N];
    C{n} = eye(mpt.Ysize(n));
    for i = 1:N-1
        j = J(i);
        U = permute(G{j}.data,[1 3 2]);
        U = reshape(U,[size(G{j}.data,1)*size(G{j}.data,3) size(G{j}.data,2)]);
        U = U*U';
        U = reshape(U,[size(G{j}.data,1) size(G{j}.data,3) size(G{j}.data,1) size(G{j}.data,3)]);
        U = permute(U,[1 3 2 4]);
        if j < n-1 || j > n+1
            U = reshape(U,[size(G{j}.data,1)*size(G{j}.data,1) size(G{j}.data,3)*size(G{j}.data,3)]);
        elseif j == n-1
            U = reshape(U,[size(G{j}.data,1)*size(G{j}.data,1) size(G{j}.data,3) size(G{j}.data,3)]);
        elseif j == n+1
            U = reshape(U,[size(G{j}.data,1) size(G{j}.data,1) size(G{j}.data,3)*size(G{j}.data,3)]);
        end
        C{j} = U;
    end
    
    if n == 1 || n == 2
        G_left = 1;
    else
        G_left = C{1};
        for i = 2:n-2
            G_left = G_left*C{i};
        end
    end

    if n == numel(C) || n == numel(C)-1 
        G_right = 1;
    else
        G_right = C{n+2};
        for i = n+3 : numel(C)
            G_right = G_right * C{i};
        end
    end 

    if n > 1
        G_left = G_left*reshape(C{n-1},[size(C{n-1},1) size(C{n-1},2)*size(C{n-1},3)]);
        G_left = reshape(G_left,[size(C{n-1},2) size(C{n-1},3)]);
    end
    if n < N
        G_right = reshape(C{n+1},[size(C{n+1},1)*size(C{n+1},2) size(C{n+1},3)])*G_right;
        G_right = reshape(G_right,[size(C{n+1},1) size(C{n+1},2)]);
    end
    utu = kron(kron(G_right',eye(mpt.Ysize(n))),G_left); 

end