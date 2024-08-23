function Ghat = updateCore(mpt,Yt,n)
%%  DESCRIPTION
%   updates the nth core while keeping all other cores fixed
%
%   INPUT
%   mptt        MPT object
%   Yt          Data tensor (Tensor object)
%   n           number of core that is being updated
%   OUTPUT
%   coreupdate  updated core
%
%   January 2020, Clara Menzen
%%
    G = mpt.cores;
    N = mpt.N;
    for j = 1:N-n
        Gn{j} = G{n+j}; 
    end
    for j = 1:n-1
       Gn{N-n+j} = G{j}; 
    end
    % Permute X to In...IN I1...In-1
    Ghat = Tensor(permute(Yt.data,[n:N 1:n-1]));
    % (n+1)th core
    Ghat = contractN(Gn{1},Ghat,2,2);
    % (n+2)th core
    if ismatrix(Gn{1}.data)
        Gdata = reshape(Ghat.data,[Ghat.dims(1) 1 Ghat.dims(2:end)]);
        Ghat = changeTdata(Ghat,Gdata);
    end
    Ghat = contractN(Gn{2},Ghat,[1 2],[2 4]);
    if n+2 == N
       if ismatrix(Gn{2})
           Gdata = reshape(Ghat.data,[1 Ghat.dims]);
           Ghat  = changeTdata(Ghat,Gdata);
       end
    end
    % (n+3)th - (n-1)th core
    for i = 3:N-1
        Ghat = contractN(Gn{i},Ghat,[1 2],[1 4]);
        if ismatrix(Gn{i}.data)
            Gdata = reshape(Ghat.data,[1 Ghat.dims]);
            Ghat  = changeTdata(Ghat,Gdata);
        end
    end
    Gdata = permute(Ghat.data,[1 3 2]);
    ranks = mpt.cores{n}.ranks;
    Ghat  = changeCdata(G{n},Gdata,ranks);

end