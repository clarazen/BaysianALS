function pmpt = shiftpMPTnorm(pmpt,n,dir)
    %% DESCRIPTION
    % Shifts norm in a MPT in site-n-mixed-canonical form
    % to the next (dir = 1) or previous (dir = -1) core.
    % Also the covariance matrices are updated
    % INPUT
    % mpt   MPT object 
    % n     not orthogonal core
    % dir   direction in which norm should be shifted (right or left)
    % OUTPUT
    % mpt   MPT object with norm shifted one core to the left or right
    %
    %%
    Gt  = pmpt.cores;
    if dir == 1
        Gl = lr_unfolding(Gt{n},'left');
        ind   = 1;
        [Q,R] = qr(Gl,0);
    elseif dir == -1
        Gr    = lr_unfolding(Gt{n},'right');
        ind   = 3;
        [Q,R] = qr(Gr',0);
        Q = Q';
    end

    %% Change mean according to the orthogonalization in the tensor framework
    newG_o    = reshape(Q, Gt{n}.dims);
    newmean_o = reshape(newG_o,[numel(newG_o),1]);
    G_next    = nmodeprod(Gt{n+dir},R,ind);
    mean_next = reshape(G_next,[numel(G_next),1]);

    %% Change covariance matrix according ot the orthogonalization in the tensor framework
    if dir == 1
        % covariance matrix of updated core
        newcova_o = covaorth1(pmpt,pmpt.cores{n}.cova,pinv(R'),pinv(R),n,1,0);
        % covariance matrix of next core to the right
        cova_next = covaorth2(pmpt,pmpt.cores{n+1}.cova,R,R',n,1,2);
    elseif dir == -1
        % covariance matrix of updated core
        newcova_o = covaorth2(pmpt,pmpt.cores{n}.cova,pinv(R'),pinv(R),n,0,1);
        % covariance matrix of next core to the left
        cova_next = covaorth1(pmpt,pmpt.cores{n-1}.cova,R,R',n,0,-1);
    end

    r1 = pmpt.cores{n}.ranks;
    r2 = pmpt.cores{n+dir}.ranks;
    pmpt.cores{n}     = changeprobCdata(pmpt.cores{n},newG_o,r1,newmean_o,newcova_o);
    pmpt.cores{n+dir} = changeprobCdata(pmpt.cores{n+dir},G_next,r2,mean_next,cova_next);

    if isempty(pmpt.canonical) == false
        pmpt = changeCanonical(pmpt,frobnorm(pmpt.cores{n+dir}),n+dir);
    end

end