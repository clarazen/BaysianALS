function mpts = shiftMPTnorm(mpt,n,dir)
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
    % June 2020, Clara Menzen
    %%
        mpts = mpt;
        Gt   = mpts.cores;
        if dir == 1
            Gl = lr_unfolding(Gt{n},'left');
            ind   = 1;
            [Q,R] = qr(Gl,0);
        elseif dir == -1
            Gr = lr_unfolding(Gt{n},'right');
            ind  = 3;
            [Q,R] = qr(Gr',0);
            Q = Q';
        end

        ranks             = mpts.cores{n}.ranks;
        mpts.cores{n}     = changeCdata(mpts.cores{n}, reshape(Q, Gt{n}.dims),ranks);
        ranks             = mpts.cores{n+dir}.ranks;
        mpts.cores{n+dir} = changeCdata(mpts.cores{n+dir}, nmodeprod(Gt{n+dir},R,ind),ranks);
        mpts              = changeMPT(mpts,mpts.cores);

        if isempty(mpts.canonical) == false
            mpts = changeCanonical(mpts,frobnorm(mpts.cores{n+dir}),n+dir);
        end
    end