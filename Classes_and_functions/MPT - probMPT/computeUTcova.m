function PUT = computeUTcova(mpt,mUT,w0,w,M,Ap,Am,eps)
    %% DESCRIPTION
    % Computes the mean with the unscented transform from the propagated
    % sigma points. Called by UT_pmpt.
    % 
    % INPUT
    % mpt       sigma point 0
    % mUT       mean computed with UT
    % w,w0      weight factors
    % M         number of elements in pmpt
    % Ap        summation of sigma points with +
    % Am        summation of sigma points with -
    % eps       tolerance for rounding
    %
    % OUTPUT
    % PUT       covariance matrix 
    %
    % November 2020, Clara Menzen
    %%
    m_mat = vectorbyvectorTT(mUT,ones(M,1));
    s     = mUT.cores{1}.dims;
    N     = mpt.N;
    % sigma point 0
    op    = addTTs({mpt,multiplyMPT(mUT,-1)});             
    op    = roundTT(op,eps);
    op    = multiplyMPT(outerprodTT(op,op),w0); 
    op    = roundTT(op,eps);
    % sigma points 1...M
    sum1  = addTTs({Ap,multiplyMPT(m_mat,-1)});             
    sum1  = roundTT(sum1,eps);
    %%
%     cm = corebymatrixTT(sum1.cores{1},sqrt(w)*eye(M),s);
%     sum1.cores{1} = cm;
%     for i = N:2
%         sum1 = shiftMPTnorm(sum1,i,-1);
%     end
%     
%     tmp = reshape(sum1.cores{1}.data,[size(sum1.cores{1}.data,2)/M M size(sum1.cores{1}.data,3)]);
%     tmp = permute(tmp,[1 3 2]);
%     C1 = reshape(tmp,[size(sum1.cores{1}.data,2)/M*size(sum1.cores{1}.data,3) M]);
%     sv = svd(C1);
%     for k = 1:numel(sv)
%         if sv(k) < 1e-5
%             break
%         end
%     end
%     [U,~,~] = svds(C1,k-1);
    
    G{1}  = corebymatrixbycoreTT(sum1.cores{1},w*eye(M),s);
    for n = 2:N
        G{n} = opCore(sum1.cores{n});
    end
    opp   = MPT(G);
    opp   = roundTT(opp,eps);
    % sigma points M+1...2M+1
    sum2  = addTTs({Am,multiplyMPT(m_mat,-1)});             
    sum2  = roundTT(sum2,eps);
    G{1}  = corebymatrixbycoreTT(sum2.cores{1},w*eye(M),s);
    for n = 2:N
        G{n} = opCore(sum2.cores{n});
    end
    opm   = MPT(G);
    opm   = roundTT(opm,eps);
    % sum of all weighted outer products
    PUT   = addTTs({op,opp});   
    PUT   = roundTT(PUT,eps);
    PUT   = addTTs({PUT,opm});   
    PUT   = roundTT(PUT,eps);
end