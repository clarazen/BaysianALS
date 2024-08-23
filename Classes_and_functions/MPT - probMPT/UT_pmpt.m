function [mUT,PUT] = UT_pmpt(pmpt,epsm,epsP,Pplease)
    %% DESCRIPTION
    %   Computes an approximation of the joint distribution of given random
    %   variables, which are the vectorized stacked TT cores. The non-linear
    %   function is the mapping from the TT cores to the tensor which the TT
    %   approximates
    %
    %   INPUT
    %   mpt        MPT object
    %   OUTPUT
    %   mUT        approximate mean of joint distribution as MPT object
    %   PUT        approximate covariance of joint distribution
    %
    %   Resources:
    %   Sarkka: Bayesian Filtering and Smoothing, 2013. p.81-84
    %
    %   August-November 2020, Clara Menzen
    %%        
    %% Create 2n+1 Sigma Points
    N = pmpt.N;
    for n = 1:N
        numelcore(n) = numel(pmpt.cores{n}.data);    
    end
    M = sum(numelcore);
    
    % parameters
    kappa  = 3-M; beta = 2; alpha = .001;
    lambda = alpha^2*(M+kappa) - M;
    % weights
    wm0 = lambda/(M+lambda);
    wP0 = lambda/(M+lambda) + (1-alpha^2+beta);
    w   = 1/(2*(M+lambda));
    
    i1 = 1;
    i2 = numelcore(1);
    for n = 1:N
        m              = pmpt.cores{n}.mean*ones(1,M);
        cholP          = zeros(size(m,1),M);
        cholP(:,i1:i2) = chol(pmpt.cores{n}.cova)';
        Xp{n,1}        = m + sqrt(M+lambda)*cholP;
        Xm{n,1}        = m - sqrt(M+lambda)*cholP;
        if n<N
            i1 = i1 + numelcore(n);
            i2 = i2 + numel(pmpt.cores{n+1}.data);
        end        
    end
       
    %% Propagate sigma points thru non-linear function
    for i = 1:M
        for n = 1:N
            ranks = pmpt.cores{n}.ranks;
            if n == 1
                e = zeros(M,1); e(i) = 1;
                sizes = [size(pmpt.cores{1}.data,1)...
                         size(pmpt.cores{1}.data,2)*M...
                         size(pmpt.cores{1}.data,3)];
                     
                Cdata = Xp{1,1}(:,i)*e';     
                Cdata = reshape(Cdata,[size(pmpt.cores{1}.data) M]);
                Cdata = permute(Cdata,[1 2 4 3]);
                Cdata = reshape(Cdata,sizes);
                Gp{1} = Core(Cdata,ranks);
                
                Cdata = Xm{1,1}(:,i)*e';
                Cdata = reshape(Cdata,[size(pmpt.cores{1}.data) M]);
                Cdata = permute(Cdata,[1 2 4 3]);
                Cdata = reshape(Cdata,sizes);
                Gm{1} = Core(Cdata,ranks);
                               
            else
                Cdata = reshape(Xp{n,1}(:,i),size(pmpt.cores{n}.data));
                Gp{n} = Core(Cdata,ranks);
                Cdata = reshape(Xm{n,1}(:,i),size(pmpt.cores{n}.data));
                Gm{n} = Core(Cdata,ranks);
            end
        end
        mptp{i} = MPT(Gp);
        mptm{i} = MPT(Gm);
    end
     
    [mUT,A1,A2] = computeUTmean(pmpt,mptp,mptm,wm0,w,M,epsm);
    
    if Pplease == true
        PUT = computeUTcova(pmpt,mUT,wP0,w,M,A1,A2,epsP);
    else
        PUT = [];
    end
end