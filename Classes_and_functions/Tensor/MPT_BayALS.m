function [pmpt,testing] = MPT_BayALS(Yt,pmpt0,sigmae,maxiter,orthog,expparams)
%%  DESCRIBTION
%   This method approximates an MPT with the alternating linear scheme in a
%   Bayesian framework.
%
%   INPUT
%   Yt          Tensor object 
%   mpt0        Prior mpt
%   sigmae      noise variance
%   maxtier     Maximum number of iterations
%   orthog      orthogonalization step true/false
%   
%   OUTPUT
%   mpt         MPT object of estimate
%
%   July 2020, Clara Menzen
%%
    N   = Yt.order;
    testing = [];
    if orthog == true && isempty(pmpt0.canonical) == true % orthoganilize core 2 to N
        for n = N:-1:2
            pmpt0 = shiftpMPTnorm(pmpt0,n,-1);
        end
        pmpt0 = changeCanonical(pmpt0,frobnorm(pmpt0.cores{1}),1);
    end
    pmpt = pmpt0;
    for iter = 1:maxiter
        switch orthog
            case true
                swipe = [1:N N-1:2];
                Dir   = [ones(1,N-1) -ones(1,N-1)];
                for k = 1:numel(swipe)
                    n   = swipe(k); 
                    dir = Dir(k);
                    pmpt.cores{n} = mP_update(pmpt,pmpt0,Yt,sigmae,n,orthog); 
                    pmpt  = shiftpMPTnorm(pmpt,n,dir);
                    pmpt0 = shiftpMPTnorm(pmpt0,n,dir);
                    if isempty(expparams) == false
                        testing = experiments(testing,pmpt,pmpt0,expparams,Yt,iter,maxiter,sigmae,n+dir);
                    end
                end
                
            case false
                for n = 1:N
                    pmpt.cores{n} = mP_update(pmpt,pmpt0,Yt,sigmae,n,orthog);
                    if isempty(expparams) == false
                        testing = experiments(testing,pmpt,pmpt,expparams,Yt,iter,maxiter,sigmae,n);
                    end
                end
        end 
    end
end