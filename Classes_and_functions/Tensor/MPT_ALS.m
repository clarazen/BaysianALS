function [mpt,testing] = MPT_ALS(Yt,mpt0,maxiter,orthog,testparams)
%%  DESCRIBTION
%   This method approximates an MPT with the alternating linear scheme.
%
%   INPUT
%   tensor      Tensor object 
%   mpt0        Initial mpt
%   maxtier     Maximum number of iterations
%   orthog      boolean which specifies if an orthogonalization step is
%               incorporated
%   
%   OUTPUT
%   mpt         MPT object
%
%   July 2020, Clara Menzen
%%
    N   = Yt.order;
    mpt = mpt0;
    if orthog == true % orthoganilize core 2 to N
        for n = N:-1:2
            mpt = shiftMPTnorm(mpt,n,-1);
        end
        mpt = changeCanonical(mpt,frobnorm(mpt.cores{1}),1);
    end
    for iter = 1:maxiter
        swipe = [1:N N-1:2];
        Dir   = [ones(1,N-1) -ones(1,N-1)];
        for k = 1:numel(swipe)
            n  = swipe(k);
            dir = Dir(k);
            if orthog == true
                mpt.cores{n} = updateCore(mpt,Yt,n);
                mpt          = shiftMPTnorm(mpt,n,dir);
            else
                utu   = UTU(mpt,n);
                als   = updateCore(mpt,Yt,n).data;
                alsr  = reshape(als,[numel(als) 1]);
                Cdata = reshape(pinv(utu)*alsr,size(als));
                mpt.cores{n} = changeCdata(mpt.cores{n},Cdata,[mpt.ranks(n) mpt.ranks(n+1)]);
            end
        end 
        if isempty(testparams) == false
            yhat       = reshape(approxTensor(mpt).data,[numel(testparams.truth) 1]);
            er         = norm(yhat-testparams.truth)/norm(testparams.truth);
            testing.er = er;
        else
            testing = [];
        end
    end
end