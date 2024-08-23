function testing = experiments(testing,pmpt,pmpt0,expparams,Yt,iter,maxiter,sigmae,n)

    if expparams.expnr == 1
        
        mUT      = UT_pmpt(pmpt,expparams.eps_m ,expparams.eps_P,false);
        yhat     = reshape(approxTensor(mUT).data,[numel(expparams.truth) 1]);
        y        = reshape(Yt.data,[numel(Yt.data) 1]);
        er_truth = norm(yhat-expparams.truth)/norm(expparams.truth);
        er_meas  = norm(yhat-y)/norm(y);
        %
        Py       = sigmae^2*eye(numel(y));
        likprior = log(Gaussian_prob(y,yhat,Py));
        for n = 1:pmpt.N
            xn        = pmpt0.cores{n}.mean;
            mn        = pmpt.cores{n}.mean;
            Pn        = pmpt0.cores{n}.cova;
            priorprob = log(Gaussian_prob(xn,mn,Pn));
            likprior  = likprior+priorprob;
        end
        
        testing.er_truth(iter) = er_truth;
        testing.er_meas(iter)  = er_meas;
        testing.Likprior(iter) = likprior;
    
    elseif expparams.expnr == 2
        
%         [~,Phat] = UT_pmpt(pmpt,expparams.eps_m,expparams.eps_P,true);
%         PUT      = approxMatrix(Phat);
%         normP    = norm(PUT,'fro');
%         traceP   = trace(PUT);  
        normPn  = norm(pmpt.cores{n}.cova,'fro');
        tracePn = trace(pmpt.cores{n}.cova);
        testing.normPn(iter,n)  = normPn;
        testing.tracePn(iter,n) = tracePn;
%         testing.normP(iter)     = normP;
%         testing.traceP(iter)    = traceP;
%         testing.P{iter}         = PUT;
        
    elseif expparams.expnr == 3.1 && iter == maxiter
        
        [mUT,PUT]   = UT_pmpt(pmpt,expparams.eps_m ,expparams.eps_P,true);
        yhat        = reshape(approxTensor(mUT).data,[numel(expparams.truth) 1]);
        er          = norm(yhat-expparams.truth)/norm(expparams.truth);
        testing.er  = er;
        testing.mUT = mUT;
        testing.PUT = PUT;
        
    elseif expparams.expnr == 3 && iter == maxiter
        
        mUT         = UT_pmpt(pmpt,expparams.eps_m,expparams.eps_P,false);
        yhat        = reshape(approxTensor(mUT).data,[numel(expparams.truth) 1]);
        er          = norm(yhat-expparams.truth)/norm(expparams.truth);
        testing.er  = er;
        testing.mUT = mUT;
        
    elseif expparams.expnr == 4 && iter == maxiter
        
        [mUT,PUT]   = UT_pmpt(pmpt,expparams.eps_m ,expparams.eps_P,true);
        testing.PUT = PUT;
        testing.mUT = mUT;  
        
    else
        testing = [];
    end
   
    
      
end