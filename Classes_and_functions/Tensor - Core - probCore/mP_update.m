function Ghat = mP_update(mpt,mpt0,Yt,sigmae,n,orthog)
%%  DESCRIPTION
%   Computes the update of the core in the classical ALS case and the
%   update of the mean and covariance in the Bayesian ALS case
%
%   INPUT
%   tt      TT object
%   X       Data tensor (core object)
%   sigmae  noise 
%   n       number of core that is being updated
%   method  Classical or Bayesian ALS
%
%   OUTPUT
%   tt      TT object with updated core and covariance matrix
%
%   References:
%   Sarkka: Bayesian Filtering and Smoothing, 2013. P.29, eqn. (3.4)
%
%   January 2020, Clara Menzen
%%
    
    G     = mpt.cores;
    G0    = mpt0.cores;
    invP0 = inv(G0{n}.cova);
    m0    = G0{n}.mean;
    
    als   = updateCore(mpt,Yt,n).data;
    if orthog == true
        new_P = pinv( invP0 + eye(size(invP0))/sigmae^2 );
        new_m = new_P* ( (reshape(als,numel(als),1)/sigmae^2 + invP0*m0) );
    else
        utu   = UTU(mpt,n);     
       	new_P = pinv( invP0 + utu/sigmae^2 );
        new_m = new_P* (reshape(als,numel(als),1)/sigmae^2 + invP0*m0);
    end
    
    Cdata = reshape(new_m, G{n}.dims);
    Ghat  = changeprobCdata(G{n},Cdata,G{n}.ranks,new_m,new_P);

end