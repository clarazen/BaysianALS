function [mUT,Ap,Am] = computeUTmean(mpt,mptp,mptm,w0,w,M,eps)
    %% DESCRIPTION
    % Computes the mean with the unscented transform from the propagated
    % sigma points. Called by UT_pmpt.
    % 
    % INPUT
    % mpt       sigma point 0
    % mptp      sigma points 1...M
    % mptm      sigma points M+1...2M+1
    % w,w0      weight factors
    % M         number of elements in pmpt
    % eps       tolerance for rounding
    %
    % OUTPUT
    % mUT       mean
    % Ap        summation of sigma points with +
    % Am        summation of sigma points with -
    %
    % November 2020, Clara Menzen
    %%
    % sigma point 0
    A      = multiplyMPT(mpt,w0);               
    % sigma points 1...M
    Ap = mptp{1};
    for i = 2:M
        Ap     = addTTs({Ap, mptp{i}});                      
        wAp     = roundTT(Ap,eps);
    end
    Awp    = matrixbyvectorTT(Ap,w*ones(M,1));
    Awp    = roundTT(Awp,eps);
    % sigma points M+1...2M+1
    Am = mptm{1};
    for i = 2:M
        Am     = addTTs({Am,mptm{i}});                      
        Am     = roundTT(Am,eps);
    end
    Awm    = matrixbyvectorTT(Am,w*ones(M,1));
    Awm    = roundTT(Awm,eps);
    % sum of all weighted sigma points
    mUT    = addTTs({A,Awp,Awm});               
    mUT    = roundTT(mUT,eps);
end