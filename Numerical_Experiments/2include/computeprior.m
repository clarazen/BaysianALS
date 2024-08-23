function pmpt0 = computeprior(sigma,a,G)
    %% Prior mean and covariance
    N = size(G,2);
    r = size(G{1}.data,3);
    ranks = [1 r r 1];
    for i = 1 : N
        M  = numel(G{i}.data);
        m0 = reshape(G{i}.data,[M 1]) + a*randn(M,1);
        m0 = m0/norm(m0)*(norm(reshape(G{i}.data,[M 1])));
        P0 = sigma^2*eye(M);
        G0 = reshape(m0,size(G{i}.data));
        prior.cores{i} = probCore(G0,[ranks(i) ranks(i+1)],m0,P0);
    end
    pmpt0 = probMPT(prior.cores);
end