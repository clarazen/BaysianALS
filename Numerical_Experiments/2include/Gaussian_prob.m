function p = Gaussian_prob(x,m,P)
    
    p = 1/(sqrt(2*pi)^numel(x)*det(P)) * exp(-0.5*(x-m)'*inv(P)*(x-m));

end