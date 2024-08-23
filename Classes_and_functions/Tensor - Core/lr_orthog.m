function orth = lr_orthog(core)
%%   DESCRIPTION
%   This method checks if a core is left/right/not orthogonal.
%   Left-orthogonal means that the left-unfolding Gl of a core satisfies 
%   Gl'*Gl = I (R_n+1 x R_n+1).  
%   Right-orthogonal means that the right-unfolding Gr of a core satisfies 
%   Gr*Gr' = I (R_n x R_n). 
%   It uses the property of an orthogonal matrix that all singular values 
%   are one, so the sum of the singular values must be equal to the size 
%   of the matrix.
%
%   INPUT: core, OUTPUT: orthogonality
%
%   Resources: 
%   Holtz et al.: The Alternating Linear Scheme for Tensor 
%   Optimisation in the TT Format, 2012. p.7 eqn. (2.13) (2.14)
%
%   July 2020, Clara Menzen
%%
    SVright = svd(lr_unfolding(core,'right'));
    SVleft  = svd(lr_unfolding(core,'left'));

    if abs(sum(SVright) - size(SVright,1)) < 1e-10
        orth = "right orthogonal"; 
    elseif abs(sum(SVleft) - size(SVleft,1)) < 1e-10
        orth = "left orthogonal";  
    else
        orth = "not orthogonal";   
    end
end