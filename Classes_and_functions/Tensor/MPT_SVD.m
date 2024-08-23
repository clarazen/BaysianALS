function mpt = MPT_SVD(tensor,dims,eps)
%%  DESCRIPTION
%   Method that computes an MPT (generalisation of MPS size(dims,1)=1 
%   and MPO size(dims,1)=2 to size(dims,1)=n) for a given tensor with given
%   accuracy eps.
%
%   INPUT
%   tensor   Tensor object
%   dims     desired dimensions for each of the MPT cores
%            e.g. [I1 J1; I2 J2; I3 J3] see example below
%   eps      accuracy
%   OUTPUT
%   mpt      Returns an object of the class MPT
%
%   Example: large matrix: N x M = (I1*I2*I3 x J1*J2*J3)
%   ___         ____              ____                 ___ 
%  /   \  ==>  /    \   ==>      /    \      ==>      /   \
%  \___/  (1)  \____/   (2)      \____/      (3)      \___/ 
%   | |        /||||\            /||||\               / | \
%   | |       / |||| \          / |||| \             /  |  \
%   N M   I1 I2 I3 J1 J2 J3  I1 J2 I2 J2 I3 J3   I1*J1 I2*J2 I3*J3  
%         ___    ___    ___      ___    ___    ___
%   ==>  /   \__/   \__/   \    /   \__/   \__/   \
%   (4)  \___/  \___/  \___/    \___/  \___/  \___/
%          |      |      |       | |    | |    | |
%        I1*J1  I2*J2  I3*J3    I1 J1  I2 J2  I3 J3
%
%   Resources:
%   [1] V. Oseledets: Approximation of 2^d x 2^d matrices using 
%       Tensor Decomposition, 2010, p.2136-2138
%
%   July 2020, Clara Menzen
%%
    Yt = tensor.data;

    if numel(dims) > numel(tensor.dims)
        % (1) reshape given obeject into tesor with the specified dimensions
        newdims = reshape(dims,[1,numel(dims)]);
        Yt  = reshape(Yt,newdims);

        % (2) permute indices for same core next to each other
        newperm = [];
        for i = 1:size(dims,1)
            newperm = [newperm i  i + size(dims,1)];
        end
        Yt = permute(Yt,newperm);

        % (3) create long index
        Yt = reshape(Yt,[prod(dims,2)]');
        tensor = Tensor(Yt);
    end

    % (4) compute TT
    mpt = TT_SVD(tensor,eps);

end