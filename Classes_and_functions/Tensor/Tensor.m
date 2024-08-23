classdef Tensor
    %%   DESCRIPTION
    %    An object of the Tensor class is an N-way array with dimensions
    %    I_1,...,I_N. 
    %
    %    A tensor is fully defined with the tensorial data from which each 
    %    dimension and the total order can be inferred.
    %                 ___
    %               /     \
    %              |   Y   |
    %               \ ___ /
    %                / | \
    %               /  |  \
    %             I_1,...,I_N
    %
    %    Methods:
    %    - constructor
    %    - n-mode-product, mode-n-matricization
    %    - Frobenius norm
    %    - inner product
    %    - MPT_SVD, MPT_ALS, MPT_BayALS
    %
    %   July 2020, Clara Menzen
%%
    properties
        data
        dims
        order
    end
    
    methods
        %% Contructor
        function tensor = Tensor(Tdata)
            tensor.data  = Tdata;
            tensor.dims  = size(Tdata);
            tensor.order = numel(tensor.dims);
        end
    end
end