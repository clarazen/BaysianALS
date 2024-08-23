function A = contractN(core1,core2,contrdimindex1,contrdimindex2)
%%  DESRIBTION
%   computes contraction over multiple indexes specified in contrdimindex1
%   for core 1 and contrdimindex2 for core2
%   INPUT
%   core 1/2            Core object
%   contrdimindex1/2    specifies which dimensions should be contracted over.
%                       The user needs to be careful that the right numbers are
%                       together
%   OUTPUT
%   A                   tensor of the remaining dimensions
%
%   EXAMPLE 
%   core1: 2 x 4 x 3
%   core2: 4 x 5 x 2 x 6
%   contrdimindex1 = [1 2]
%   contrdimindex2 = [3 1]
%   First, core 1 requires permuting to 3 x 2 x 4 
%   and core 2 to 2 x 4 x 5 x 6
%   Second, both cores are reshaped such that the dimensions of the 
%   matrix multiplication are (3 x 8) x (8 x 30). 
%   Third, the resulting matrix is reshaped into the remaining dimensions
%   which are in this case 3 x 5 x 6
%                 ___                ___
%               /     \            /     \ 
%           ---|       |---    ---|       |---              
%               \ ___ /            \ ___ /   
%                | | |             / | | \
%                2 4 2            4  5 2  6
%   Permute:
%                 ___                ___
%               /     \            /     \ 
%           ---|       |---    ---|       |---              
%               \ ___ /            \ ___ /   
%                | | |             / | | \
%                3 2 4            2  4 5  6 
%   Reshape and contract:
%                  3
%                 _|_               
%               /     \           
%           ---|       |---                
%               \ ___ /                  ___
%                  |                   /     \
%                 2*4          ->  ---|       |--- 
%                 _|_                  \ ___ /
%               /     \                  /|\
%           ---|       |---             3 5 6  
%               \ ___ /             
%                  |             
%                 6*5 
% 
%   February 2020, Clara Menzen
%%

    dim1 = core1.dims;
    order1 = core1.order;
    dim2 = core2.dims;
    order2 = core2.order;
    contrdim1 = dim1(contrdimindex1);
    contrdim2 = dim2(contrdimindex2);
    otherdim1 = dim1;
    otherdim1(contrdimindex1) = [];
    otherdim2 = dim2;
    otherdim2(contrdimindex2) = [];

    % permute if necesary
    newperm1tmp = 1:order1;
    newperm1tmp(contrdimindex1) = [];
    newperm1 = [newperm1tmp contrdimindex1];
    A1 = permute(core1.data,newperm1);

    newperm2tmp = 1:order2;
    newperm2tmp(contrdimindex2) = [];
    newperm2 = [contrdimindex2 newperm2tmp];
    A2 = permute(core2.data,newperm2);

    % reshape 
    A1 = reshape(A1, [numel(A1)/prod(contrdim1) prod(contrdim1)] );
    A2 = reshape(A2, [prod(contrdim2) numel(A2)/prod(contrdim2)]);

    % summation over middle indexes
    A = A1*A2;

    % release the outer dimensions
    A = Tensor(reshape(A,[otherdim1 otherdim2]));

end