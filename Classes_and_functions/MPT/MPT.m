classdef MPT 
%%  DESCRIBTION
%   An objetct of the class MPT a collection of N cores of the class Core 
%   that form a MPT (generalization of MPS and MPO). For an MPS, the
%   edges are I1,...,IN, for an MPO these represent a double index.
%                 ___         ___                 ___
%         R_1=1 /     \ R_2 /     \ R_3   R_N-1 /     \ R_N=1
%           ---|   G1  |---|   G2  |--- ... ---|   GN  |---
%               \ ___ /     \ ___ /             \ ___ /
%                  |           |                   |
%                  I1          I2                  IN
%   July 2020, Clara Menzen
 %%   
    properties
        cores     % Core objects
        N         % number of cores (order of approximated tensor)
        ranks     % ranks of the TT
        Ysize     % size of tensor that is decomposed
        canonical % struct (if TT is in site-k-mixed canonical format) of 
                  % - norm of the TT
                  % - number of core which contains the norm
    end
    
    methods
        function mpt = MPT(varargin)
            mpt = changeMPT(mpt,varargin{1});
            if nargin == 2
                norm     = varargin{2}.norm;
                normcore = varargin{2}.normcore;
                mpt      = changeCanonical(mpt,norm,normcore);
            end
        end
    end
end