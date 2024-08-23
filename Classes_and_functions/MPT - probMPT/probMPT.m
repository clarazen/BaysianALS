classdef probMPT < MPT 
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
    end
    
    methods
        function pmpt = probMPT(varargin)
            pmpt@MPT(varargin{1})
        end
    end
end