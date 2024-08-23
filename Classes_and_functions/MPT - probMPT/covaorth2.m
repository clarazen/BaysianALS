function P = covaorth2(mpt,P,Rleft,Rright,n,dir1,dir2)
%% DESCRIPTION
%  Computes the covariance matrix after orthogonalization of
%  the current core and also computes the covariance matrix of
%  the next core where the norm is moved to
%
%   INPUT
%   mpt         MPT object
%
%   OUTPUT
%
%   August 2020, Clara Menzen
%%
    r   = mpt.ranks;
    sY  = mpt.Ysize;
    Pin = P; % for testing
    
    sP = size(P);
    P  = reshape(P,[sP(1) r(n+dir1) sY(n+dir1) r(n+dir2)]);
    P  = permute(P,[1 3 4 2]);
    P  = reshape(P,[numel(P)/r(n+dir1) r(n+dir1) ])*Rright;
    P  = reshape(P,[sP(1) sY(n+dir1) r(n+dir2) r(n+dir1)]);
    P  = permute(P,[1 4 2 3]);
    P  = Rleft*reshape(P,[r(n+dir1) numel(P)/r(n+dir1) ]);
    P  = reshape(P,[sP(1) sP(2)]);
    
    %%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     se = size(Pin,1)/r(n+dir1);
%     test = kron(eye(se),Rleft)*Pin*kron(eye(se),Rright);
%     max(max(abs(test-P)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

