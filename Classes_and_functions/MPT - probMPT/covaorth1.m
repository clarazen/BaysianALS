function P = covaorth1(pmpt,P,Rleft,Rright,n,dir1,dir2)
%% DESCRIPTION
%  Computes the covariance matrix after orthogonalization of
%  the current core and also computes the covariance matrix of
%  the next core where the norm is moved to
%
%   INPUT
%   pmpt         probMPT object
%   OUTPUT
%
%   August 2020, Clara Menzen  
%%
    r   = pmpt.ranks;
    sY  = pmpt.Ysize;
    Pin = P; %for testing
    
    sP  = size(P);
    P   = reshape(P,[ numel(P)/r(n+dir1) r(n+dir1)])*Rright;
    P   = reshape(P,[ r(n+dir2) sY(n+dir2) r(n+dir1) sP(2) ]);
    P   = permute(P,[3 1 2 4]);
    P   = Rleft*reshape(P,[ r(n+dir1) numel(P)/r(n+dir1) ]);
    P   = reshape(P,[ r(n+dir1) r(n+dir2) sY(n+dir2) sP(2) ]);
    P   = permute(P,[2 3 1 4]);
    P   = reshape(P,[ sP(1) sP(2)]);
    
    %%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     se = size(Pin,1)/r(n+dir1);
%     test = kron(Rleft,eye(se))*Pin*kron(Rright,eye(se));
%     max(max(abs(test - P)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
end


