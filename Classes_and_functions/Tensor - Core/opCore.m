function op = opCore(core)
%% DESCRIPTION
% Computes the outer product between two cores
% used in computeUTcova.m
    
    Sn = core.data;
    sn = core.dims;
    rn = core.ranks;
    Sn = reshape(Sn,[rn(1)*sn(2)*rn(2) 1]);
    Sn = Sn*Sn';
    Sn = reshape(Sn,[rn(1) sn(2) rn(2) rn(1) sn(2) rn(2)]);
    Sn = permute(Sn,[1 4 2 5 3 6]);
    op = Core(reshape(Sn,[rn(1)^2 sn(2)^2 rn(2)^2]),[rn(1)^2 rn(2)^2]);

end
