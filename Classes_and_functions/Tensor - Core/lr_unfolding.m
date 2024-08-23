function G = lr_unfolding(core,rl)
%%  DESCRIPTION
%   This method computes the left- or right-unfolding of a threeway
%   Core object.
%   Left-unfolding:
%     ___                   ___
%   -|___|-  ->  R_n I_n --|___|-- R_n+1 
%      |
%
%   Right-unfolding:
%     ___               ___
%   -|___|-  ->  R_n --|___|-- I_n R_n+1 
%     |||
%
%   July 2020, Clara Menzen
%%
    Gt = core.data;
    D  = [size(Gt,1) size(Gt,3)];
    
    if     rl == "right"
        G = reshape(Gt,[D(1) numel(Gt)/D(1)]);
    elseif rl == "left"
        G = reshape(Gt,[numel(Gt)/D(2) D(2)]);
    end
end