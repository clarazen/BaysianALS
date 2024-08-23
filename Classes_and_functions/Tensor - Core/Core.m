classdef Core < Tensor
    properties
        ranks
        %ortho
    end
    methods
        function core = Core(Tdata,ranks)
            core@Tensor(Tdata)
            core.ranks = [ranks(1) ranks(2)];
            %core.ortho = lr_orthog(core);
        end
    end
end