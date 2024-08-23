classdef probCore < Core
    properties
        mean
        cova
    end
    methods
        function pcore = probCore(Tdata,ranks,mean,cova)
            pcore@Core(Tdata,ranks)
            pcore.mean = mean;
            pcore.cova = cova;
        end
    end
end