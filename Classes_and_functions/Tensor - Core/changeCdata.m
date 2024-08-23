function core = changeCdata(core,Cdata,ranks)

    core       = changeTdata(core,Cdata);
    core.ranks = [ranks(1) ranks(2)];
    %core.ortho = lr_orthog(core);
    
end