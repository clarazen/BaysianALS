function pcore = changeprobCdata(pcore,Cdata,ranks,mean,cova)

    pcore      = changeCdata(pcore,Cdata,ranks);
    pcore.mean = mean;
    pcore.cova = cova;
    
end