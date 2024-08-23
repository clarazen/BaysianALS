function mpt = changeMPT(mpt,cores)
    mpt.cores = cores;
    mpt.N     = size(mpt.cores,2);
    mpt.ranks = 1;
    mpt.Ysize = [];
    for i = 1:mpt.N
        r = mpt.cores{i}.ranks(2);
        mpt.ranks = [mpt.ranks r];
        s = mpt.cores{i}.dims(2);
        mpt.Ysize = [mpt.Ysize s];
    end
end