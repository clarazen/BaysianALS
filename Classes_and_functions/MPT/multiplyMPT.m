function mpt = multiplyMPT(mpt,a)
%% Description
%  Function that multiplies a tensor train by a scalar
%  INPUT
%  mpt      MPT object
%  a        scalar
%%
    mpt.cores{1}.data = mpt.cores{1}.data*a;
    mpt = changeMPT(mpt,mpt.cores);
end