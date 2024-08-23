function tensor = changeTdata(tensor,Tdata)
    tensor.data  = Tdata;
    tensor.dims  = size(Tdata);
    tensor.order = numel(tensor.dims);
end