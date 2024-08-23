function ttsum = addTTs(ttsummands)
%   Computes the cores of the sum of multiple tensors in TT format given by the
%   cores
%
%   INPUT
%   ttsummands  TT array of multiple TTs for summation
%
%   OUTPUT
%   ttsum       TT object which is the sum of tensors
%
%   Resources:
%   V. Oseledets: Tensor-Train Decomposition, 2011, p.2308
%
%   August 2020, Clara Menzen
%%   
    nos = numel(ttsummands); % number of summands
    N   = ttsummands{1}.N;
    for n = 1:N
        if n == 1
            Cdata = ttsummands{1}.cores{1}.data;
            for j = 2:nos
                smd   = ttsummands{j}.cores{1}.data;
                Cdata = cat(3, Cdata, smd);
            end
            C{n} = Core(Cdata,[1,size(Cdata,3)]);
        elseif n == N 
            Cdata = ttsummands{1}.cores{N}.data;
            for j = 2:nos
                Cdata = cat(1, Cdata, ttsummands{j}.cores{N}.data);
            end
            C{n} = Core(Cdata,[size(Cdata,1),size(Cdata,3)]);
        else
            sz = ttsummands{1}.cores{n}.dims;
            for j = 2:nos
                s = ttsummands{j}.cores{n}.dims;
                sz = [sz(1)+s(1),s(2),sz(3)+s(3)];
            end
            Cdata = zeros(sz);
            i1 = 1; j1 = 1; 
            i2 = ttsummands{1}.cores{n}.dims(1);
            j2 = size(ttsummands{1}.cores{n}.data,3);
            for j = 1:nos
                Cdata(i1:i2,:,j1:j2) = ttsummands{j}.cores{n}.data;  
                if j < nos
                    i1 = i2+1; i2 = i2 + ttsummands{j+1}.cores{n}.dims(1);
                    j1 = j2+1; j2 = j2 + size(ttsummands{j+1}.cores{n}.data,3);
                end
            end
            C{n} = Core(Cdata,[size(Cdata,1),size(Cdata,ndims(Cdata))]);
        end
    end
    ttsum = MPT(C);
end