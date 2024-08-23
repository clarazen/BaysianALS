clear; close all; clc;

r = 3; n = 5;
for i = 1:1000   
    G{1} = Core(randn(1,n,r),[1 r]); 
    G{2} = Core(randn(r,n,r),[r r]); 
    G{3} = Core(randn(r,n)  ,[r 1]); 
    Truthtt{i} = MPT(G);
    Truth{i}   = approxTensor(Truthtt{i});
    truth{i}   = reshape(Truth{i}.data, [numel(Truth{i}.data),1]);
end

clear r n G i
save('truth')