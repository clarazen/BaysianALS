function ip = innerprodTT(ttA,ttB)
%%  DESCRIBTION
%   Computes the inner product for two same sized tensors in TT-format
%
%   INPUT
%   A        TT objetc
%   B        TT object
%   OUTPUT
%   AB       inner product of tt 1 and 2
%
%   December 2020, Clara Menzen
%%
    A = ttA.cores;
    B = ttB.cores;
    N = ttA.N;
    ip = 1;
    
    for n = 1:N
        Bn = permute(B{n}.data,[1 3 2]);
        Bn = reshape(Bn,[B{n}.dims(1)*size(B{n}.data,3) B{n}.dims(2)]);
        An = permute(A{n}.data,[1 3 2]);
        An = reshape(An,[A{n}.dims(1)*size(A{n}.data,3) A{n}.dims(2)]);
        BA = reshape(Bn*An',[B{n}.dims(1) size(B{n}.data,3) A{n}.dims(1) size(B{n}.data,3)]);
        BA = permute(BA,[3 1 4 2]);
        BA = reshape(BA,[B{n}.dims(1)*size(A{n}.data,1) size(B{n}.data,3)*size(A{n}.data,3)]);
        ip = ip*BA;
    end
    
end