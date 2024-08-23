clear; format compact; format short; clc; close all;


Xt = randn(3,5,3);
A  = randn(3,3);

B_true1 = reshape( A*reshape(Xt,[3 15]) ,[3 5 3] );
B_test1 = nmodeprod( Tensor(Xt),A,1 );
B_true1-B_test1

B_true2 = reshape( reshape(Xt,[15 3])*A ,[3 5 3] );
B_test2 = nmodeprod( Tensor(Xt),A',3 );
B_true2-B_test2