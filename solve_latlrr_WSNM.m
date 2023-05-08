function [Z,L,E] = solve_latlrr_WSNM(X,beta,alpha,rho,C1,C2,p)
%LATLRR Summary of this function goes here
% This routine solves the optmization problem of WSN-LLRR (noisy model)
% min |Z|^{p}_w,sp + alpha*|L|^{p}_w,sp + beta*|E|_1
% s.t., X = XZ + LX + E  
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        alpha -- usually alpha = 1
%        beta  -- this parameter depends on the noise level of data
if nargin<4
    rho = 1.12;
end
if nargin<3
    alpha = 1.0;
end
if nargin<2
    beta = 1.0; 
end
%% 
Q1 = orth(X');
Q2 = orth(X);
A = X*Q1;
B = Q2'*X;
[Z,L,E] = solve_latlrra_WSNM(X,A,B,beta,alpha,rho,C1,C2,p);
Z = Q1*Z;
L = L*Q2';