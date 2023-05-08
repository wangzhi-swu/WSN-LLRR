function [Z,L,E] = solve_latlrra_WSNM(X,A,B,beta,alpha,rho,C1,C2,p)
%LATLRR Summary of this function goes here
% This routine solves the optmization problem of WSN-LLRR (noisy model)
% min |Z|^{p}_w,sp + alpha*|L|^{p}_w,sp + beta*|E|_1
% s.t., X = AZ + LB + E  
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        alpha -- usually alpha = 1
%        beta  -- this parameter depends on the noise level of data
%%
tol = 1e-6; 
maxIter = 500;

[d n] = size(X);
nA = size(A,2);
dB = size(B,1);
max_mu = 1e10; 
mu = 1e-8;   

ata = A'*A;
bbt = B*B';
inv_ata = inv(ata+eye(nA));
inv_bbt = inv(bbt+eye(dB));
%% Initializing optimization variables
% intialize
J = zeros(nA,n);
Z = zeros(nA,n);
S = zeros(d,dB);
L = zeros(d,dB);
E = sparse(d,n);

Y1 = zeros(d,n);
Y2 = zeros(nA,n);
Y3 = zeros(d,dB);
%% Opt.X
[U1, S1, V1] = svd(X, 'econ');
siagS1 = diag(S1);
    
Opt.siagS1 = siagS1;
Opt.d = d;
Opt.n = n;

% Opt.m = m;
%% Start main loop
display = true;
iter = 0;
if display
    disp(['initial,r(Z)=' num2str(rank(J)) ',r(L)=' num2str(rank(L)) ',|E|_1=' num2str(sum(sum(abs(E))))]);
end
while iter<maxIter
    iter = iter + 1;
    %update J
    temp = Z + Y2/mu;
	[U, sigma, V] = svd(temp, 'econ');
    diagS = diag(sigma);
    
    [tempDiagS, svp] = WSNM(Opt, diagS, C1/mu, p);                %%WSNM: wi, p
    J = U(:,1:svp)*diag(tempDiagS(1:svp))*V(:,1:svp)';
	rZ = svp; %rank(J);

    %update S
    temp = L + Y3/mu;
	[U, sigma, V] = svd(temp, 'econ');
    diagS = diag(sigma);

	[tempDiagS, svp] = WSNM(Opt, diagS, C2/mu, p);     %%WSNM: wi, p
    S = U(:,1:svp)*diag(tempDiagS(1:svp))*V(:,1:svp)';
    rL = svp; %rank(S);

    %udpate Z
    Z = inv_ata*(A'*X-A'*L*B-A'*E+J+(A'*Y1-Y2)/mu);
    
    %update L
    L = (X*B'-A*Z*B'-E*B'+S+(Y1*B'-Y3)/mu)*inv_bbt; 
    
    %update E
    temp = X-A*Z-L*B+Y1/mu;
    E = max(0,temp - beta/mu)+min(0,temp + beta/mu);

    %update the multiplies
    leq1 = X-A*Z-L*B-E;
    leq2 = Z-J;
    leq3 = L-S;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
    stopC = max(max(max(abs(leq3))),stopC);
    if display && (iter==1 || mod(iter,50)==0 || stopC<tol)
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ',stopALM=' num2str(stopC,'%2.3e') ...
            ',rank(Z)=' num2str(rZ) ',rank(L)=' num2str(rL) ',|E|_1=' num2str(sum(sum(abs(E))))]);
    end
    
    if stopC<tol
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        Y3 = Y3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
end
