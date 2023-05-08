function [] = WSNLLRR_yaleb10()
load('yaleb10.mat');
X = obj.X;
gnd = obj.cids;
K = max(gnd);

%run WSN-LLRR
for p = 0.65
    for lambda =0.2
        for C1 = 0.01
            for C2 = C1
                for rho = 1.5
                    for alpha = 4
                        tic;
                        Z = solve_latlrr_WSNM(X, lambda, 1, rho, C1, C2, p);
                        
                        %post processing
                        [U,S,V] = svd(Z,'econ');
                        S = diag(S);
                        r = sum(S>1e-4*S(1));
                        U = U(:,1:r);S = S(1:r);
                        U = U*diag(sqrt(S));
                        U = normr(U);
                        L = (U*U').^(2*alpha);
                        
                        % spectral clustering
                        D = diag(1./sqrt(sum(L,2)));
                        L = D*L*D;
                        [U,S,V] = svd(L);
                        V = U(:,1:K);
                        V = D*V;
                        
                        n = size(V,1);
                        M = zeros(K,K,20);
                        rand('state',123456789);
                        for i=1:size(M,3)
                            inds = false(n,1);
                            while sum(inds)<K
                                j = ceil(rand()*n);
                                inds(j) = true;
                            end
                            M(:,:,i) = V(inds,:);
                        end
                        idx = kmeans(V,K,'emptyaction','singleton','start',M,'display','off');
                        
                        grps = bestMap(gnd,idx);
                        missrate = sum(gnd(:) ~= grps(:)) / length(gnd);
                        acc = 1 - missrate;
                        
                        NMI   = nmi(gnd(:), grps(:));
                        [AR]=RandIndex(gnd(:),grps(:));
                        
                        disp(['lambda=' num2str(lambda) ', p=' num2str(p)  ', C1=' num2str(C1) ', C2=' num2str(C2)  ', rho=' num2str(rho) ',alpha=' num2str(alpha) ', seg acc=' num2str(acc)]);
                        dlmwrite('latlrr_Yale10_WSNM.txt', [dim lambda p C1 C2 rho alpha acc NMI AR] , '-append', 'delimiter', '\t', 'newline', 'pc');
                    end
                end
            end
        end
    end
end
