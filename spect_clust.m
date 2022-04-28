function [Lambda,V,flag] = spect_clust(S,K)
    DinvS = diag(sum(S,2).^(-0.5));   
    L = DinvS * S * DinvS;
    clear DinvS
    clear S
    L = eye(size(L)) - L;
    
    if K==-1
        [V,Lambda] = eig(L);  %% V is eigenvector and diagonal matrix Lambda of generalized eigenvalues
        flag = 'All';
    else
        [V,Lambda,flag] = eigs(L, K, 'smallestabs');
    end
    clear L
    
    
    Lambda = diag(Lambda);
    
    s_lambda = sort(Lambda);
    
    plot(diff(s_lambda(2:40)));
    
end