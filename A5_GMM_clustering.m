function GMModel = A5_GMM_clustering(Spectral_Data, k, load_data)
    % Input: Spectral_Data
    % Output: Gaussian Mixture Model
    if load_data==1
        load('GMM_50ms_AllTrueRSUs_sharedcov_true_diagonal.mat','GMModel')
    else
        rng(0)

        Replicates = 5000;
        CovarianceType = 'diagonal';
        SharedCovariance = false;
        MaxIter = 700;
        

        %% GMM
        
        Results = struct();
        
        for k = 1:18

        GMModel = fitgmdist(Spectral_Data.spectral_Data,k,...
            'Replicates',Replicates,...
            'CovarianceType',CovarianceType,...
            'SharedCovariance',SharedCovariance,...
            'Options',...
            statset('MaxIter',MaxIter,...
            'UseParallel',true));
        
        Results.Model{k}=GMModel;
        Results.AIC{k}=GMModel.AIC;
        Results.BIC{k}=GMModel.BIC;
        Results.k{k}=k;
        
        Results.Rep{k}=Replicates;
        Results.Covtype{k}=CovarianceType;
        Results.MaxIte{k}=MaxIter;
        Results.Reg{k}=0;
        
        fprintf('Number of clusters in GMM:%d\n',k)
        
        end

        save('GMM_50ms_AllTrueRSUs_sharedcov_false_diagonal.mat','Results')
        
    end
end

