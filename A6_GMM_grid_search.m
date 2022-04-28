function GMModel = A6_GMM_grid_search(Spectral_Data, load_data)

% Input: Spectral_Data
% Output: Best Gaussian Mixture Model

    if load_data==1
        load('GMM_50ms_AllTrueRSUs_Results.mat','GMModel')
    else
        
        rng(0)

        k = 5:30;
        nK = numel(k);
        Sigma = {'diagonal','full'};
        nSigma = numel(Sigma);
        SharedCovariance = {true,false};
        SCtext = {'true','false'};
        nSC = numel(SharedCovariance);
        RegularizationValue = 0;
        MaxIter = 800;
        options = statset('MaxIter',MaxIter,'UseParallel',true);
        Replicates = 5000;

        % Preallocation
        gm = cell(nK,nSigma,nSC);
        aic = zeros(nK,nSigma,nSC);
        bic = zeros(nK,nSigma,nSC);
        converged = false(nK,nSigma,nSC);
        
        Results = struct();
        Results.Model = cell(nK,nSigma,nSC);
        Results.AIC = zeros(nK,nSigma,nSC);
        Results.BIC = zeros(nK,nSigma,nSC);
        Results.SharedCov = cell(nK,nSigma,nSC);
        Results.Rep = zeros(nK,nSigma,nSC);
        Results.Covtype = cell(nK,nSigma,nSC);
        Results.MaxIte = zeros(nK,nSigma,nSC);
        Results.Reg = zeros(nK,nSigma,nSC);
        Results.k = zeros(nK,nSigma,nSC);
        
        X = Spectral_Data.spectral_Data;

        % Fit all models
        for m = 1:nSC
            for j = 1:nSigma
                for i = 1:nK
                    gm{i,j,m} = fitgmdist(X,k(i),...
                        'Replicates',Replicates,...
                        'CovarianceType',Sigma{j},...
                        'SharedCovariance',SharedCovariance{m},...
                        'RegularizationValue',RegularizationValue,...
                        'Options',options);

                    aic(i,j,m) = gm{i,j,m}.AIC;
                    bic(i,j,m) = gm{i,j,m}.BIC;
                    converged(i,j,m) = gm{i,j,m}.Converged;

                    Results.Model{i,j,m} = gm{i,j,m};
                    Results.AIC(i,j,m) = aic(i,j,m);
                    Results.BIC(i,j,m) = bic(i,j,m);
                    Results.SharedCov{i,j,m} = SCtext{m};
                    Results.Rep(i,j,m) = Replicates;
                    Results.Covtype{i,j,m} = Sigma{j};
                    Results.MaxIte(i,j,m) = MaxIter;
                    Results.Reg(i,j,m) = RegularizationValue;
                    Results.k(i,j,m) = i;

                    fprintf('CovarianceType:%s\nSharedCovariance:%s\nNumber of clusters in GMM:%d\n',Sigma{j},SCtext{m},k(i))

                end
            end
        end
        
        save('GMM_100ms_AllTrueRSUs_Expert_quality matrix_Histo_CoilGap.mat','Results')

        allConverge = (sum(converged(:)) == nK*nSigma*nSC);
        fprintf('All of GMM models Converged (binary):%d',allConverge)

        figure
        bar(reshape(bic,nK,nSigma*nSC))
        title('BIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex')
        xlabel('$c$','Interpreter','Latex')
        ylabel('BIC')
        legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
            'Full-unshared'})

        [~,q] = min(bic(:));
        GMModel = gm{q};
    end
end

