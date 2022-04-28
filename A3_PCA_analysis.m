function N_PCA = A3_PCA_analysis(Neural_Data_5T2P, load_data, plot_ind,...
                                 expert_ind)
    % Input: Neural Data
    % Output: Number of significant PCs
    
    if nargin < 4                          
        expert_ind = -1;
    end
    
    %% PCA
    
    Centered = false;
    
    switch expert_ind
        case -1
            temp_ND = Neural_Data_5T2P.pdf;   %% All mices
            temp_label = '';
        case 1
            temp_ND = Neural_Data_5T2P.pdf(Neural_Data_5T2P.expert~=0,:);
            temp_label = 'Expert';
        case 0
            temp_ND = Neural_Data_5T2P.pdf(Neural_Data_5T2P.expert==0,:);
            temp_label = 'Novice';
    end    

    for k = 1:size(temp_ND,1)
        temp_ND(k,:) = temp_ND(k,:) - mean(temp_ND(k,:));
    end
    
    [~, Score, Latent] = pca(temp_ND','Centered',Centered);   %% latent:principal component variances: the eigenvalues of the covariance matrix of temp_ND'.
    
    if plot_ind==1
        figure
        subplot(1,2,1)
        plot(cumsum(Latent)/sum(Latent))
        grid on
        xlim([0,250])
        title("Cumilative variance ratio")

        subplot(1,2,2)
        plot(Latent/sum(Latent))
        grid on
        title("Variance ratio")
        xlim([0,250])
        
        figure
        for i = 1:36
           subplot(6,6,i) 
           plot(Score(:,i))
           hold on
           grid on
           title(join(["Kernel",num2str(i)]))
           set (gca, 'XTick' , [], 'YTick' , [])
        end

        figure
        for i = 37:72
           subplot(6,6,i-36) 
           plot(Score(:,i))
           hold on
           grid on
           title(join(["Kernel",num2str(i)]))
        end

        figure
        for i = 73:80
           subplot(6,6,i-72) 
           plot(Score(:,i))
           hold on
           grid on
           title(join(["Kernel",num2str(i)]))
        end
    end
    
    if load_data == 1
        load(['PCA_Analysis_Data',temp_label,'.mat'],'Null_Sigma','Latent')
    else
        rng(0)
        Null_Sample = 1e3;                                 %% number of permutations
        Null_Sigma = zeros(Null_Sample,length(Latent));    %% 1000*125
        
        parfor j=1:Null_Sample                              
            temp = temp_ND;
            for i = 1:size(temp,1)
                p = randperm(size(temp,2));
                temp(i,:) = temp(i,p);
            end
            [~, ~, Latent_temp] = pca(temp','Centered',Centered);
            Null_Sigma(j,:) = Latent_temp;
        end
        
        save(['PCA_Analysis_Data',temp_label,'.mat'],'Null_Sigma','Latent')
    end
    
    
    y = mean(Null_Sigma,1);
    dy = std(Null_Sigma);

    if plot_ind==1
        figure
        plot(Latent,'LineWidth',2,'Color','#4B3F72')
        hold on
        plot(y,'LineWidth',2,'Color','#0072BD');
        plot([y-dy;y+dy]','LineStyle','--','Color','#0072BD');
        grid on
        legend('Variance ratio','Null variance ratio')
    end
    
    p_vals = zeros(size(Latent));    %% latent is the null hyp
    for i = 1:length(Latent)
        [~,p] = ttest(Null_Sigma(:,1),Latent(i));
        p_vals(i) = p;
    end

    ind = 1:length(Latent);
    N_PCA = max(ind((Latent> y(1))&(p_vals<(0.05/length(Latent)))));  %% bonferroni correction
end