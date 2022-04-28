function Spectral_Data = A4_spectral_embedding(Neural_Data_5T2P_preprocessed, N_PCA, load_data, plot_ind, expert_ind)
    % Input: Neural Data and N_PCA
    % Output: Trasnformed version of data in spectral space
    
    if nargin < 5
        expert_ind = -1;
    end
    
    switch expert_ind
        case -1
            temp_ND_pdf = Neural_Data_5T2P_preprocessed.pdf;
            temp_ND_area = Neural_Data_5T2P_preprocessed.Area;
            temp_label = '';
            sig = 0.135;      
            K_Spect = 1+8;   %% From elbow method
        case 1
            temp_ND_pdf = Neural_Data_5T2P_preprocessed.pdf(Neural_Data_5T2P_preprocessed.expert~=0,:);
            temp_ND_area = Neural_Data_5T2P_preprocessed.Area(Neural_Data_5T2P_preprocessed.expert~=0);
            temp_label = 'Expert';
            sig = 0.115;
            K_Spect = 1+13;
        case 0
            temp_ND_pdf = Neural_Data_5T2P_preprocessed.pdf(Neural_Data_5T2P_preprocessed.expert==0,:);
            temp_ND_area = Neural_Data_5T2P_preprocessed.Area(Neural_Data_5T2P_preprocessed.expert==0);
            temp_label = 'Novice';
            sig = 0.165;
            K_Spect = 1+10;
    end  
   
    
    if load_data == 1
        load(['spectral_Data_smaller', temp_label, '.mat'], "Spectral_Data")
    else
        %% PCA
        Centered = false;
        for k = 1:size(temp_ND_pdf,1)
            temp_ND_pdf(k,:) = temp_ND_pdf(k,:) - mean(temp_ND_pdf(k,:));
        end
        
        [Coeff, ~, ~] = pca(temp_ND_pdf','Centered',Centered);
        PCA_Data = Coeff(:,1:N_PCA);
        
        D_Both = squareform(pdist(PCA_Data,'euclidean'));   
        
        S_Both = exp(-D_Both/sig);
        
        average = mean(S_Both(:));
        
        [Lambda_Both,V_Both,flag] = spect_clust(S_Both,-1);
        disp("Spectraling embedding flag:")
        disp(flag)

        spectral_Data = V_Both(:,2:K_Spect); %% after excluding the very first
        
        Spectral_Data = struct();
        Spectral_Data.spectral_Data =spectral_Data;
        Spectral_Data.ClusterCounter = Neural_Data_5T2P_preprocessed.ClusterCounter;
        
        save(['spectral_Data', temp_label, '.mat'], "Spectral_Data", "K_Spect", "Lambda_Both","V_Both","sig")
        save(['spectral_Data_smaller', temp_label, '.mat'], "Spectral_Data","sig")
        
        if plot_ind==1
            figure
            histogram(S_Both(:))
            grid on
            xlim([0,1])
            title('Similarity Histogram')
            
            N_Neurons = size(temp_ND_pdf,1);
            areas = unique(temp_ND_area);
            Ind = 1:N_Neurons;
            Ind2 = zeros(size(areas));

            figure
            imagesc(S_Both)

            hold on
            for i = 1:length(areas)
                Area = string(areas{i});
                temp = max(Ind(temp_ND_area==Area));
                Ind2(i) = round(median(Ind(temp_ND_area==Area)));
                plot([0,N_Neurons+1], [1,1]*temp, 'black','LineWidth',0.5)
                plot([1,1]*temp, [0,N_Neurons+1], 'black','LineWidth',0.5)
            end

            ax = gca;
            ax.YTick = Ind2;
            ax.YTickLabel = areas;
            ax.XTick = Ind2;
            ax.XTickLabel = areas;
            xtickangle(90)

            colorbar
            title("Similarity heatmap")
            
            figure
            plot(Lambda_Both(2:end),'*')
            hold on
            plot(Lambda_Both(2:end))
            axis([1,50,0.96,1])
            grid on
            title("Eigenvalues")
            
            figure
            plot3(V_Both(:,2), V_Both(:,3), V_Both(:,4), '.')
            grid on
            title("Neurons in the 1st 3 dimensions of the spectral space")
            pbaspect ([1 1 1])
        end
    end
end