function Neural_Data = data_2_pdf_cdf(Data,dT_base,dT_after,dt_res_pdf,dt_res_cdf,F_Name,verbose,j)

    % Given the data for a trial type, the function computes the firing
    % rate (pdf) and cumulative of all neurons
    
    N_Exp = length(Data.(F_Name));
    
    N_Neurons = 0;
    for Exp = 1:N_Exp
        N_Neurons = N_Neurons + length(Data.Area{Exp});  %% in the end we have number of all neurons in all 35 exp
    end    
    
    Neural_Data = struct();
    Neural_Data.MDS = zeros(1,N_Neurons);
    Neural_Data.clusterID = zeros(1,N_Neurons);
    Neural_Data.ClusterCounter = zeros(1,N_Neurons);
    Neural_Data.Area = cell(1,N_Neurons);
    Neural_Data.type = cell(1,N_Neurons);
    Neural_Data.width = zeros(1,N_Neurons);
    
    Neural_Data.NeuronNumber = zeros(1,N_Neurons);
    
    if j ==1 %%for baseline
        t_pdf_base = dt_res_pdf:dt_res_pdf:dT_base;
        t_pdf_base = length(t_pdf_base);
        Neural_Data.pdf = zeros(N_Neurons,t_pdf_base);  %% har radif yek neuron va har soton yek bin bedoone overlapp
        
        t_cdf_base = 0:dt_res_cdf:dT_base;                   %% cumolative khob qatan az sefr shoro mishe
        t_cdf_base = length(t_cdf_base);
        Neural_Data.cdf = zeros(N_Neurons,t_cdf_base);
    elseif j ==2 %%for after
        t_pdf_after = dt_res_pdf:dt_res_pdf:dT_after;
        t_pdf_after = length(t_pdf_after);
        Neural_Data.pdf = zeros(N_Neurons,t_pdf_after);
        
        t_cdf_after = 0:dt_res_cdf:dT_after;                   
        t_cdf_after = length(t_cdf_after);
        Neural_Data.cdf = zeros(N_Neurons,t_cdf_after);
        
    end
    
    
    Neural_Data.N = zeros(N_Neurons,1);        
    Neural_Data.N_trial = zeros(N_Neurons,1); 
    
    CR_whisking = 0;
    names = fieldnames(Data);
    for i = 1:length(names)
        if string(names{i})== "Whisking_AbsTimes_CR"
            CR_whisking = 1;
        end
    end
    
    FA_ITI = 0;
    names = fieldnames(Data);
    for i = 1:length(names)
        if string(names{i})== "LickOnsets_ITI_FA_LickTrace"
            FA_ITI = 1;
        end
    end
    
    k = 1;
    for Exp = 1:N_Exp
        if verbose==1
            disp(join([string(100*Exp/N_Exp), "%"],""))
        end
        N_Neuron_temp = length(Data.Area{Exp});   %% number of neuron in each exp
        
        for Neuron = 1:N_Neuron_temp
            if CR_whisking
                N_Trials = length(Data.Whisking_AbsTimes_CR{Exp});
            elseif FA_ITI
                N_Trials = length(Data.LickOnsets_ITI_FA_LickTrace{Exp});
            else
                 N_Trials = length(Data.TrialOnsets_All{Exp}); %% number of trials in each exp
            end
              
            Trials = 1:N_Trials;
            
%             if (~CR_whisking) && (~FA_ITI)
%                  Trials = Trials(Data.ReactionTimes{Exp}<=1);
%             end
           
           
            N_Trials = length(Trials);
            
            if N_Trials >=1
                temp = Data.(F_Name){Exp}{Neuron}{Trials(1)};
                for i = 2:N_Trials
                    temp = [temp, Data.(F_Name){Exp}{Neuron}{Trials(i)}];   %% yek list az zamane spike haye yek neuron masalan dar base dar ex1 va hame trial haye ex1
                end

                [CDF,PDF,N] =...
                    pdf_cdf_spikes(temp,dT_base,dT_after,dt_res_pdf,dt_res_cdf,j,0,0,0);

                PDF = PDF(:)';
                Neural_Data.pdf(k,:) = PDF./N_Trials;  %% average ro tedede trial haye yek exp

                CDF = CDF(:)';
                Neural_Data.cdf(k,:) = CDF;
                
                Neural_Data.MDS(k) = Data.MDS(Exp);
                Neural_Data.clusterID(k) = Data.clusterID{Exp}(Neuron);
                Neural_Data.ClusterCounter(k) = Data.ClusterCounter{Exp}(Neuron);
                Neural_Data.Area(k) = Data.Area{Exp}(Neuron);
                Neural_Data.width(k) = Data.width{Exp}(Neuron);
                Neural_Data.type(k) = Data.type{Exp}(Neuron);
                
                Neural_Data.NeuronNumber(k) = k; 
                Neural_Data.N(k) = N;              %% tedad kole spike haye yek neuron dar masalan hit dar masalan base
                Neural_Data.N_trial(k) = N_Trials;

                k = k + 1;
            else
                Neural_Data.pdf(k,:) = zeros(size(Neural_Data.pdf(k,:)));
                Neural_Data.cdf(k,:) = ones(size(Neural_Data.cdf(k,:)));
                
                Neural_Data.MDS(k) = Data.MDS(Exp);
                Neural_Data.clusterID(k) = Data.clusterID{Exp}(Neuron);
                Neural_Data.ClusterCounter(k) = Data.ClusterCounter{Exp}(Neuron);
                Neural_Data.Area(k) = Data.Area{Exp}(Neuron);
                Neural_Data.width(k) = Data.width{Exp}(Neuron);
                Neural_Data.type(k) = Data.type{Exp}(Neuron);
                
                Neural_Data.NeuronNumber(k) = k;
                Neural_Data.N(k) = 0;
                Neural_Data.N_trial(k) = 0;

                k = k + 1;
            end
        end
    end
end