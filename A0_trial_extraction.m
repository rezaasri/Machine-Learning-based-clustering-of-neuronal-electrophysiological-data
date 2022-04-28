function Data_5T2P = A0_trial_extraction(spikeData, load_data)
% Input: spikeData
% Output: A similar struct whose unnecessary fields are removed and the
% spike times of its trials are extracted
% If load_data = 1, then function does not do anything and just load
% the prviously processed data

if load_data == 1
    load('SpikeData_5T2P.mat','Data_5T2P')
else
    %% Selecting sub feilds
    Data = struct();
    
    Data.MDS = spikeData.MDS;
    Data.TrialOnsets_All = spikeData.TrialOnsets_All;
    Data.StimIndices = spikeData.StimIndices;
    Data.StimAmps = spikeData.StimAmps;
    Data.HitIndices = spikeData.HitIndices;
    Data.MissIndices = spikeData.MissIndices;
    Data.FAIndices = spikeData.FAIndices;
    Data.CRIndices = spikeData.CRIndices;
    Data.ReactionTimes = spikeData.ReactionTimes;
    Data.spikets = spikeData.spikets;
    Data.JawOnsetsTms = spikeData.JawOnsetsTms;
    
    %taking just the start of whisking time
    for Exp = 1:length(Data.spikets)
        mat = (spikeData.Whisking_AbsTimes_CR{Exp});
        mat = mat(:,1)';
        Data.Whisking_AbsTimes_CR{Exp} = mat;
    end
    
    %Data.SpontaneousLicks = spikeData.SpontaneousLicks; %no more exist
    Data.LickOnsets_ITI_FA_LickTrace = spikeData.SpontLicksPiezo_LickOn; %% Ana said use this, the last one had baseline problems!
    Data.clusterID = spikeData.clusterID;
    Data.ClusterCounter = spikeData.ClusterCounter;
    Data.Area = spikeData.Area;
    Data.width = spikeData.width;
    Data.type = spikeData.type;
    
    
    %% Trial-specific spike-timings part by part
    Data.spikets_trials_base_whisker = cell(size(Data.spikets));    % Baseline of whisker stim
    Data.spikets_trials_whisker = cell(size(Data.spikets));         % After Whisker stim
    
    Data.spikets_trials_base_jaw = cell(size(Data.spikets));   % Baseline of jawonset
    Data.spikets_trials_jaw = cell(size(Data.spikets));        % after jawonset
    
    Data.spikets_trials_base_whisking = cell(size(Data.spikets));   % Baseline of whisking onset
    Data.spikets_trials_whisking = cell(size(Data.spikets));        % after whisking onset
    
    %         h = 0; % counter for modulated trials in terms of producing spikes in 5ms gap
    %         N_neuron_times_N_trials = 0;
    for Exp = 1:length(Data.spikets)
        
        Data.spikets_trials_base_whisker{Exp} = cell(size(Data.spikets{Exp}));
        Data.spikets_trials_whisker{Exp} = cell(size(Data.spikets{Exp}));
        
        Data.spikets_trials_base_jaw{Exp} = cell(size(Data.spikets{Exp}));
        Data.spikets_trials_jaw{Exp} = cell(size(Data.spikets{Exp}));
        
        Data.spikets_trials_base_FA_ITI{Exp} = cell(size(Data.spikets{Exp}));
        Data.spikets_trials_FA_ITI{Exp} = cell(size(Data.spikets{Exp}));
        
        Data.spikets_trials_base_whisking{Exp} = cell(size(Data.spikets{Exp}));
        Data.spikets_trials_whisking{Exp} = cell(size(Data.spikets{Exp}));
        
        for Neuron = 1:length(Data.spikets{Exp})
            Data.spikets_trials_base_whisker{Exp}{Neuron} = ...
                cell(size(Data.TrialOnsets_All{Exp}));
            
            Data.spikets_trials_whisker{Exp}{Neuron} = ...
                cell(size(Data.TrialOnsets_All{Exp}));
            
            Data.spikets_trials_base_jaw{Exp}{Neuron} = ...
                cell(size(Data.JawOnsetsTms{Exp}));
            
            Data.spikets_trials_jaw{Exp}{Neuron} = ...
                cell(size(Data.JawOnsetsTms{Exp}));
            
            Data.spikets_trials_base_FA_ITI{Exp}{Neuron} = ...
                cell(size(Data.LickOnsets_ITI_FA_LickTrace{Exp}));
            
            Data.spikets_trials_FA_ITI{Exp}{Neuron} = ...
                cell(size(Data.LickOnsets_ITI_FA_LickTrace{Exp}));
            
            
            Data.spikets_trials_base_whisking{Exp}{Neuron} = ...
                cell(size(Data.Whisking_AbsTimes_CR{Exp}));
            
            Data.spikets_trials_whisking{Exp}{Neuron} = ...
                cell(size(Data.Whisking_AbsTimes_CR{Exp}));
            
            %                 area = Data.Area{Exp}{Neuron};
            
            for Trial = 1:length(Data.TrialOnsets_All{Exp})
                T0 = Data.TrialOnsets_All{Exp}(Trial);
                
                T_sta = T0-1;
                T_end = T0;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_base_whisker{Exp}{Neuron}{Trial} = Spikes;
                
                %                     if area == "mPFC"
                %                         n_spikes_in_5ms = floor(length(Spikes)*0.07); %% numbers are the probability of spiking in 5ms
                %                     elseif area == "tjM1"
                %                         n_spikes_in_5ms = floor(length(Spikes)*0.045);
                %                     elseif area == "wS1"
                %                         n_spikes_in_5ms = floor(length(Spikes)*0.03);
                %                     end
                %
                %                     if n_spikes_in_5ms ~=0
                %                        h = h+1;
                %                     end
                %
                T_sta = T0;
                T_end = T0+1.5;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                %
                %                     if n_spikes_in_5ms ~= 0
                %                         Spikes = cat(2,sort(randsample(linspace(0,0.005,50),n_spikes_in_5ms,true)),Spikes); %% adding random spikes in 5ms coil artifact
                %                     end
                %
                Data.spikets_trials_whisker{Exp}{Neuron}{Trial} = Spikes;
                
            end
            
            for Trial = 1:length(Data.JawOnsetsTms{Exp})
                T1 = Data.JawOnsetsTms{Exp}(Trial);
                
                T_sta = T1-1;
                T_end = T1;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_base_jaw{Exp}{Neuron}{Trial} = Spikes;
                
                T_sta = T1;
                T_end = T1+1.5;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_jaw{Exp}{Neuron}{Trial} = Spikes;
            end
            
            for Trial = 1:length(Data.LickOnsets_ITI_FA_LickTrace{Exp})
                T2 = Data.LickOnsets_ITI_FA_LickTrace{Exp}(Trial);
                
                T_sta = T2-1;
                T_end = T2;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_base_FA_ITI{Exp}{Neuron}{Trial} = Spikes;
                
                T_sta = T2;
                T_end = T2+1.5;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_FA_ITI{Exp}{Neuron}{Trial} = Spikes;
            end
            
            for Trial = 1:length(Data.Whisking_AbsTimes_CR{Exp})
                T2 = Data.Whisking_AbsTimes_CR{Exp}(Trial);
                
                T_sta = T2-1;
                T_end = T2;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_base_whisking{Exp}{Neuron}{Trial} = Spikes;
                
                T_sta = T2;
                T_end = T2+1.5;
                Ind = (Data.spikets{Exp}{Neuron} >= T_sta) &...
                    (Data.spikets{Exp}{Neuron} <  T_end);
                Spikes = Data.spikets{Exp}{Neuron}(Ind) - T_sta;
                Data.spikets_trials_whisking{Exp}{Neuron}{Trial} = Spikes;
            end
            
        end
        
        %             N_neuron_times_N_trials = N_neuron_times_N_trials + length(Data.TrialOnsets_All{Exp})*length(Data.spikets{Exp});
        
    end
    
    %         fprintf('N_neuron_times_N_trials:%d\nN_neuron_times_N_trials modulated in 5ms coil artefact:%d\nFraction:%f%%\n',N_neuron_times_N_trials,h,h/N_neuron_times_N_trials*100)
    
    Data.spikets = [];
    Data = rmfield(Data,'spikets');
    
    %% 6 kind of trials
    Data_5T2P = struct();
    
    Data_5T2P.Hit_whisker = select_data_trials(Data,'Hit_whisker');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_base_jaw');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_jaw');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_base_whisking');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_whisking');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_base_FA_ITI');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'spikets_trials_FA_ITI');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'Whisking_AbsTimes_CR');
    Data_5T2P.Hit_whisker = rmfield(Data_5T2P.Hit_whisker,'LickOnsets_ITI_FA_LickTrace');
    
    Data_5T2P.Miss_whisker = select_data_trials(Data,'Miss_whisker');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_base_jaw');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_jaw');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_base_whisking');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_whisking');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_base_FA_ITI');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'spikets_trials_FA_ITI');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'Whisking_AbsTimes_CR');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'JawOnsetsTms');
    Data_5T2P.Miss_whisker = rmfield(Data_5T2P.Miss_whisker,'LickOnsets_ITI_FA_LickTrace');
    
    Data_5T2P.Hit_jaw = select_data_trials(Data,'Hit_jaw');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_base_whisker');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_whisker');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_base_whisking');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_whisking');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_base_FA_ITI');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'spikets_trials_FA_ITI');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'Whisking_AbsTimes_CR');
    Data_5T2P.Hit_jaw = rmfield(Data_5T2P.Hit_jaw,'LickOnsets_ITI_FA_LickTrace');
    
    Data_5T2P.FA_ITI.MDS=Data_5T2P.Hit_whisker.MDS;
    Data_5T2P.FA_ITI.clusterID=Data_5T2P.Hit_whisker.clusterID;
    Data_5T2P.FA_ITI.ClusterCounter=Data_5T2P.Hit_whisker.ClusterCounter;
    Data_5T2P.FA_ITI.Area=Data_5T2P.Hit_whisker.Area;
    Data_5T2P.FA_ITI.width=Data_5T2P.Hit_whisker.width;
    Data_5T2P.FA_ITI.type=Data_5T2P.Hit_whisker.type;
    Data_5T2P.FA_ITI.LickOnsets_ITI_FA_LickTrace=Data.LickOnsets_ITI_FA_LickTrace;
    Data_5T2P.FA_ITI.spikets_trials_base_FA_ITI=Data.spikets_trials_base_FA_ITI;
    Data_5T2P.FA_ITI.spikets_trials_FA_ITI=Data.spikets_trials_FA_ITI;
    
    Data_5T2P.CR_whisking.MDS=Data_5T2P.Hit_whisker.MDS;
    Data_5T2P.CR_whisking.clusterID=Data_5T2P.Hit_whisker.clusterID;
    Data_5T2P.CR_whisking.ClusterCounter=Data_5T2P.Hit_whisker.ClusterCounter;
    Data_5T2P.CR_whisking.Area=Data_5T2P.Hit_whisker.Area;
    Data_5T2P.CR_whisking.width=Data_5T2P.Hit_whisker.width;
    Data_5T2P.CR_whisking.type=Data_5T2P.Hit_whisker.type;
    Data_5T2P.CR_whisking.Whisking_AbsTimes_CR=Data.Whisking_AbsTimes_CR;
    Data_5T2P.CR_whisking.spikets_trials_base_whisking=Data.spikets_trials_base_whisking;
    Data_5T2P.CR_whisking.spikets_trials_whisking=Data.spikets_trials_whisking;
    %Data_5T2P.CR_whisking = select_data_trials(Data,'CR_whisking');
    
    %% Save Data
    save('SpikeData_5T2P.mat','Data_5T2P')
end
end