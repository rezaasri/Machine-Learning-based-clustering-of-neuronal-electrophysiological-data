function [Neural_Data_5T2P,dt_res_pdf,F_names_Parts,F_names_Trials,t_pdf_base,t_pdf_after] = A1_firing_rate_computation(Data_5T2P, load_data, verbose)
    % Input: The output of A0_trial_extraction
    % Output: A structure containing the firing rate (pdf) and cumulative
    % firing rate (cdf) of each neuron for each trial-type
    
    if load_data == 1
        load('SpikeData_5T2P_Neural_50ms.mat','Neural_Data_5T2P','dt_res_pdf','F_names_Parts','F_names_Trials','t_pdf_base','t_pdf_after')
    else
        %% Setting specification
        dt_res_pdf = 0.1;                     %% bin size
        dt_res_cdf = 0.1;
        
        dT_base = 1;
        dT_after = 1.5;
        
        t_cdf_base = 0:dt_res_cdf:dT_base;
        t_cdf_after = 0:dt_res_cdf:dT_after;
        
        t_pdf_base = dt_res_pdf:dt_res_pdf:dT_base;
        t_pdf_after = dt_res_pdf:dt_res_pdf:dT_after;

        
        %% Firing rate (pdf) and cumulative firing rate (cdf) computation
        Neural_Data_5T2P = struct();
        
        F_names_Trials = fieldnames(Data_5T2P);

        F_names_Parts = cell(2,1);
        F_names_Parts{1} = 'base_whisker';
        F_names_Parts{2} = 'whisker';

        Trial_Types = fieldnames(Data_5T2P.Hit_whisker);
        Trial_Types = Trial_Types(end-1:end);       %% spikets_trials...
        
        for i = 1:(length(F_names_Trials)-3)       %% for whisker trials
            trial_name = F_names_Trials{i};
            for j = 1:length(F_names_Parts)
                if verbose==1
                    part_name = F_names_Parts{j};
                    disp('***********************************')
                    disp(trial_name)
                    disp('***********************************')
                    disp(part_name)
                    disp('***********************************')
                end
                Neural_Data_5T2P.(trial_name).(part_name) = data_2_pdf_cdf(...
                                                                Data_5T2P.(trial_name),dT_base,dT_after,dt_res_pdf,dt_res_cdf,...
                                                                Trial_Types{j},verbose,j);
            end
        end
        
        F_names_Parts{1} = 'base_jaw';
        F_names_Parts{2} = 'jaw';
        
        Trial_Types = fieldnames(Data_5T2P.Hit_jaw);
        Trial_Types = Trial_Types(end-1:end);
        
        for i = (length(F_names_Trials)-2)     %% for Hit_jaw trial
            trial_name = F_names_Trials{i};
            for j = 1:length(F_names_Parts)
                if verbose==1
                    part_name = F_names_Parts{j};
                    disp('***********************************')
                    disp(trial_name)
                    disp('***********************************')
                    disp(part_name)
                    disp('***********************************')
                end
                Neural_Data_5T2P.(trial_name).(part_name) = data_2_pdf_cdf(...
                                                                Data_5T2P.(trial_name),dT_base,dT_after,dt_res_pdf,dt_res_cdf,...
                                                                Trial_Types{j},verbose,j);
            end
        end
        
         F_names_Parts{1} = 'base_FA_ITI';
        F_names_Parts{2} = 'FA_ITI';
        
        Trial_Types = fieldnames(Data_5T2P.FA_ITI);  
        Trial_Types = Trial_Types(end-1:end);
        
        for i = (length(F_names_Trials)-1)     %%for FA_ITI trial
            trial_name = F_names_Trials{i};
            for j = 1:length(F_names_Parts)
                if verbose==1
                    part_name = F_names_Parts{j};
                    disp('***********************************')
                    disp(trial_name)
                    disp('***********************************')
                    disp(part_name)
                    disp('***********************************')
                end
                Neural_Data_5T2P.(trial_name).(part_name) = data_2_pdf_cdf(...
                                                                Data_5T2P.(trial_name),dT_base,dT_after,dt_res_pdf,dt_res_cdf,...
                                                                Trial_Types{j},verbose,j);
            end
        end
        
        F_names_Parts{1} = 'base_whisking';
        F_names_Parts{2} = 'whisking';
        
        Trial_Types = fieldnames(Data_5T2P.CR_whisking);
        Trial_Types = Trial_Types(end-1:end);
        
        for i = length(F_names_Trials)     %% for the only CR_whisking trial
            trial_name = F_names_Trials{i};
            for j = 1:length(F_names_Parts)
                if verbose==1
                    part_name = F_names_Parts{j};
                    disp('***********************************')
                    disp(trial_name)
                    disp('***********************************')
                    disp(part_name)
                    disp('***********************************')
                end
                Neural_Data_5T2P.(trial_name).(part_name) = data_2_pdf_cdf(...
                                                                Data_5T2P.(trial_name),dT_base,dT_after,dt_res_pdf,dt_res_cdf,...
                                                                Trial_Types{j},verbose,j);
            end
        end


        %% Save Neural_Data_5T2P
        save('SpikeData_5T2P_Neural_100ms.mat','Neural_Data_5T2P',...
                                    'dt_res_cdf','dt_res_pdf','F_names_Parts','F_names_Trials',...
                                    't_cdf_base','t_cdf_after','t_pdf_base','t_pdf_after')
    end

end