function [Neural_Data_5T2P,areas,T_pdf_len] = A2_data_preprocessing(Neural_Data_5T2P,F_names_Parts,F_names_Trials,t_pdf_base,t_pdf_after, load_data, verbose, plot_ind)
    % Input: The output of A1_firing_rate_computation
    % Output: preprocessed firing rates
    if load_data == 1
        load('SpikeData_5T2P_Neural_Normal_pdf_50ms_norm2_RSUs.mat','Neural_Data_5T2P','areas','T_pdf_len')
    else
        %% Removing baseline mean and which neuron is responsive in which part (baseline,..)
        if verbose == 1
            disp("Removing Baseline: Done")
        end
        N_Neurons = size(Neural_Data_5T2P.Hit_whisker.base_whisker.pdf,1);
        
        for i = 1:(length(F_names_Trials)-3)                                     % for whisker
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).Responsive = zeros(N_Neurons,2);         %% 2 parts
            for k = 1:N_Neurons
                baseline = mean(Neural_Data_5T2P.(trial_name).base_whisker.pdf(k,:));
                baseline_v = std(Neural_Data_5T2P.(trial_name).base_whisker.pdf(k,:));
                
                F_names_Parts{1} = 'base_whisker';
                F_names_Parts{2} = 'whisker';

                for j = 1:length(F_names_Parts)
                    part_name = F_names_Parts{j};
                    Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) = ...
                        Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) - baseline;

                    Neural_Data_5T2P.(trial_name).Responsive(k,j) = ...
                        max(abs( Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) )) > (2*baseline_v);      % WARNING: at least 2 std of baselin
                end
            end
        end
        
        for i = (length(F_names_Trials)-2)                                      % for Hit_jaw
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).Responsive = zeros(N_Neurons,2);         %% 2 parts
            for k = 1:N_Neurons
                
                baseline = mean(Neural_Data_5T2P.(trial_name).base_jaw.pdf(k,1:5));   %taking -1 to -0.5
                baseline_v = std(Neural_Data_5T2P.(trial_name).base_jaw.pdf(k,1:5));
                
                F_names_Parts{1} = 'base_jaw';
                F_names_Parts{2} = 'jaw';

                for j = 1:length(F_names_Parts)
                    part_name = F_names_Parts{j};
                    Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) = ...
                        Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) - baseline;

                    Neural_Data_5T2P.(trial_name).Responsive(k,j) = ...
                        max(abs( Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) )) > (2*baseline_v);      % WARNING: at least 2 std of baselin
                end
            end
        end
        
        for i = (length(F_names_Trials)-1)                                         % for FA_ITI
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).Responsive = zeros(N_Neurons,2);         %% 2 parts
            for k = 1:N_Neurons
                
                baseline = mean(Neural_Data_5T2P.(trial_name).base_FA_ITI.pdf(k,1:5));   %taking -1 to -0.5
                baseline_v = std(Neural_Data_5T2P.(trial_name).base_FA_ITI.pdf(k,1:5));
                
                F_names_Parts{1} = 'base_FA_ITI';
                F_names_Parts{2} = 'FA_ITI';

                for j = 1:length(F_names_Parts)
                    part_name = F_names_Parts{j};
                    Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) = ...
                        Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) - baseline;

                    Neural_Data_5T2P.(trial_name).Responsive(k,j) = ...
                        max(abs( Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) )) > (2*baseline_v);      % WARNING: at least 2 std of baselin
                end
            end
        end
        
        for i = length(F_names_Trials)                                             % for whisking
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).Responsive = zeros(N_Neurons,2);         %% 2 parts
            for k = 1:N_Neurons
                baseline = mean(Neural_Data_5T2P.(trial_name).base_whisking.pdf(k,1:5));   %taking -1 to -0.5
                baseline_v = std(Neural_Data_5T2P.(trial_name).base_whisking.pdf(k,1:5));
                
                F_names_Parts{1} = 'base_whisking';
                F_names_Parts{2} = 'whisking';

                for j = 1:length(F_names_Parts)
                    part_name = F_names_Parts{j};
                    Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) = ...
                        Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) - baseline;

                    Neural_Data_5T2P.(trial_name).Responsive(k,j) = ...
                        max(abs( Neural_Data_5T2P.(trial_name).(part_name).pdf(k,:) )) > (2*baseline_v);      % WARNING: at least 2 std of baselin
                end
            end
        end
        

        %% Combining parts (baseline and whisker/jaw/ITI/whisking) together
        if verbose == 1
            disp("Combining Time-Slots: Done")
        end
        T_pdf_len_base = length(t_pdf_base);
        T_pdf_len_after = length(t_pdf_after);

        N_Neurons = size(Neural_Data_5T2P.Hit_whisker.base_whisker.pdf,1);   
        for i = 1:(length(F_names_Trials)-3)
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).pdf = zeros(N_Neurons,T_pdf_len_base+T_pdf_len_after);  %% for combining 2 parts
            Neural_Data_5T2P.(trial_name).N = zeros(N_Neurons,2);

            Neural_Data_5T2P.(trial_name).pdf(:,0*T_pdf_len_base+(1:T_pdf_len_base)) = ...
                Neural_Data_5T2P.(trial_name).base_whisker.pdf;

            Neural_Data_5T2P.(trial_name).pdf(:,1*T_pdf_len_base+(1:T_pdf_len_after)) = ...
                Neural_Data_5T2P.(trial_name).whisker.pdf;
         
            Neural_Data_5T2P.(trial_name).N(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_whisker.N;
            Neural_Data_5T2P.(trial_name).N(:,2) = ...
                Neural_Data_5T2P.(trial_name).whisker.N;
            
            Neural_Data_5T2P.(trial_name).N_trial(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_whisker.N_trial;

        end
        
        for i = (length(F_names_Trials)-2)
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).pdf = zeros(N_Neurons,T_pdf_len_base+T_pdf_len_after);  %%  for combining 2 parts
            Neural_Data_5T2P.(trial_name).N = zeros(N_Neurons,2);

            Neural_Data_5T2P.(trial_name).pdf(:,0*T_pdf_len_base+(1:T_pdf_len_base)) = ...
                Neural_Data_5T2P.(trial_name).base_jaw.pdf;

            Neural_Data_5T2P.(trial_name).pdf(:,1*T_pdf_len_base+(1:T_pdf_len_after)) = ...
                Neural_Data_5T2P.(trial_name).jaw.pdf;
         
            Neural_Data_5T2P.(trial_name).N(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_jaw.N;
            Neural_Data_5T2P.(trial_name).N(:,2) = ...
                Neural_Data_5T2P.(trial_name).jaw.N;
            
            Neural_Data_5T2P.(trial_name).N_trial(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_jaw.N_trial;

        end
        
        for i = (length(F_names_Trials)-1)
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).pdf = zeros(N_Neurons,T_pdf_len_base+T_pdf_len_after);  %%  for combining 2 parts
            Neural_Data_5T2P.(trial_name).N = zeros(N_Neurons,2);

            Neural_Data_5T2P.(trial_name).pdf(:,0*T_pdf_len_base+(1:T_pdf_len_base)) = ...
                Neural_Data_5T2P.(trial_name).base_FA_ITI.pdf;

            Neural_Data_5T2P.(trial_name).pdf(:,1*T_pdf_len_base+(1:T_pdf_len_after)) = ...
                Neural_Data_5T2P.(trial_name).FA_ITI.pdf;
         
            Neural_Data_5T2P.(trial_name).N(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_FA_ITI.N;
            Neural_Data_5T2P.(trial_name).N(:,2) = ...
                Neural_Data_5T2P.(trial_name).FA_ITI.N;
            
            Neural_Data_5T2P.(trial_name).N_trial(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_FA_ITI.N_trial;

        end
        
        for i = length(F_names_Trials)
            trial_name = F_names_Trials{i};
            Neural_Data_5T2P.(trial_name).pdf = zeros(N_Neurons,T_pdf_len_base+T_pdf_len_after);  
            Neural_Data_5T2P.(trial_name).N = zeros(N_Neurons,2);

            Neural_Data_5T2P.(trial_name).pdf(:,0*T_pdf_len_base+(1:T_pdf_len_base)) = ...
                Neural_Data_5T2P.(trial_name).base_whisking.pdf;

            Neural_Data_5T2P.(trial_name).pdf(:,1*T_pdf_len_base+(1:T_pdf_len_after)) = ...
                Neural_Data_5T2P.(trial_name).whisking.pdf;
         
            Neural_Data_5T2P.(trial_name).N(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_whisking.N;
            Neural_Data_5T2P.(trial_name).N(:,2) = ...
                Neural_Data_5T2P.(trial_name).whisking.N;
            
            Neural_Data_5T2P.(trial_name).N_trial(:,1) = ...
                Neural_Data_5T2P.(trial_name).base_whisking.N_trial;

        end


        %% Making the general matrix   %% all 5 trials together
        if verbose == 1
            disp("Making General Matrix: Done")
        end
        T_pdf_len = size(Neural_Data_5T2P.Hit_whisker.pdf,2); %% inja 25->10+15
        N_Neurons = size(Neural_Data_5T2P.Hit_whisker.pdf,1);
        Neural_Data_5T2P.pdf = zeros(N_Neurons,5*(T_pdf_len_base+T_pdf_len_after));
        Neural_Data_5T2P.N = zeros(N_Neurons,5*2);
        Neural_Data_5T2P.N_trial = zeros(N_Neurons,5);
        Neural_Data_5T2P.ClusterCounter = zeros(N_Neurons,1);

        for i = 1:length(F_names_Trials)
            
            trial_name = F_names_Trials{i};
            
            Neural_Data_5T2P.pdf(:,(i-1)*T_pdf_len+(1:T_pdf_len)) = ...
                Neural_Data_5T2P.(trial_name).pdf;
            
            Neural_Data_5T2P.N(:,(i-1)*2+(1:2)) = ... 
                Neural_Data_5T2P.(trial_name).N;

            Neural_Data_5T2P.N_trial(:,i) = ... 
                Neural_Data_5T2P.(trial_name).N_trial;
        end
        
        Neural_Data_5T2P.ClusterCounter(:,1) = ... 
                (Neural_Data_5T2P.Hit_whisker.base_whisker.ClusterCounter)';

        temp_name = fieldnames(Neural_Data_5T2P.Hit_whisker.whisker);
        temp_name = temp_name([1:2 4:(end-4)]);     %% except ClusterCounter, pdf, cdf, N, N_trial

        for i = 1:length(temp_name)
            f_name = temp_name{i};
            Neural_Data_5T2P.(f_name) = Neural_Data_5T2P.Hit_whisker.whisker.(f_name);
        end

        %% Normalizing firing rates
        if verbose == 1
            disp("Normalization: Done")
        end
        N_Neurons = size(Neural_Data_5T2P.pdf,1);
        for k = 1:N_Neurons
            scale = max(Neural_Data_5T2P.pdf(k,:)) - ...
                    min(Neural_Data_5T2P.pdf(k,:));
            if scale ~= 0
                Neural_Data_5T2P.pdf(k,:) = Neural_Data_5T2P.pdf(k,:)/scale;
            end
        end

        %% Removing "bad" neurons
        if verbose == 1
            disp("Removing Neurons: Done")
        end
        N_Tot = sum(Neural_Data_5T2P.N,2);                
        Cond1 = N_Tot>=1;                              
        Cond2 = (min(Neural_Data_5T2P.N_trial,[],2)>0);  %%more than 0 spikes throughout the recording, and with more than 1 trials for each trial-type Hit, Miss..   

        Ind_2P = Cond1 & Cond2;

        Ind_2P = Ind_2P & (Neural_Data_5T2P.type=="RSU")';                            % type == RSU
        
        if verbose == 1
            disp('Number of total neurons:')
            disp(N_Neurons)
            disp('Number of removed neurons:')
            disp(N_Neurons - sum(Ind_2P))
        end

        Neural_Data_5T2P.pdf = Neural_Data_5T2P.pdf(Ind_2P,:);
        Neural_Data_5T2P.N = Neural_Data_5T2P.N(Ind_2P,:);
        Neural_Data_5T2P.N_trial = Neural_Data_5T2P.N_trial(Ind_2P,:);
        Neural_Data_5T2P.ClusterCounter = Neural_Data_5T2P.ClusterCounter(Ind_2P,:);

        temp_name = fieldnames(Neural_Data_5T2P);
        temp_name = temp_name(10:end);  

        for i = 1:length(temp_name)
            f_name = temp_name{i};
            Neural_Data_5T2P.(f_name) = Neural_Data_5T2P.(f_name)(Ind_2P);
        end

        %% On the average PDF plots (grand average PSTH)
        areas = unique(Neural_Data_5T2P.Area);
        %t_pdf_new = (0:2:(1.99*T_pdf_len-1))/T_pdf_len - 1;  %% from -1 to 0.98 (200 bins)
        if plot_ind == 1
            for Area_ind = 1:length(areas)
                Area = string(areas{Area_ind});
                figure

                Areas_trial = Neural_Data_5T2P.Area;
                %Experts_trial = Neural_Data_5T2P.expert;
                Ind = 1:length(Areas_trial);

                y = mean(Neural_Data_5T2P.pdf(Ind((Areas_trial==Area)),:),1);

                dy = std(Neural_Data_5T2P.pdf(Ind((Areas_trial==Area)),:),1)...
                    /sqrt(sum((Areas_trial==Area)));
                
                p1 = plot(y,'LineWidth',2,'Color','#0072BD');
                hold on
                plot([1:length(y);1:length(y)]',[y-dy;y+dy]','LineStyle','--','Color','#0072BD');
                
                
                plot([10,10],[-1,1],'--red');
                plot([25,25],[-1,1],'--black');
                plot([35,35],[-1,1],'--red');
                plot([50,50],[-1,1],'--black');
                plot([60,60],[-1,1],'--red');
                plot([75,75],[-1,1],'--black');
                plot([85,85],[-1,1],'--red');
                plot([100,100],[-1,1],'--black');
                plot([110,110],[-1,1],'--red');
                
                ax = gca;
                ax.XTick = [0,0.4,1.4,2.4,3.4,4.4,5]*T_pdf_len;
                ax.XTickLabel = {'-1','Hit.Stim.Aligned','Miss.Stim.Aligned','Hit.JawOnset.Aligned','Spontaneous.Licks.Piezo.Aligned','SpontaneousWhiskingAligned','1.5'};
                ax.FontSize = 12;
                xlim([0,125])
    
                ax.XLabel.String = sprintf('\n%s', 'All trials duration: -1 to 1.5s');

                grid on
                ylim([-0.2,0.6])
                title(join(["Area:",Area]))
                legend(p1,"Expert")
                fig = gcf;
                set(fig,'WindowState','fullscreen');
                saveas(fig,join([Area,'_pdf.png'],""))
                close(fig)
                pause(0.5)
            end
        end

        %% Sorting based on areas
        [Neural_Data_5T2P.Area, ind_temp] = sort(Neural_Data_5T2P.Area);
        Neural_Data_5T2P.pdf = Neural_Data_5T2P.pdf(ind_temp,:);
        Neural_Data_5T2P.N = Neural_Data_5T2P.N(ind_temp,:);
        Neural_Data_5T2P.N_trial = Neural_Data_5T2P.N_trial(ind_temp,:);
        Neural_Data_5T2P.ClusterCounter = Neural_Data_5T2P.ClusterCounter(ind_temp,:);

        temp_name = fieldnames(Neural_Data_5T2P);
        temp_name = temp_name(10:end);

        for i = 1:length(temp_name)
            f_name = temp_name{i};
            if f_name ~= "Area"
                Neural_Data_5T2P.(f_name) = Neural_Data_5T2P.(f_name)(ind_temp);
            end
        end

        %% HeatMap of firing rate
        N_Neurons = size(Neural_Data_5T2P.pdf,1);
        areas = unique(Neural_Data_5T2P.Area);
        Ind = 1:N_Neurons;
        Ind2 = zeros(size(areas));
        if plot_ind == 1
            figure
            imagesc(Neural_Data_5T2P.pdf)
            hold on
            %%khat vasate heatmap
            line([1,1]*(T_pdf_len_base/(T_pdf_len_base+T_pdf_len_after))*(T_pdf_len)+0.5, [0 N_Neurons], 'Color', 'w', 'LineStyle', '--')
            plot([1,1]*(T_pdf_len)+0.5, [0,N_Neurons], 'black','LineWidth',0.5)   
            line([1,1]*(1+(T_pdf_len_base/(T_pdf_len_base+T_pdf_len_after)))*(T_pdf_len)+0.5, [0 N_Neurons], 'Color', 'w', 'LineStyle', '--')
            plot([1,1]*2*T_pdf_len+0.5, [0,N_Neurons], 'black','LineWidth',0.5)
            line([1,1]*(2+(T_pdf_len_base/(T_pdf_len_base+T_pdf_len_after)))*(T_pdf_len)+0.5, [0 N_Neurons], 'Color', 'w', 'LineStyle', '--')
            plot([1,1]*3*T_pdf_len+0.5, [0,N_Neurons], 'black','LineWidth',0.5)
            line([1,1]*(3+(T_pdf_len_base/(T_pdf_len_base+T_pdf_len_after)))*(T_pdf_len)+0.5, [0 N_Neurons], 'Color', 'w', 'LineStyle', '--')
            plot([1,1]*4*T_pdf_len+0.5, [0,N_Neurons], 'black','LineWidth',0.5)
            line([1,1]*(4+(T_pdf_len_base/(T_pdf_len_base+T_pdf_len_after)))*(T_pdf_len)+0.5, [0 N_Neurons], 'Color', 'w', 'LineStyle', '--')
            plot([1,1]*5*T_pdf_len+0.5, [0,N_Neurons], 'black','LineWidth',0.5)
           
            for i = 1:length(areas)
                Area = string(areas{i});
                temp = max(Ind(Neural_Data_5T2P.Area==Area));
                Ind2(i) = round(median(Ind(Neural_Data_5T2P.Area==Area)));
                plot([0,5*T_pdf_len+1], [1,1]*temp, 'black','LineWidth',0.5)
            end

            ax = gca;
            ax.YTick = Ind2;
            ax.YTickLabel = areas;
            ax.XTick = [0.4,1.4,2.4,3.4,4.4]*T_pdf_len;
            ax.XTickLabel = {'Hit.Stim.Aligned','Miss.Stim.Aligned','Hit.JawOnset.Al','Spon.Licks.Piezo.Al','Spon.Whisking.Al'};
            ax.FontSize = 12.5;

                ax.XLabel.String = sprintf('\n%s', 'All trials duration: -1 to 1.5s');

            title('Area Heatmap. All 3 Stimuli combined (StimAmp 1 excluded), Hit and Miss balanced in number of trials')
            pbaspect ([1 1 1])
            colorbar

            fig = gcf;
            set(fig,'WindowState','fullscreen');
            saveas(fig,join(['area_HeatMap.png'],""))
        end
        
        %% Save Data
        save('SpikeData_5T2P_Neural_Normal_pdf_100ms_norm2_RSUs.mat','Neural_Data_5T2P','areas','T_pdf_len')

    end
end