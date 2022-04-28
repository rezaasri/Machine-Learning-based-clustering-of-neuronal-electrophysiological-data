function Data_selected = select_data_trials(Data,Type)

names = fieldnames(Data);


switch Type
    case 'Hit_whisker'
        Data_selected = Data;
        for Exp = 1:length(Data.(names{1}))
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==4); %% balancing
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==4);
            if a<b
                c = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==4);
                Ind1 = sort(c(randperm(size(c,1), a)));
            else
                Ind1 = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==4);
            end
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            if a<b
                c = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==3);
                Ind2 = sort(c(randperm(size(c,1), a)));
            else
                Ind2 = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            end
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            if a<b
                c = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==2);
                Ind3 = sort(c(randperm(size(c,1), a)));
            else
                Ind3 = find(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            end
            
            Ind = sort(cat(1,Ind1,Ind2,Ind3));
            
            for i = 1:length(names)
                if (string(names{i}) == "spikets_trials_base_whisker")||...
                        (string(names{i}) == "spikets_trials_whisker")
                    
                    for Neuron = 1:length(Data.(names{i}){Exp})
                        Data_selected.(names{i}){Exp}{Neuron} =...
                            Data.(names{i}){Exp}{Neuron}(Ind);
                    end
                    
                elseif (string(names{i}) ~= "clusterID")&&...
                        (string(names{i}) ~= "ClusterCounter")&&...
                        (string(names{i}) ~= "Area")&&...
                        (string(names{i}) ~= "width")&&...
                        (string(names{i}) ~= "type")&&...
                        (string(names{i}) ~= "MDS")&&...
                        (string(names{i}) ~= "spikets_trials_base_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_base_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_base_FA_ITI")&&...
                        (string(names{i}) ~= "spikets_trials_FA_ITI")&&...
                        (string(names{i}) ~= "Whisking_AbsTimes_CR")&&...
                        (string(names{i}) ~= "LickOnsets_ITI_FA_LickTrace")
                    
                    Data_selected.(names{i}){Exp} = Data.(names{i}){Exp}(Ind);
                    
                end
            end
        end
        
    case 'Miss_whisker'
        Data_selected = Data;
        for Exp = 1:length(Data.(names{1}))
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==4); %% balancing
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==4);
            if a>b
                c = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==4);
                Ind1 = sort(c(randperm(size(c,1), b)));
            else
                Ind1 = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==4);
            end
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            if a>b
                c = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==3);
                Ind2 = sort(c(randperm(size(c,1), b)));
            else
                Ind2 = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==3);
            end
            
            a = sum(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            b = sum(Data.HitIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            if a>b
                c = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==2);
                Ind3 = sort(c(randperm(size(c,1), b)));
            else
                Ind3 = find(Data.MissIndices{Exp}==1 & Data.StimAmps{Exp}==2);
            end
            
            Ind = sort(cat(1,Ind1,Ind2,Ind3));
            
            for i = 1:length(names)
                if (string(names{i}) == "spikets_trials_base_whisker")||...
                        (string(names{i}) == "spikets_trials_whisker")
                    
                    for Neuron = 1:length(Data.(names{i}){Exp})
                        Data_selected.(names{i}){Exp}{Neuron} =...
                            Data.(names{i}){Exp}{Neuron}(Ind);
                    end
                    
                elseif (string(names{i}) ~= "clusterID")&&...
                        (string(names{i}) ~= "ClusterCounter")&&...
                        (string(names{i}) ~= "Area")&&...
                        (string(names{i}) ~= "width")&&...
                        (string(names{i}) ~= "type")&&...
                        (string(names{i}) ~= "MDS")&&...
                        (string(names{i}) ~= "spikets_trials_base_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_base_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_base_FA_ITI")&&...
                        (string(names{i}) ~= "spikets_trials_FA_ITI")&&...
                        (string(names{i}) ~= "Whisking_AbsTimes_CR")&&...
                        (string(names{i}) ~= "JawOnsetsTms")&&...
                        (string(names{i}) ~= "LickOnsets_ITI_FA_LickTrace")
                    
                    Data_selected.(names{i}){Exp} = Data.(names{i}){Exp}(Ind);
                    
                end
            end
        end
        
    case 'FA_whisker'
        Data_selected = Data;
        for Exp = 1:length(Data.(names{1}))
            Ind = 1:length(Data.StimIndices{Exp});
            Ind = Ind(Data.FAIndices{Exp}==1);
            
            for i = 1:length(names)
                if (string(names{i}) == "spikets_trials_base_whisker")||...
                        (string(names{i}) == "spikets_trials_whisker")
                    
                    
                    for Neuron = 1:length(Data.(names{i}){Exp})
                        Data_selected.(names{i}){Exp}{Neuron} =...
                            Data.(names{i}){Exp}{Neuron}(Ind);
                    end
                    
                elseif (string(names{i}) ~= "clusterID")&&...
                        (string(names{i}) ~= "ClusterCounter")&&...
                        (string(names{i}) ~= "Area")&&...
                        (string(names{i}) ~= "width")&&...
                        (string(names{i}) ~= "type")&&...
                        (string(names{i}) ~= "MDS")&&...
                        (string(names{i}) ~= "spikets_trials_base_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_jaw")&&...
                        (string(names{i}) ~= "spikets_trials_base_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_base_FA_ITI")&&...
                        (string(names{i}) ~= "spikets_trials_FA_ITI")&&...
                        (string(names{i}) ~= "Whisking_AbsTimes_CR")&&...
                        (string(names{i}) ~= "LickOnsets_ITI_FA_LickTrace")
                    
                    Data_selected.(names{i}){Exp} = Data.(names{i}){Exp}(Ind);
                    
                end
            end
        end
        

    case 'Hit_jaw'
        Data_selected = Data;
        for Exp = 1:length(Data.(names{1}))
            Ind = 1:length(Data.StimIndices{Exp});
            Ind = Ind(Data.HitIndices{Exp}==1 & (~isnan(Data.JawOnsetsTms{Exp}))' & Data.StimAmps{Exp}~=1);
            
            for i = 1:length(names)
                if (string(names{i}) == "spikets_trials_base_jaw")||...
                        (string(names{i}) == "spikets_trials_jaw")
                    
                    
                    for Neuron = 1:length(Data.(names{i}){Exp})
                        Data_selected.(names{i}){Exp}{Neuron} =...
                            Data.(names{i}){Exp}{Neuron}(Ind);
                    end
                    
                elseif (string(names{i}) ~= "clusterID")&&...
                        (string(names{i}) ~= "ClusterCounter")&&...
                        (string(names{i}) ~= "Area")&&...
                        (string(names{i}) ~= "width")&&...
                        (string(names{i}) ~= "type")&&...
                        (string(names{i}) ~= "MDS")&&...
                        (string(names{i}) ~= "spikets_trials_base_whisker")&&...
                        (string(names{i}) ~= "spikets_trials_whisker")&&...
                        (string(names{i}) ~= "spikets_trials_base_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_whisking")&&...
                        (string(names{i}) ~= "spikets_trials_base_FA_ITI")&&...
                        (string(names{i}) ~= "spikets_trials_FA_ITI")&&...
                        (string(names{i}) ~= "Whisking_AbsTimes_CR")&&...
                        (string(names{i}) ~= "LickOnsets_ITI_FA_LickTrace")
                    
                    Data_selected.(names{i}){Exp} = Data.(names{i}){Exp}(Ind);
                    
                end
            end
        end
        
        
end

end