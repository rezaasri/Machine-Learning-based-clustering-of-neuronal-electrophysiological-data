function [CDF,PDF,N] = pdf_cdf_spikes(Spikes,dT_base,dT_after,...
                                                  dt_res_pdf,dt_res_cdf,j,...
                                                  pdf_plot_ind,...
                                                  cdf_plot_ind,...
                                                  spikes_plot_ind)
    
    %% CDF for the sequence
    N = length(Spikes);
    
    if j ==1  %%for baseline
        t_cdf_base = 0:dt_res_cdf:dT_base;
        CDF = zeros(size(t_cdf_base));
        for s = Spikes
            CDF = CDF + (t_cdf_base>=s)/N;
        end
    elseif j ==2  %% for after
        t_cdf_after = 0:dt_res_cdf:dT_after;
        CDF = zeros(size(t_cdf_after));
        for s = Spikes
            CDF = CDF + (t_cdf_after>=s)/N;
        end
    end
    
    

    %% PDF for the sequence
    
    if j ==1
        t_pdf_base = dt_res_pdf:dt_res_pdf:dT_base;
        PDF = zeros(size(t_pdf_base));
        for s = Spikes
            PDF = PDF + ((t_pdf_base>=s) - (t_pdf_base>=(s+ dt_res_pdf)));
        end
    elseif j ==2
        t_pdf_after = dt_res_pdf:dt_res_pdf:dT_after;
        PDF = zeros(size(t_pdf_after));
        for s = Spikes
            PDF = PDF + ((t_pdf_after>=s) - (t_pdf_after>=(s+ dt_res_pdf)));
        end
    end
    

    %% Plot
    if (pdf_plot_ind==1)||(spikes_plot_ind==1)||(cdf_plot_ind==1)
        figure
        if (pdf_plot_ind==1)
            bar(t_pdf-(dt_res_pdf/2), PDF/max(PDF), 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.5);
        end
        hold on
        if (cdf_plot_ind==1)
            plot(t_cdf, CDF, 'black');
        end
        hold on
        if (spikes_plot_ind==1)
            stem(Spikes, ones(size(Spikes)), 'Color','#D95319')
        end
        axis([0,dT,0,1])
        grid on
    end

end