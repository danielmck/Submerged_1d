function calc_flux_thresh
    flux_contours = zeros(20,50);
    flux_con = load('Results/flux_conditions.txt');
    max_pe = zeros(50,1);
    pe_time = zeros(50,1);
    for l=1:50
        flux=flux_con(l,1);
        theta_init=flux_con(l,2);
        init_phi=flux_con(l,3);
        unit = floor(flux);
        tenth = mod(l,10);
        if (tenth == 0)
            fname = "Ive_comp_4_deep_"+num2str(unit)+"_flux.txt";
        else
            fname = "Ive_comp_4_deep_"+num2str(unit)+"_"+num2str(tenth)+"_flux.txt";
        end
        cd Results
        data = load(fname);
        cd ..
        len_data = size(data,1);
        flux_vals = zeros(len_data,1);
%         for j = 1:len_data
%             flux_vals(j) = calc_flux(data(j,402:601));
%         end
%         for j = 1:19
%             time_index = sum(flux_vals>j*0.05*flux);
%             flux_contours(j+1,l) = data(time_index,1);
%         end
%         static_index = sum(flux_vals>1e-4);
%         flux_contours(1,l) = data(static_index,1);
        max_pe(l) = max(max(data(:,2:201)));
        [x,~]=find(data(:,2:201)==max_pe(l));
        pe_time(l) = data(x(end),1);
    end
    save("max_pe.txt", 'max_pe','-ascii')
    save("max_pe_time.txt", 'pe_time','-ascii')
end