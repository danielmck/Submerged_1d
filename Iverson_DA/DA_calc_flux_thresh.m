function DA_calc_flux_thresh
    flux_contours = zeros(20,50);
    flux_con = load('../Iverson_Closure/Results/flux_conditions.txt');
    for l=1:50
        flux=flux_con(l,1);
        theta_init=flux_con(l,2);
        init_phi=flux_con(l,3);
        unit = floor(flux);
        tenth = mod(l,10);
        if (tenth == 0)
            fname = "Ive_DA_4_deep_"+num2str(unit)+"_flux.txt";
        else
            fname = "Ive_DA_4_deep_"+num2str(unit)+"_"+num2str(tenth)+"_flux.txt";
        end
        cd Results
        data = load(fname);
        cd ..
        for j = 1:19
            time_index = sum(data(:,4)>j*0.05*flux);
            flux_contours(j+1,l) = data(time_index,1);
        end
        static_index = sum(data(:,4)>1e-4);
        flux_contours(1,l) = data(static_index,1);
    end
    save("DA_flux_contours.txt", 'flux_contours','-ascii')
end