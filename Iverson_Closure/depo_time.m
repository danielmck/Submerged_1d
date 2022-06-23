all_depo = zeros(23,1);
for l=1:23
    theta_init = 8.6+0.1*l;
%         theta = theta_init;
    unit = floor(theta_init);
    tenth = mod(l+6,10);
    if (tenth == 0)
        fname = "Ive_comp_4_deep_"+num2str(unit)+"_start.txt";
    else
        fname = "Ive_comp_4_deep_"+num2str(unit)+"_"+num2str(tenth)+"_start.txt";
    end
    cd Results
    data = load(fname);
    times = data(:,1);
    uf_top = data(:,600);
    flow_times = sum(uf_top > 1e-4);
    all_depo(l) = times(flow_times+1);
    cd ..
end
all_depo
%% 
SetPaperSize(10,10)
semilogy(linspace(8.7,11,23),all_depo)
xlabel("Slope Angle of ICs")
ylabel("Time for top of flow to come to rest")
title("Time of Final Deposition for different Slope Angles")
PrintFig('Ive_deposit_time_angle')
