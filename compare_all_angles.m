function compare_all_angles
    Rauter_ave = load("Rauter_Closure/EqnOfState_Results/depth_ave_angles.txt");
    DA_ave = load("Iverson_DA/DA_Results/DA_pe_base_angles_v3.txt");
    Iverson_ave = load("Iverson_Closure/Results/Ive_basal_pe_angles.txt");
    
    Iverson_h = load("Iverson_Closure/Results/Ive_flow_height_angles.txt");
    Ive_Is_eroding = Iverson_h>0.001;
    Ive_depo_start = sum(Iverson_h>0.999,1);
    Ive_depo_end = sum(Iverson_h>0.001,1);
    
    Rauter_h = load("Rauter_Closure/EqnOfState_Results/flow_height_angles.txt");
    Rau_Is_eroding = Rauter_h>0.999;
    Rau_Is_flowing = Rauter_h>0.001;
    Rau_depo_start = sum(Rauter_h>0.999,1);
    
    Ive_DA_diff = (DA_ave(:,2:end))./Iverson_ave*100.*Iverson_h;
    Rau_Ive_diff = (Iverson_ave)./max(Rauter_ave(:,2:end),1e-2)*100;
    
    Ive_DA_abs_diff = (DA_ave(:,2:end)-Iverson_ave);
    Rau_Ive_abs_diff = (Rauter_ave(:,2:end)-Iverson_ave);
    
    SetPaperSize(12,10);
    x = linspace(0.1,7,70);
    y = linspace(1,15000,1501);
    [Y,X] = meshgrid(y,x);
%     SetPaperSize(12,10);
    hold on
    [C,h] = contourf(Y,X,log2(Ive_DA_diff(1:1501,1:end)'/100),100,"HandleVisibility",'off');
    plot(10*Ive_depo_start,linspace(0.1,7,70),'black','LineWidth',1.5,'DisplayName',"Start of Deposition");
    plot(10*Ive_depo_end,linspace(0.1,7,70),'red','LineWidth',1.5,'DisplayName',"End of Deposition");
%     plot(10*Rau_Erosion_start(2:end),linspace(0.1,7,70),'red','LineWidth',1.5,'DisplayName',"Start of Rauter Deposition");
    legend('Location', "southeast",'UserData', 12)
    set(h,'LineColor','none')
    c = colorbar;
    c.Label.String = "Log Base 2 of Ratio of $\Bar{u}$";
%     caxis([0,200])
    ylabel("$\theta$")
    xlabel("$t$")
    title("Log of Ratio between $\Bar{u}$ in the Depth Averaged and Full Model")
    ax = gca;
    ax.XAxis.Exponent = 0;
    xlim([0,15000])
%     PrintFig('Ive_da_depth_ave_u_log_comp')
end