function plot_2d_stab
    phase = 'water';
    var = 'd';
    crit_Fr = load("Results/crit_Fr_"+var+"_theta_"+phase+".txt");
%     
    if strcmp(var, 'alpha')
        var_list = [1e-6 5e-6 1e-5 5e-5 1e-4];
        var_name = "Compressibilies";
        slash = "\";
    else
        var_list = [2e-5 2e-4 2e-3];
        var_name = "Particle Size";
        slash = "";
    end
    if strcmp(phase, 'air')
        x_min = 17.8;
        x_max = 35;
    else
        x_min = 8.8;
        x_max = 30;
    end
    theta_list = (linspace(x_min,x_max,100));
%     Iv_list = zeros(100,1);
%     for i = 1:100
%         Iv_list(i) = newt_solve_crit_Iv(theta_list(i),2500,1);
%     end
%     xlim([x_min x_max])
    SetPaperSize(10,10)
    hold on
    for i=1:size(var_list,2)
        plot(theta_list,crit_Fr(i,:),'DisplayName',"$"+slash+var+"=$ "+num2str(var_list(i)))
    end
    ylim([0 0.6])
%     xlim([0 1e-3])
    ax = gca;
    ax.YAxis.Exponent = 0;
    YL = get(gca, 'YLim');
    plot([1 1]*x_min, [0 0.6],'--','Color',[0.2 0.2 0.2 0.5],'DisplayName','Min Slope Angle')
    lines1 = get(gca, 'Children');
    lines2 = vertcat(lines1(2:end),lines1(1));
    set(gca, 'Children', lines2 )
%     xtickformat('%5.1e');
    legend('Location', "best",'UserData', 9)
    xlabel("Slope angle $\theta$")
%     xlabel("Viscous Number $I_v$")
    ylabel("Froude Number $Fr$")
%     ylabel("Flow Height $h$")
    
    title("Critical $Fr$ for Different "+var_name+" and Slope Angles for "+strcat(upper(phase(1)),phase(2:end)))
    plot_name = strcat('Crit_Fr_',var,'_',phase);
    PrintFig(plot_name)
    movefile(plot_name+".pdf","Figures/StabilityLinePlots")
end