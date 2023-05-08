function plot_2d_stab_I
% Plots stability criteria for the linear stability analysis for different
% diameters and compressibilities.

    crit_Fr = load("crit_Fr_alpha_theta_water_I.txt");
    d_list = [1e-6 5e-6 1e-5 5e-5 1e-4];
%     d_list = [1e-5 5e-5 1e-4 5e-4 1e-3];
%     theta_list = (linspace(20,32,100));
    theta_list = (linspace(10,14,100));
    SetPaperSize(10,10)
    hold on
    for i=1:size(d_list,2)
        plot(theta_list,crit_Fr(i,:),'DisplayName',"$\alpha^{\ast}=$ "+num2str(d_list(i)))
    end
    ylim([0 0.7])
    legend('Location', "best",'UserData', 6)
    xlabel("Slope angle $\theta$")
    ylabel("Froude Number $Fr$")
    title("Critical $Fr$ for Different Compressibilities and Slope Angles for Water with $\mu(I)$")
    PrintFig('Crit_Fr_alpha_water_I')
end