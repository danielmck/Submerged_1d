function plot_crit_h
    crit_h = load("Results/multi_comp_theta_Fr.txt");
    npts = size(crit_h,1);
    theta_start = 9;
    theta_stop = 20;
    theta_list = linspace(theta_start,theta_stop,npts);
    Fr_start = 0.6;
    Fr_stop = 5;
    Fr_list = linspace(Fr_start,Fr_stop,npts);
    SetPaperSize(10,10)
    contour(theta_list,Fr_list,crit_h);
    c = colorbar();
    c.Label.String = "$h_{crit}/h_0$";
    xlabel("Slope angle ($^\circ$)")
    ylabel("Froude number")
    title("$\tau_0 = 0Pa$, $\lambda/h_0 = 50$")
    exp_graph(gcf,"crit_h_theta_Fr.pdf")
end