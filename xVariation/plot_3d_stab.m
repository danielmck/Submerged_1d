% Creates a 3D plot of stability variation with slope angle, Froude number
% and relativity density
crit_Fr = load("crit_Fr_rho_theta.txt");
crit_Fr(crit_Fr<0) = nan;
rho_list = (linspace(2.5,25,100));
theta_list = (linspace(10,30,100));
SetPaperSize(10,10)
h = surf(theta_list,rho_list,crit_Fr);
xlabel("Slope angle $\theta$")
ylh = ylabel("Relative Density $\frac{\rho_p}{\rho_f}$");
% ylh.Position(1) = ylh.Position(1) - 5;
% ylh.Position(2) = ylh.Position(2) + 6;
zlabel("Froude Number $Fr$")
title("Surface of Critical $Fr$ for different $\theta$ and $\frac{\rho_p}{\rho_f}$")
rotate(h, [0 0 1],15)
set(gcf,'Renderer','Painter')
PrintFig('Crit_Fr_3D_plot')
% hgexport(gcf,'Crit_Fr_3D_plot.pdf');