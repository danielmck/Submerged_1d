depth_ave_pe = load('Results/Ive_basal_pe_angles.txt');

depth_ave_h = load('Results/Ive_flow_height_angles.txt');
% deposit_time = sum(depth_ave_h>2e-3);
shear_lim_time = sum(depth_ave_h>0.999);
x = linspace(0.1,7,70);
y = linspace(1,500,50);
[Y,X] = meshgrid(y,x);
SetPaperSize(12,10);
hold on
[C,h] = contourf(Y,X,depth_ave_pe(1:50,:)',100,"HandleVisibility",'off');
% plot(10*shear_lim_time,linspace(0.1,7,70),'black','LineWidth',1.5,'DisplayName',"Start of Deposition");
% plot(10*deposit_time(2:end),linspace(0.1,7,70),'black','LineWidth',1.5,'DisplayName',"Start of Deposition");
set(h,'LineColor','none');
c = colorbar;
c.Label.String = "$p_e$";
ylabel("$\theta$")
xlabel("$t$")
ax = gca;
ax.XAxis.Exponent = 0;
% legend('Location', "southeast",'UserData', 10)
title("Initial Increase in $p_e$ in the Iverson Model");
PrintFig('Iverson_pe_base_early_angle')