depth_ave_h = load('EqnOfState_Results/flow_height_angles.txt');
% deposit_time = sum(depth_ave_h>2e-3);
% shear_lim_time = sum(depth_ave_h>shear_lim_dl);
x = linspace(0,7,71);
y = linspace(1,50000,5000);
[Y,X] = meshgrid(y,x);
SetPaperSize(12,10);
hold on
% xlim([0,15000]);
% colormap(jet(256));
c = colorbar;
c.Label.String = "$h$";
% caxis([0,10])
[C,h] = contourf(Y,X,depth_ave_h(1:5000,:)',50,"HandleVisibility",'off');
% plot(10*shear_lim_time,linspace(0,7,71),'red','LineWidth',1.5,'DisplayName',"Reaches Shear Limit");
% plot(10*deposit_time,linspace(0,7,71),'black','LineWidth',1.5,'DisplayName',"Start of Deposition");
set(h,'LineColor','none')
% legend('Location', "southeast",'UserData', 14)

ylabel("$\theta$")
xlabel("$t$")
ax = gca;
ax.XAxis.Exponent = 0;

title("Decrease of Flow Height in Rauter model");
PrintFig('Rauter_h_decay_angle')