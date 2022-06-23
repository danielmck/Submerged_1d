depth_ave_h = load('EqnOfState_Results/depth_ave_angles.txt');
% deposit_time = sum(depth_ave_h>2e-3);
% shear_lim_time = sum(depth_ave_h>shear_lim_dl);
x = linspace(0,7,71);
y = linspace(1,20000,2001);
[Y,X] = meshgrid(y,x);
SetPaperSize(12,10);
hold on
xlim([500,20000]);
% colormap(jet(256));
c = colorbar;
c.Label.String = "$h$";
caxis([0,2.75])
[C,h] = contourf(Y,X,depth_ave_h(1:2001,:)',500,"HandleVisibility",'off');
% plot(10*shear_lim_time,linspace(0,7,71),'red','LineWidth',1.5,'DisplayName',"Reaches Shear Limit");
% plot(10*deposit_time,linspace(0,7,71),'black','LineWidth',1.5,'DisplayName',"Start of Deposition");
set(h,'LineColor','none')
% legend('Location', "southeast",'UserData', 14)

ylabel("$u_p$")
xlabel("$t$")
ax = gca;
ax.XAxis.Exponent = 0;

title("Decrease of depth averaged velocity in the Rauter Model");
PrintFig('Rauter_u_decay_angle')