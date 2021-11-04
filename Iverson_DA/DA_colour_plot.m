depth_ave_h = load('DA_Results/DA_pe_base_angles_v2.txt');
x = linspace(0.1,7,70);
y = linspace(1,500,50);
[Y,X] = meshgrid(y,x);
SetPaperSize(12,10);
[C,h] = contourf(Y,X,depth_ave_h(1:50,2:71)',100);
set(h,'LineColor','none')
c = colorbar;
c.Label.String = "$p_e$";
ylabel("$\theta$")
xlabel("$t$")
ax = gca;
ax.XAxis.Exponent = 0;
title("Initial Increase in $p_e$ in the Iverson Depth Averaged Model");
PrintFig('Ive_DA_pe_base_early_v2')