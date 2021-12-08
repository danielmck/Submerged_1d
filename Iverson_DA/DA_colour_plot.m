depth_ave_h = load('DA_Results/DA_pe_base_angles_v3.txt');
x = linspace(0.1,7,70);
y = linspace(0,10000,1001);
[Y,X] = meshgrid(y,x);
SetPaperSize(12,10);
[C,h] = contourf(Y,X,depth_ave_h(1:1001,2:71)',100);
set(h,'LineColor','none')
c = colorbar;
c.Label.String = "Basal $p_e$";
ylabel("$\theta$")
xlabel("$t$")
ax = gca;
ax.XAxis.Exponent = 0;
title("Evolution of Basal $p_e$ in the Iverson Depth Averaged Model");
PrintFig('Ive_DA_pe_base_v3')