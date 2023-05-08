% Creates a plot of the value of the fiction coeff for different I
SetPaperSize(10,10);
points = logspace(-6,2,1000);
semilogx(points,mu_I_fn(points))
ylabel('$\mu(I)$')
xlabel('$I$')
% ax = gca;
% ax.YAxis.Exponent = 0;
% ax.XAxis.Exponent = 0;
% ytickformat('%.1e');
PrintFig('mu_profile')    