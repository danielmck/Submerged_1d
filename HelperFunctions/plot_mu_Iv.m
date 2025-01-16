SetPaperSize(10,10)
Iv_in = linspace(0,1e-1,100);
mu_in = mu_Iv_fn(Iv_in);
mu_in = max(0.32,mu_in);
plot(Iv_in,mu_in, 'color', '#8da0cb','LineWidth',1.3);
xlabel('$I_v$')
ylabel('$\mu(I_v)$')
exp_graph(gcf,'mu_Iv_plot.pdf')
