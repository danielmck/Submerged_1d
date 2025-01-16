SetPaperSize(10,10)
Iv_in = linspace(0,5e-4,100);
mu_in = mu_Iv_fn(Iv_in);
mu_in = max(0.34,mu_in);
plot(Iv_in,mu_in)
xlabel('$I_v$')
ylabel('$\mu(I_v)$')
exp_graph(gcf,'mu_Iv_plot.pdf')