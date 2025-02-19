crit_vals = load("crit_h_tau0_theta_water.txt");
n_pts1 = size(crit_vals,1);
n_pts2 = size(crit_vals,2);

theta_start = 8.6;
theta_stop = 30;
theta_vals = linspace(theta_start,theta_stop,n_pts2);
C = viridis(n_pts1+1);
coeff_vals = [1e-5 1e-4 1e-3]; %[0,5,10,20,25]; %; % [1e-6 5e-6 1e-5 5e-5 1e-4];
% coeff_str = ["1\times 10^{-6}", "5\times 10^{-6}", "1\times 10^{-5}", "5\times 10^{-5}","1\times 10^{-4}"];
% coeff_str = ["1\times 10^{-5}", "1\times 10^{-4}", "1\times 10^{-3}"];
coeff_str = ["0","5","10","20","25"];

SetPaperSize(15,10)
% 
subplot(1,2,1)
hold on
for i=1:n_pts1
    plot(theta_vals,crit_vals(i,:),"DisplayName","$\tau_0="+coeff_str(i)+"$Pa",'color',C(i,:))
end
xlabel("$\theta$ ($^{\circ}$)")
ylabel("Critical $Fr$")
xlim([theta_start,theta_stop])
leg=legend('location','east');

ylim([0.0,0.6])
subplot(1,2,2)
hold on
for i=2:n_pts1
    plot(theta_vals,crit_vals(i,:),"DisplayName","$\tau_0="+coeff_str(i)+"$Pa",'color',C(i,:))
end
xlabel("$\theta$ ($^{\circ}$)")
ylabel("Critical $Fr$")
xlim([theta_start,theta_stop])
sgtitle("Fluid -- water, $\alpha =10^{-5}$ Pa$^{-1}$, $d =10^{-4}$m") %,, , $\tau_0 = 0$Pa

% exp_graph(gcf,"crit_Fr_tau0_water.pdf")