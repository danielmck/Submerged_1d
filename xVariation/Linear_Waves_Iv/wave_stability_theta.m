function wave_stability_theta

    d = 1e-3;
    phi_c=0.585; % Volume fraction
    
    g=9.81; % m/s^2

    rho_p = 2500; % kg/m^3
    
    rho_f = 1000; % kg/m^3
    eta_f = 0.0010016; 
    
%     rho_f = 1; % kg/m^3
%     eta_f = 1.18e-5; % Pa s
    
%     theta = 20; % deg
    alpha = 1e-4; % 1/Pa
    
      
    n_pts = 100;

    max_sig = zeros(n_pts);
    num_unstab = zeros(n_pts);
    A_mat = zeros(4);

    Fr_list = linspace(0.01,1.0,n_pts);
%         rho_list = (linspace(2.5,250,n_pts));
    theta = 15;
    theta_list = linspace(10,25,n_pts); %();()logspace(-5,log10(0.005),n_pts)
    for j = 1:n_pts
        Fr = Fr_list(j);
        for l = 1:n_pts  
            theta = theta_list(l);
            [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, 0,false,true);
%             d = d_dl*h0;
            stab_out = single_Fr_stab(Fr,crit_Iv,theta, rho_p, rho_f, d, eta_f, alpha);
            num_unstab(l,j) = stab_out(1);
            k_unstab(l,j) = stab_out(2);
        end
    end
      
%                 plot(0.25.*2.^((1:num_k)/2),sigma_mat(1:3,:))
        stability = (num_unstab>0);
        contourf(Fr_list,theta_list,stability,1)
               
        SetPaperSize(10,10)
        colormap(winter)
        xlabel('Froude Number')
%         ylabel('Particle size ($m$)')
%         set(gca, 'YScale', 'log')
        ylabel('Slope Angle ($^\circ$)')
%         c = colorbar('Ticks',[0,1,2]); 
%         c.Limits = ([0,2]);
%         c.Label.String = 'Number of Unstable Eigenmodes';
        title("$d="+num2str(d,'%.2e')+"$, $\alpha = "+num2str(alpha,'%.2e')+"$")
%         colorbar('Ticks',[0,0.5],'TickLabels',{'Stable','Unstable'});
        
        fig_name = 'StabCrit_Fr_theta_Norm';
        full_fig = strcat(fig_name,'.pdf');
        exp_graph(gcf,full_fig)
        movefile(full_fig, '../Figures/StabilityPlots');

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end