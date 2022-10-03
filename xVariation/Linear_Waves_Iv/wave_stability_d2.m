function wave_stability_d2
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    
    g=9.81; % m/s^2
    
%     rho_f = 1; % kg/m^3
%     eta_f = 1.18e-5; % Pa s

    rho_f = 1000; % kg/m^3
    eta_f = 0.0010016;
    
    rho_p = 2500; % kg/m^3
    theta = 9; % deg
    rho = phi_c*rho_p + (1-phi_c)*rho_f;

    d = 1e-3;
    alpha = 1e-4; % 1/Pa
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);

    if (crit_Iv>0)
        crit_phi = phi_c./(1+sqrt(crit_Iv));

        n_pts = 500;

        max_sig = zeros(n_pts);
        num_unstab = zeros(n_pts);
        k_unstab = zeros(n_pts);
        A_mat = zeros(4);

        Fr_list = linspace(0.005,5,n_pts);
        d_list = logspace(-6,log10(5e-4),n_pts);

        for j = 1:n_pts
            Fr = Fr_list(j);
            for l = 1:n_pts  
                alpha = d_list(l);
                stab_out = single_Fr_stab(Fr,crit_Iv,theta, rho_p, rho_f, d, eta_f, alpha);
                num_unstab(l,j) = stab_out(1);
                k_unstab(l,j) = stab_out(2);
            end
        end
        
        stability = (max_sig>0);
        end_str = num2str(theta)+"deg_water_big_part_2mode.txt";
%         save("Results/nu_alpha_Fr_"+end_str, 'num_unstab','-ascii')
        save("Results/k_alpha_Fr_"+end_str, 'k_unstab','-ascii')
               
    end        
%%      
%         num_unstab = load("Results/nu_alpha_Fr_18deg_air_big_part.txt");
        f = figure;
        width = 10;
        height = 10;
        set(f, 'PaperUnits', 'centimeters');
        set(f, 'PaperSize', [width height]);
%         set(f,'Position',[15 15 width height])
        contourf(Fr_list,d_list,num_unstab,2)
        ax = gca;
        ax.YAxis.Exponent = 0;
        ytickformat('%5.1e');
%         SetPaperSize(10,10)
        colormap(winter)
        xlabel('Froude Number')
        ylabel('Particle Compressibility')
        c = colorbar;
        c.Ticks = [0 1 2];
        c.Limits = [0 2];
        c.Label.String = 'Number of Unstable Eigenmodes';
        title(strcat("$\theta = ",num2str(theta),"$, $ d =",num2str(d,'%5.1e'),"$, Fluid = Water"))
%         legend()
        set(gca,'YScale', 'log')
%         set(gca,'TickLabelInterpreter','latex')
        fig_name = 'UnstabMode_9deg_water.png';
%         full_fig = strcat(fig_name,'.pdf');
%         exportgraphics(f,full_fig,'Resolution',300)
%         PrintFig(fig_name)
%         exp_graph(f,fig_name)
%         movefile(fig_name, 'Figures/StabilityPlots');
        
%         [Fr1,alpha1] = meshgrid(Fr_list,d_list);
%         [x1,y1] = contour(Fr1,alpha1,num_unstab,[0.5 1.5],'*k');

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end