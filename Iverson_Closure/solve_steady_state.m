function solve_steady_state
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 10^8;

    phi_c=0.585; % Volume fraction
    eta_f_dl = 3.9973e-5; % Pa s
    rho_r = 2.5;
    flux_list = linspace(0.1,5,50);
    theta_list = zeros(50,1);
    phi_list = zeros(50,1);
    for i=1:50
        psi = 0.1*i;
        u_bar = 1.5*psi;
        theta = newton_raph(8.7,u_bar,1e-4,0,1000);
        p_grad = (rho_r-1)*phi_c*cosd(theta);
        Iv = 2*u_bar.*eta_f_dl./p_grad;
        phi = phi_c./(1+sqrt(abs(Iv)));
        theta_list(i) = theta;
        phi_list(i) = phi;
    end
    out_mat = horzcat(flux_list',theta_list,phi_list);
    save('flux_conditions.txt', 'out_mat','-ascii')
    
    function theta_n = newton_raph(theta_o, u_bar, tol, it_counter, max_it)
        theta_del = 0.01;
        dudt_val = acc_rate(theta_o, u_bar);
        dudt_diff = (acc_rate(theta_o+theta_del, u_bar)-dudt_val)./theta_del;
        theta_n = theta_o - dudt_val./dudt_diff;
        dudt_val_n = acc_rate(theta_n, u_bar);
        if (it_counter < max_it)
            if (abs(dudt_val_n) > tol)
                theta_n = newton_raph(theta_n, u_bar, tol, it_counter+1, max_it);
            end
        end
    end
    
    function dudt = acc_rate(theta, u_bar)
        p_grad = -(rho_r-1)*phi_c*cosd(theta);
        Iv = -2*u_bar*eta_f_dl./p_grad;
        dudt = mu_Iv_fn(Iv)*p_grad - eta_f_dl.* 2*u_bar + sind(theta) * (rho_r*phi_c + (1-phi_c));
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end