function critical_Fr_val
    n_pts = 100;
    rho_p = 2500;
    eta_f = 0.0010016; % Pa s
    d_dl = 2e-3;
    alpha = 1e-4; % 1/Pa
    phi_c = 0.585;
    
    rho_list = (linspace(2.5,25,n_pts));
    theta_list = (linspace(10,30,n_pts));
    crit_Fr = zeros(100);
    for i = 1:n_pts
        rho_f = rho_p/rho_list(i);
        for j = 1:n_pts
            theta = theta_list(j);
            rho = phi_c*rho_p + (1-phi_c)*rho_f;
            crit_mu = rho/(rho-rho_f)*tand(theta);
            crit_Fr(i,j) = get_critical_Fr(theta, rho_p, rho_f, d_dl, eta_f, alpha);
        end
    end
    save("crit_Fr_rho_theta.txt", 'crit_Fr','-ascii')
end