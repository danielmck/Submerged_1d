function Fr = get_Fr_from_flux(flux,theta,tau0)
    g=9.81;

    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    
    rho_p = 2500;
    phi_c=0.585; % Volume fraction
    rho = rho_p*phi_c+rho_f*(1-phi_c);

    pp_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
    
    Fr_min = 0;
    Fr_max = 5;
    Fr = 1;
    flux_diff = 1;
    while (abs(flux_diff) > 1e-6 && Fr_max-Fr_min>1e-7)
        [h,~] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0 ,0);
        flux_diff = flux-Fr*sqrt(g*cosd(theta))*h^(3/2);
        if (flux_diff<0)
            Fr_max = Fr;
        else
            Fr_min = Fr;
        end
        Fr = (Fr_max+Fr_min)/2;
    end
end