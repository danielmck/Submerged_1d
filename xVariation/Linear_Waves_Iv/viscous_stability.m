function viscous_stability
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    theta = 12;
    rho_p = 2500;
    rho_f = 1000;
    phi_c = 0.585;
    g=9.81;
    eta_f = 0.0010016;
    nu = 1.13e-4;
    F=0.8;
    
    rho = rho_p*phi_c+(1-phi_c)*rho_f;
    crit_Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
    u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
    h0 = ((F*sqrt(g*cosd(theta)))./u_const)^(2/3);  
    u_eq = u_const.*h0^2;
    nu_dl = nu/(u_eq*h0);
    
    P = (rho-rho_f)/rho;
    
    dIvdu = crit_Iv;
    dIvdh = -2*crit_Iv;

    dmudu = dmudIv_fn(crit_Iv).*dIvdu;
    dmudh = dmudIv_fn(crit_Iv).*dIvdh;
    
    npt = 100;
    max_imag = zeros(1,npt);
    lamda_val = linspace(8,12,npt);
    for i = 1:npt
        lambda = lamda_val(i);
        k = 2*pi/lambda;
        c1 = -1;
        c2 = 2*k-1i/F^2*(P*dmudu)-1i*k^2*nu_dl;
        c3 = k^2*(1/F^2-1)+1i*k/F^2*(P*dmudu-P*dmudh)+1i*k^3*nu_dl;
        r = roots([c1, c2, c3]);
        max_imag(i) = max(imag(r));
    end
    plot(lamda_val,max_imag)
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end