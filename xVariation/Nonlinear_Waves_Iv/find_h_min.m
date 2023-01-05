function h_min = find_h_min(h_crit,Q1,u_w,theta,tau0_dl,crit_Iv)
    mu1_Iv = 0.32;
    
    g=9.81;

    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    
    rho_p = 2500;
    phi_c=0.585; % Volume fraction
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    h_stop_dl = (tau0_dl*rho_f/rho)/(tand(theta)-P*mu1_Iv);
    
    if h_stop_dl>Q1/u_w
        h_min = Q1/u_w;
    else
        h_min = (Q1/u_w+h_crit)/2;
        fb = 1;
        tol = 1e-6;
        while (abs(fb)>tol)
            fb = force_bal(h_min);
            fb_deriv = force_bal_deriv(h_min);
            h_min = h_min - fb/fb_deriv;
        end
    end

    function ud = get_u_deriv(h)
        ud = (-Q1 + h.*u_w)./h;
    end

    function u = get_u(h)
        u = Q1./h.^2;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = crit_Iv.*2.*u/h.^2;
        num_val = tand(theta)-mu_Iv_fn(Iv).*(rho-rho_f)/rho-tau0_dl*rho_f/rho/h;
    end

    function num_val = force_bal_deriv(h)
        u = get_u(h);
        ud = get_u_deriv(h);
        Iv = crit_Iv.*2.*u/h.^2;
        Iv_deriv = crit_Iv.*2.*ud./h.^2-2.*Iv/h;
        num_val = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
    end
end