 function [crit_h, crit_Iv] = crit_Iv_rauter(theta, rho_p, rho_f, eta_f, Fr_eq, a, phi_rlp, phi_rcp, tau0, dl, delta)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    
    if ~exist("dl","var")
        dl = false;
    end
    if ~exist("delta","var")
        delta = 1;
    end
    
    reg_param = 1*10^7;
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    max_tol = 1e-7;
    resid = 1;   
    del_phi = 1e-8;
    max_phi = phi_c-1e-8;
    min_phi = 1e-8;
    phi = phi_c/(1+sqrt(newt_solve_crit_Iv(theta, rho_p, rho_f, true)));
    while (abs(resid)>max_tol)
        resid = force_bal_fn(phi);
        if phi_c-phi >  1e-7
            force_deriv = force_bal_deriv(phi);
        else
            force_deriv = (force_bal_fn(phi+del_phi)-resid)/del_phi;
        end
        if resid > 0
            min_phi = phi;
        else
            max_phi = phi;
        end
        if ((phi-resid/force_deriv<min_phi) || (phi-resid/force_deriv>max_phi))
            phi = (max_phi+min_phi)/2;
        else
            phi = phi-resid/force_deriv;
        end
        crit_h = get_h(phi);
%             crit_h = (3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
    end
    crit_u = Fr_eq*sqrt(g*cosd(theta)*crit_h);
    crit_rho = get_rho(phi);
    crit_Iv = 3*eta_f*crit_u/crit_h/((crit_rho-rho_f)*g*cosd(theta)*crit_h);

    function fb = force_bal_fn(phi)
        rho = get_rho(phi);
        h = get_h(phi);
        p_p = (rho-rho_f)*g*cosd(theta)*(h*(1-dl)+dl/(rho_f*g*cosd(theta)));
%         h=(3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
        u = Fr_eq*sqrt(g*cosd(theta)*h);
        crit_mu = rho/(rho-rho_f)*tand(theta)-tau0/p_p;
        Iv = 3*u*eta_f/h/p_p;
        fb = mu_Iv_fn(Iv)-crit_mu;
    end

    function rho = get_rho(phi)
        rho = delta*(phi*rho_p + (1-phi)*rho_f)+(1-delta)*(phi_c*rho_p + (1-phi_c)*rho_f);
    end
    
    function h=get_h(phi)
        rho = get_rho(phi);
        Iv_phi = (phi_c/phi-1)^2;
        mult_term = Iv_phi/3/eta_f/sqrt(g*cosd(theta));
        pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);
        rt_h = roots([mult_term*(rho-rho_f)*g*cosd(theta), 0, -mult_term*delta*pp_contact, -Fr_eq]);
        rt_h = rt_h(imag(rt_h)==0);
        h = rt_h.^2;
    end

    function hd=get_h_deriv(phi)
        rho = get_rho(phi);
        h = get_h(phi);
        pp_eq = (rho-rho_f)*g*cosd(theta)*h;
        pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);
        Iv_phi = (phi_c/phi-1)^2;
        phi_part = -2*(phi_c/phi-1)*phi_c/phi^2*sqrt(h)*(pp_eq-delta*pp_contact)+Iv_phi*sqrt(h)*((rho_p-rho_f)*g*cosd(theta)*h-delta*a*(phi_rcp-phi_rlp)/(phi_rcp-phi)^2);
        h_part = Iv_phi*(1/2/sqrt(h)*delta*pp_contact-3/2*sqrt(h)*(rho-rho_f)*g*cosd(theta));
        hd = phi_part/h_part;
    end

    function ud=get_u_deriv(phi)
        h = get_h(phi);
        ud = 1/2*Fr_eq*sqrt(g*cosd(theta)/h)*get_h_deriv(phi);
    end

    function deriv = force_bal_deriv(phi)
        h=get_h(phi);
        rho = get_rho(phi);
        p_p = (rho-rho_f)*g*cosd(theta)*(h*(1-dl)+dl/(rho_f*g*cosd(theta)));
        u = Fr_eq*sqrt(g*cosd(theta)*h);
        hd=get_h_deriv(phi);
        ud = get_u_deriv(phi);
        ppd = g*cosd(theta)*((rho_p-rho_f)*h+(rho-rho_f)*hd);
        Iv = 3*u*eta_f/h/p_p;
        dIvdphi = 3*ud*eta_f/h/p_p-Iv/h*hd-Iv/p_p*ppd;
        deriv = dmudIv_fn(Iv)*dIvdphi-tau0/p_p^2*ppd+rho_f/(rho_p-rho_f)/phi^2*tand(theta);
    end
end