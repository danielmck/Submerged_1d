function [Fr_eq, crit_Iv] = crit_Iv_rauter_h(theta, rho_p, rho_f, eta_f, h0, a, phi_rlp, phi_rcp, tau0, dl, delta)
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
%             crit_h = (3*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
    end
    crit_rho = get_rho(phi);
    p_p = (crit_rho-rho_f)*g*cosd(theta)*(h0*(1-dl)+dl/(rho_f*g*cosd(theta)));
    iv_phi = (phi_c/phi-1)^2;
    if delta>0
        pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);   
    else
        pp_contact = 0;
    end
    pp_shear = p_p-pp_contact;
    u = pp_shear*iv_phi*h0/3/eta_f; 
    Fr_eq = u/sqrt(g*cosd(theta)*h0);
    crit_Iv = 3*u*eta_f/h0/p_p;

    function fb = force_bal_fn(phi)
        rho = get_rho(phi);
        p_p = (rho-rho_f)*g*cosd(theta)*(h0*(1-dl)+dl/(rho_f*g*cosd(theta)));
        iv_phi = (phi_c/phi-1)^2;
        if delta>0
            pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);   
        else
            pp_contact = 0;
        end
        pp_shear = max(p_p-pp_contact,0);
        u = pp_shear*iv_phi*h0/3/eta_f;    
        crit_mu = rho/(rho-rho_f)*tand(theta)-tau0/p_p;
        Iv = 3*u*eta_f/h0/p_p;
        fb = mu_Iv_fn(Iv)-crit_mu;
    end

    function rho = get_rho(phi)
        rho = delta*(phi*rho_p + (1-phi)*rho_f)+(1-delta)*(phi_c*rho_p + (1-phi_c)*rho_f);
    end

    function ud=get_u_deriv(phi)
        rho = get_rho(phi);
        Iv_phi = (phi_c/phi-1)^2;
        p_p = (rho-rho_f)*g*cosd(theta)*(h0*(1-dl)+dl/(rho_f*g*cosd(theta)));
        ppd = g*cosd(theta)*((rho_p-rho_f)*h0);
        if delta>0
            pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);
            pp_contact_d = a*(phi_rcp-phi_rlp)/(phi_rcp-phi)^2;
        else
            pp_contact = 0;
            pp_contact_d = 0;
        end
        pp_shear = p_p-pp_contact;
        pp_shear_d = ppd-pp_contact_d;
        Iv_phi_d = -2*(phi_c/phi-1)*phi_c/phi^2;
        ud = (pp_shear_d*Iv_phi+Iv_phi_d*pp_shear)*h0/3/eta_f;
    end

    function deriv = force_bal_deriv(phi)
        rho = get_rho(phi);
        p_p = (rho-rho_f)*g*cosd(theta)*(h0*(1-dl)+dl/(rho_f*g*cosd(theta)));
        iv_phi = (phi_c/phi-1)^2;
        if delta>0
            pp_contact = a*(phi-phi_rlp)/(phi_rcp-phi);   
        else
            pp_contact = 0;
        end
        pp_shear = max(p_p-pp_contact,0);
        u = pp_shear*iv_phi*h0/3/eta_f;  
        ud = get_u_deriv(phi);
        ppd = g*cosd(theta)*((rho_p-rho_f)*h0);
        Iv = 3*u*eta_f/h0/p_p;
        dIvdphi = 3*ud*eta_f/h0/p_p-Iv/p_p*ppd;
        deriv = dmudIv_fn(Iv)*dIvdphi-tau0/p_p^2*ppd+rho_f/(rho_p-rho_f)/phi^2*tand(theta);
    end
end