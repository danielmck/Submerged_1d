 function [crit_h, Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_eq, tau0)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;
    count = 0;
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    max_tol = 1e-7;
    resid = 1;   
    del_Iv = 1e-8;
    Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
    while (abs(resid)>max_tol)
        resid = force_bal_fn(Iv);
        if Iv > 1e-7
            force_deriv = force_bal_deriv(Iv);
        else
            force_deriv = (force_bal_fn(Iv+del_Iv)-resid)/del_Iv;
        end
        Iv = max(Iv-resid/force_deriv,1e-8);
        crit_h = get_h(Iv);
        count=count+1;
    end

    
    function fb = force_bal_fn(Iv)
        h=get_h(Iv);
        pp = (rho-rho_f)*g*h*cosd(theta);
        crit_mu = rho/(rho-rho_f)*tand(theta)-tau0/pp;
        fb = mu_Iv_fn(Iv)-crit_mu;
    end
    
    function h=get_h(Iv)
        h = (2*Fr_eq*eta_f/(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))^(2/3);
    end

    function deriv = force_bal_deriv(Iv)
        deriv = dmudIv_fn(Iv) + tau0*2/3/(2*Fr_eq*eta_f/((rho-rho_f)*sqrt(g*cosd(theta))))^(2/3)/Iv^(1/3);
    end
        
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end