function crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;

    crit_mu = rho/(rho-rho_f)*tand(theta);
    crit_Iv = 2e-6;
    resid = 1;
    max_tol = 1e-7;
    del_Iv = 1e-8;
    if (crit_mu < mu1_Iv)
        "Slope is too shallow to sustain flow, try with a larger angle"
        crit_Iv = -1;
    else
        while (abs(resid)>max_tol)
            resid = mu_Iv_fn(crit_Iv)-crit_mu;
            if (crit_Iv>1e-6)
                crit_Iv = crit_Iv - resid./dmudIv_fn(crit_Iv);
            else
                dmudIv = (mu_Iv_fn(crit_Iv+del_Iv) - mu_Iv_fn(crit_Iv))/del_Iv;
                crit_Iv = crit_Iv - resid./dmudIv;
            end
            if (crit_Iv<0)
                "Newton solver has failed to converge, try a larger initial value"
                crit_Iv = -1;
                break
            end
        end
    end
        
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end