function crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f, var_rho, phi_param)
    mu1_Iv = 0.32;
    if ~exist("var_rho","var")
        var_rho = false;
    end
    if ~exist("phi_param","var")
        phi_param = var_rho;
    end

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    pp = (rho-rho_f);

    crit_mu_c = rho/pp*tand(theta);
    crit_Iv = 2e-6;
    resid = 1;
    max_tol = 1e-7;
    del_Iv = 1e-8;
    if (crit_mu_c < mu1_Iv)
        warning("Slope is too shallow to sustain flow, try with a larger angle");
        crit_Iv = -1;
    else
        while (abs(resid)>max_tol)
            if var_rho
                phi = phi_c/(1+sqrt(crit_Iv));
                rho = phi_param*(phi*rho_p + (1-phi)*rho_f)+(1-phi_param)*(phi_c*rho_p + (1-phi_c)*rho_f);
                phi_deriv = -phi_c/2/(1+sqrt(crit_Iv))^2/sqrt(crit_Iv);
                rho_deriv = phi_param*(rho_p-rho_f)*phi_deriv;
                pp = (rho-rho_f);
                crit_mu = rho/pp*tand(theta);
                mu_deriv = rho_f/(rho-rho_f)^2*tand(theta)*rho_deriv;
            else
                crit_mu = crit_mu_c;
                mu_deriv = 0;
            end
            resid = mu_Iv_fn(crit_Iv)-crit_mu;
            if (crit_Iv>1e-6)
                crit_Iv = crit_Iv - resid./(dmudIv_fn(crit_Iv)+mu_deriv);
            else
                dmudIv = (mu_Iv_fn(crit_Iv+del_Iv) - mu_Iv_fn(crit_Iv))/del_Iv;
                crit_Iv = crit_Iv - resid./(dmudIv+mu_deriv);
            end
            if (crit_Iv<0)
                warning("Newton solver has failed to converge, try a larger initial value");
                crit_Iv = -1;
                break
            end
        end
    end
end