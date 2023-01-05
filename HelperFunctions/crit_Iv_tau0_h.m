function [crit_Fr, Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h, tau0, dl)
    mu1_Iv = 0.32;

    count = 0;
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    if tau0>0
        max_tol = 1e-7;
        resid = 1;   
        del_Iv = 1e-8;
        Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
        Iv_max = Iv;
        Iv_min = 1e-8;
        pp = (rho-rho_f)*g*h*cosd(theta);
        if dl
            tau0 = tau0*rho_f*g*cosd(theta)*h;
        end
        crit_mu = rho/(rho-rho_f)*tand(theta)-tau0/pp;
        if crit_mu>mu1_Iv
            while (abs(resid)>max_tol)
                resid = mu_Iv_fn(Iv)-crit_mu;
                if Iv > 1e-7
                    force_deriv = dmudIv_fn(Iv);
                else
                    force_deriv = ((mu_Iv_fn(Iv+del_Iv)-crit_mu)-resid)/del_Iv;
                end
                next_Iv = Iv-resid/force_deriv;
                if next_Iv<Iv_min
                    next_Iv = (Iv_min+Iv_max)/2;
                elseif next_Iv>Iv_max
                    next_Iv = (Iv_min+Iv_max)/2;
                end
                if resid<0
                    Iv_min = Iv;
                else
                    Iv_max = Iv;
                end
                Iv = next_Iv;
                count=count+1;
            end
            crit_Fr=get_Fr(Iv);
        else
            Iv=0;
            crit_Fr=0;
        end
    else
        Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
        crit_Fr=get_Fr(Iv);
    end
    
    function Fr=get_Fr(Iv)
        Fr = (1/3/eta_f*(Iv*(rho-rho_f)*sqrt(g*cosd(theta))))*h^(3/2);
    end
end