function crit_I = newt_solve_crit_I(theta, rho_p, rho_f)
%     mu1_I=0.342; 
%     mu2_I=0.557;  
%     I_0 = 0.069;
    
    mu1_I=tand(20); 
    mu2_I=tand(33);  
    I_0 = 0.3;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;

    crit_mu = rho/(rho-rho_f)*tand(theta);
    crit_I = 2e-6;
    resid = 1;
    max_tol = 1e-7;
    del_I = 1e-8;
    if (crit_mu < mu1_I)
        "Slope is too shallow to sustain flow, try with a larger angle"
        crit_I = -1;
    else
        while (abs(resid)>max_tol)
            resid = mu_I_fn(crit_I)-crit_mu;
            if (crit_I>1e-6)
                crit_I = crit_I - resid./dmudI_fn(crit_I);
            else
                dmudIv = (mu_I_fn(crit_I+del_I) - mu_I_fn(crit_I))/del_I;
                crit_I = crit_I - resid./dmudIv;
            end
            if (crit_I<0)
                "Newton solver has failed to converge, try a larger initial value"
                crit_I = -1;
                break
            end
        end
    end
        
    function mu_val = mu_I_fn(I)
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end

    function dmudI = dmudI_fn(I)
        dmudI = (mu2_I-mu1_I)*I_0./(I_0+abs(I)).^2;
    end
end