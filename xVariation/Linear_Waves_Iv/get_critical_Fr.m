function h_temp = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha, tau0)
% Solves for the critical Froude number for linear instability for waves to
% form with single_Fr_stab and halving the potential span.

%     [h0, crit_Iv] = crit_Iv = 
%     u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
%     if (crit_Iv > 0)
        
    Fr_min = 1e-5;
    Fr_max = 1.5;
    Fr_tol = 1e-6;

    while ((Fr_max-Fr_min)>Fr_tol)
        Fr_mid = (Fr_max+Fr_min)/2;
        [h_temp, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_mid, tau0,false,true);
        stab_out = single_Fr_stab(Fr_mid,crit_Iv,theta, rho_p, rho_f, d, eta_f, alpha,tau0);
        % stab_out = single_stab_no_pe(Fr_mid,crit_Iv,theta, rho_p, rho_f, eta_f, tau0);
        num_unstab = stab_out(1);
        stable = (num_unstab==0);
        if stable
            Fr_min = Fr_mid;
        else
            Fr_max = Fr_mid;
        end
    end
    Fr_out = (Fr_max+Fr_min)/2;
    
        
%         h_out = ((Fr_out*sqrt(g*cosd(theta)))./u_const)^(2/3);
end