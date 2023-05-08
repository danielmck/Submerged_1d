function Fr_out = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha)
% Solves for the critical Froude number for linear instability for waves to
% form with single_Fr_stab and halving the potential span.
    g=9.8;
    phi_c = 0.585;

    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    if (crit_Iv > 0)
        
        Fr_min = 1e-4;
        Fr_max = 1.5;
        Fr_tol = 1e-5;

        while ((Fr_max-Fr_min)>Fr_tol)
            Fr_mid = (Fr_max+Fr_min)/2;
            stab_out = single_Fr_stab(Fr_mid,crit_Iv,theta, rho_p, rho_f, d, eta_f, alpha);
            num_unstab = stab_out(1);
            stable = (num_unstab==0);
            if stable
                Fr_min = Fr_mid;
            else
                Fr_max = Fr_mid;
            end
        end
        Fr_out = (Fr_max+Fr_min)/2;
        h_out = ((Fr_out*sqrt(g*cosd(theta)))./u_const)^(2/3);
    else
        Fr_out = -1;
        h_out = -1;
    end
end