function [xi_vals, h_vals, u_vals, n_vals] = roll_waves_lambda(lambda,theta,rho_f_dl,rho_p_dl,eta_f_dl,Fr_min)
% Creates roll waves of the two equation model that are close to a 
% prescribed value of lambda. Use the ode method of solving the roll waves.
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g_dl = 1/rho_f_dl/cosd(theta);

%     Fr_min = get_critical_Fr(theta, rho_p, rho_f, d_dl, eta_f, alpha);
    rho_dl = phi_c*rho_p_dl+(1-phi_c)*rho_f_dl;
    pp_grad_dl = (rho_p_dl-rho_f_dl)*g_dl*phi_c*cosd(theta);
    crit_Iv = newt_solve_crit_Iv(theta, rho_p_dl, rho_f_dl);
    crit_u = crit_Iv.*pp_grad_dl./eta_f_dl./2;
%     
%     v_scale = crit_u;
%     p_scale = rho_f.*g.*cosd(theta).*h_crit;
%     t_scale = sqrt(h./g);
%     z_scale = h_crit;
%     rho_f_dl = rho_f*v_scale^2/p_scale;
%     rho_p_dl = rho_p*v_scale^2/p_scale;
% %     rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
%     eta_f_dl = eta_f/(p_scale*t_scale);
    
    h_crit_min = (Fr_min*sqrt(g_dl*cosd(theta))/crit_u)^(2/3);
    if (1 < h_crit_min)
        error("Flow not deep enough for travelling waves to occur, must be above "+num2str(h_crit_min))  
    else
        
        Q1 = sqrt(g_dl.*cosd(theta));
        u_w = 1 + sqrt(g_dl.*cosd(theta));

%         determ = -4*u_w^3+27*Q1^2*crit_Iv*pb_grad./eta_f./2;
        eq_heights = roots([crit_Iv.*pp_grad_dl./eta_f_dl./2 0 -u_w Q1]);
        h_min = max(eq_heights(eq_heights<(1-1e-6) & eq_heights>0));
%         h_max = get_h_plus(h_min);
%         hold on
        n_line = 100;
        if (size(h_min,1) == 0)
            h_low = 0.1;
        else
            h_low = max(0.1,h_min+1e-3);
        end
        h_high = 0.999;
        h1_list = linspace(h_low,h_high,n_line);
%         SetPaperSize(12,10)
        wave_lambda = lambda+1;
        i = 1;
        while (wave_lambda>lambda) && (i<=n_line)
            h1 = h1_list(i);
            h2 = get_h_plus(h1);
            [xi_vals,h_vals,u_vals,n_vals] = single_roll_wave(h1,theta,rho_dl,rho_f_dl,g_dl,eta_f_dl,pp_grad_dl,crit_Iv);
            if h_vals(end)>0
                wave_lambda = xi_vals(end);
            end
            i=i+1;
        end
%         xlabel('$\xi(m)$')
%         ylabel('h(m)')
%         title('Roll waves on a 15 degree slope')
%         PrintFig('Ive_15deg_roll_wave')
    end

    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(g_dl.*h_minus.*cosd(theta))))/2;
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dhdxi = h_deriv_fn(xi,h)
        if (abs(h-1)>1e-6)
            dhdxi = g_dl.*h.^3.*force_bal(h)./(g_dl.*h.^3.*cosd(theta)-Q1.^2);
        else
            dhdxi = (g_dl.*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(g_dl.*3.*h.^2.*cosd(theta));
        end
    end

    function u = get_u(h)
        u = (-Q1 + h.*u_w)./h;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = 3.*u.*eta_f_dl./(pp_grad_dl.*h.^2);
        num_val = tand(theta)-mu_Iv_fn(Iv).*(rho_dl-rho_f_dl)/rho_dl;
    end
end