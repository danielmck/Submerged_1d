function roll_waves
    h_crit = 0.8; % layer height (m)
    d=1.43e-5; % grain diameter (m)

    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^3
    theta = 15; % deg

    Fr_min = 0.5;

    alpha = 0.0001; % 1/Pa


    % v_scale = sqrt(g.*h);
    % p_scale = rho_f.*g.*h;
    % t_scale = sqrt(h./g);
    % z_scale = h;
    rho = phi_c*rho_p+(1-phi_c)*rho_f;

    pb_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);

    crit_mu = rho./(rho-rho_f).*tand(theta);
    crit_Iv = 1e-4;
    del_Iv = 1e-6;
    resid = 1;
    max_tol = 1e-7;
    while (abs(resid)>max_tol)
        resid = mu_Iv_fn(crit_Iv)-crit_mu;
        diff_mu = (mu_Iv_fn(crit_Iv+del_Iv)-mu_Iv_fn(crit_Iv))./del_Iv;
        crit_Iv = crit_Iv - resid./diff_mu;
    end
    h_crit_min = (Fr_min*sqrt(g*cosd(theta))/(crit_Iv*pb_grad))^(2/3);
    if (h_crit < h_crit_min)
        'Flow not deep enough for travelling waves to occur'
        h_crit_min
    else
        crit_u = crit_Iv.*h_crit^2.*pb_grad./eta_f./2;
        Q1 = h_crit^(3/2).*sqrt(g.*cosd(theta));
        u_w = crit_u + sqrt(h_crit.*g.*cosd(theta));

        -4*u_w^3+27*Q1^2*crit_Iv*pb_grad./eta_f./2
        roots([crit_Iv*pb_grad./eta_f./2 0 -u_w Q1])
        h_min = h_crit_min;
        del_h = 1e-6;
        h_min_resid = 1;
        h_min_max_tol = 1e-7;
        while (abs(h_min_resid)>h_min_max_tol)
            h_min_resid = force_bal(h_min);
            diff_bal = (force_bal(h_min+del_h)-force_bal(h_min))./del_h;
            h_min = h_min - h_min_resid./diff_bal;
        end
        h_max = (-h_min + sqrt(h_min^2+8*Q1^2/(g*h_min*cosd(theta))))/2;
        hold on
        n_line = 10;
        h_low = 0.01;
        h_high = 0.1;
        p = h_crit*crit_u;
        fluxes = zeros(n_line,1);
%         SetPaperSize(12,10)
        for i = 1:n_line
            h1 = h_low+(i-1)/n_line*(h_high-h_low);
            h2 = get_h_plus(h1);
            [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,500,10001),h1);
            u_vals = get_u(h_vals);
            n_point = sum(h_vals<h2);
            plot(xi_vals(h_vals<h2),h_vals(h_vals<h2),'k')
            flux = (h_vals(1).*u_vals(1)+h_vals(n_point).*u_vals(n_point) + 2*sum(h_vals(2:n_point-1).*u_vals(2:n_point-1)))./(2*n_point);
            fluxes(i) = flux;
        end
        xlabel('$\xi(m)$')
        ylabel('h(m)')
        title('Roll waves on a 15 degree slope')
%         PrintFig('Ive_15deg_roll_wave')
    end

    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(g.*h_minus.*cosd(theta))))/2;
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dhdxi = h_deriv_fn(xi,h)
        if (abs(h-h_crit)>1e-6)
            dhdxi = g.*h.^3.*force_bal(h)./(g.*h.^3.*cosd(theta)-Q1.^2);
        else
            dhdxi = (g.*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(g.*3.*h.^2.*cosd(theta));
        end
    end

    function u = get_u(h)
        u = (-Q1 + h.*u_w)./h;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = 2.*u.*eta_f./(pb_grad.*h.^2);
        num_val = tand(theta)-mu_Iv_fn(Iv).*(rho-rho_f)/rho;
    end
end