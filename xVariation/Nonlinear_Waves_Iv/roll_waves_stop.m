function [xi_vals, h_vals] = roll_waves_stop
    theta = 12;
    Fr = 0.8;

    phi_c=0.585; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
%     h_crit = 0.8; % layer height (m)
    tau0 = 20;
    
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    nu = 5e-4;

    reg_param = 1*10^7;
    del_h = 1e-7;
    
    rho = phi_c*rho_p+(1-phi_c)*rho_f;

    pp_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
    
    
    [h_crit, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_const = crit_Iv.*pp_grad./eta_f./2;
    u_crit = u_const.*h_crit^2;
    pb_crit = rho_f*g*cosd(theta)*h_crit;
    
    z_scale = h_crit;
    v_scale = u_crit;
    p_scale = pb_crit;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    g_dl = g*t_scale/v_scale; 
    nu_dl = nu/v_scale/z_scale;
    
    pp_grad_dl = pp_grad/p_scale*z_scale;
    tau0_dl = tau0/p_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    
%     h_crit = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);
    
    
    h_crit_min = (0.5*sqrt(g*cosd(theta))/u_const)^(2/3);
    if (h_crit < h_crit_min)
        error('Flow not deep enough for travelling waves to occur');
        h_crit_min
    else
        crit_u = crit_Iv.*h_crit^2.*pp_grad./eta_f./2;
        Q1_max = sqrt(g_dl.*cosd(theta));
        u_w_max = 1 + sqrt(g_dl.*cosd(theta));
%         Q1 = h_stop^(3/2).*sqrt(g.*cosd(theta));
%         u_w = sqrt(h_stop.*g.*cosd(theta));
%         u_w= 2.3536;
%         Q1 = 1.2759;
        u_w = u_w_max-1e-4;
        Q1=u_w-1;
%         h_min = 0.8*h_stop;
%         h_max = get_h_plus(h_min);
        h_stop = tau0_dl*rho_f_dl/rho_dl/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);
        h_plus = Q1/u_w;

        opts = odeset('Events',@WaveEndFn);
        [xi_vals,h_vals]=ode15s(@h_deriv_fn,[0,100],[h_plus,0.1],opts);
        u_vals = (-Q1./h_vals+u_w);
        plot(xi_vals,h_vals)
%         xi_vals = xi_vals(h_vals<h2)/h_crit;
%         h_vals = h_vals(h_vals<h2)/h_crit;
%         xlabel('$\xi(m)$')
%         ylabel('h(m)')
%         title('Roll waves on a 15 degree slope')
%         PrintFig('Ive_15deg_roll_wave')
    end

    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(g_dl.*h_minus.*cosd(theta))))/2;
    end

    function h_minus = get_h_minus(h_plus)
        h_minus = (-h_plus + sqrt(h_plus.^2+8*Q1^2./(g_dl.*h_plus.*cosd(theta))))/2;
    end

    function [position,isterminal,direction] = WaveEndFn(xi,y)
        position = y(1) - h_plus;
        isterminal = 0;
        direction = 0;
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dydxi = h_deriv_fn(xi,y)
        h = y(1);
        u = u_w - Q1/h;
        n = y(2);

        dhdxi = n;
        n_coeff = 1-Q1^2.*Fr^2/h^3;
%             Fr = Fr_eq*abs(u)/h;
        Iv = crit_Iv*abs(u)/h^2;
        force_bal = tand(theta)-sign(u).*P*mu_Iv_fn(Iv)-tau0_dl*rho_f/rho/h;
        n_eq = (force_bal)./n_coeff;
        dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr^2/nu_dl/Q1*n_coeff*(n-n_eq);
        dydxi=[dhdxi,dndxi]';
    end

    function dhdxi = h_deriv_fn_back(xi,h)
        if (abs(h-h_crit)>1e-6)
            dhdxi = -g_dl.*h.^3.*force_bal(h)./(g_dl.*h.^3.*cosd(theta)-Q1.^2);
        else
            dhdxi = -(g_dl.*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(g_dl.*3.*h.^2.*cosd(theta));
        end
    end

    function u = get_u(h)
        u = (-Q1 + h.*u_w)./h;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = 2.*abs(u).*eta_f_dl./(pp_grad_dl.*h.^2);
        num_val = tand(theta)-sign(u)*mu_Iv_fn(Iv).*(rho-rho_f)/rho-sign(u)*tau0_dl*rho_f_dl/rho_dl/h;
    end
end