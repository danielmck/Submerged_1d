function [xi_vals,h_vals,u_vals,n_vals]=single_roll_wave(h1,theta,rho_dl,rho_f_dl,g_dl,eta_f_dl,pp_grad_dl,crit_Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    
%     if ~exist("h1","var")
%         h1 = 0.0023;
%         theta=12;
%         rho_f_dl = 1000;
%         rho_dl = rho_f_dl*(1-0.585)+2500*0.585;
%         
%     end

    reg_param = 1*10^7;
    
    del_h = 1e-8;

    phi_c=0.585; % Volume fraction

    Q1 = sqrt(g_dl.*cosd(theta));
    u_w = 1 + sqrt(g_dl.*cosd(theta));
    eq_heights = roots([crit_Iv.*pp_grad_dl./eta_f_dl./2 0 -u_w Q1]);
    h_min = max(eq_heights(eq_heights<(1-1e-6) & eq_heights>0));
    h2 = get_h_plus(h1);
    if (h1>h_min)
        opts = odeset('Events',@WaveEndFn);
        [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,100,201),h1,opts);
        u_vals = get_u(h_vals);
        n_vals = zeros(size(xi_vals));
        for i = 1:size(xi_vals,1)
            n_vals(i) = h_deriv_fn(xi_vals(i),h_vals(i));
        end
    else
        warning("Initial height is below minimum wave height")
        xi_vals = -1;
        h_vals = -1;
        u_vals = -1;
        n_vals = -1;
    end

    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = y - h2;
        isterminal = 1;
        direction = 0;
    end
    
    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(g_dl.*h_minus.*cosd(theta))))/2;
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dhdxi = h_deriv_fn(xi,h)
        if (abs(h-1)>1e-6)
            dhdxi = g_dl.*cosd(theta).*h.^3.*force_bal(h)./(g_dl.*h.^3.*cosd(theta)-Q1.^2);
        else
            dhdxi = (g_dl.*cosd(theta).*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(g_dl.*3.*h.^2.*cosd(theta));
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