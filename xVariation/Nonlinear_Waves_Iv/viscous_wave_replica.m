function viscous_wave_replica
%  A replica of the visous wave from the viroulet paper in order to get the
%  process working.
    mu1 = tand(20.9);
    mu2 = tand(32.76);
    beta = 0.136;
    L = 8.25e-4;
    

    g=9.81; % m/s^2
    theta = 29;
    gamma = (mu2-tand(theta))/(tand(theta)-mu1);
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
%    crit_Fr = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
    Fr_eq = 0.8; 
%     Fr_eq = h0*beta/gamma/L;
    h0 = Fr_eq*L*gamma/beta;
    
    u_eq = Fr_eq*sqrt(g*h0*cosd(theta));
    nu = 1.13e-3;
    
    z_scale = h0;
    v_scale = u_eq;
    t_scale = z_scale/v_scale;

    g_dl = g*t_scale/v_scale; 
    nu_dl = nu/z_scale/v_scale;
    L_dl = L/z_scale;
    R = u_eq*sqrt(h0)/nu;
    
%     u_w = 1+sqrt(g_dl*cosd(theta));
%     Q1 = u_w-1;
    
%     Q1 = sqrt(h0^3*g*cosd(theta));
    u_w_dl = 1+1/Fr_eq-1e-6;
    Q1_dl = u_w_dl-1;

%     u_w = u_eq*u_w_dl;
%     Q1 = h0*(u_w-u_eq);
    
    
    
%     init_h = 1.001*h0;
%     h_max = (-init_h + sqrt(init_h.^2+8*Q1^2./(g.*init_h.*cosd(theta))))/2;
%     init_u = u_w+(1-u_w)/init_h;
%     init_n = 0.0;
%     init_w1 = Q1;
%     init_w3 = (1-u_w)/sqrt(init_h)*init_n;
% %     init_vals = [init_w1, init_h, init_w3];
%     init_vals = [init_w1, init_h, init_n];
    
%     h_list = linspace(0.5,1.5,101)*h0;
%     u_list = u_w - Q1./h_list;
%     Fr_list = u_list./sqrt(g*h_list*cosd(theta));
%     mu_list = mu1 + (mu2-mu1)./(beta.*h_list./(L.*Fr_list)+1);
%     n_eq_list = (tand(theta)-mu_list)./(1-Q1^2./(g*h_list.^3*cosd(theta)));
%     plot(h_list,mu_list)
    
    opts = odeset('Events',@WaveEndFn);
%     [xi_vals,out_vals]=ode15s(@full_system_orig,[0, 5e-2],init_vals);
    [xi_vals,out_vals]=ode15s(@full_system_orig_dl,[0, 500],[1.001,0]);
    [min_h, h_ind] = min(out_vals(:,1));
    n_min = out_vals(h_ind,2);
%     [xi_wave,out_wave] = ode15s(@full_system_orig_dl,[0, 50],out_vals(h_ind,:),opts);
    plot(xi_vals(:,1),out_vals(:,1))
    max_h = (-min_h + sqrt(min_h.^2+8*Q1_dl^2.*Fr_eq.^2./min_h))/2;
%     hold on
%     opts_nv = odeset('Events',@NV_WaveEndFn);
%     [xi_nv,out_nv] = ode15s(@non_viscous_system,[0, 50],min_h,opts_nv);
%     plot(xi_nv(:,1),out_nv(:,1))
    
    function dydx = full_system_orig(x,y)
        h = y(2);
        u = u_w - Q1/h;
        n = y(3);

        dy1dx = 0;
        dy2dx = n;
        n_coeff = 1-Q1^2/(g*h^3*cosd(theta));
        Fr = u/sqrt(g*h*cosd(theta));
        mu_b = mu1 + (mu2-mu1)/(beta*h/(L*Fr)+1);
        n_eq = (tand(theta)-mu_b)./n_coeff;
        dy3dx = 1/(2*h)*n^2 + g*h^(3/2)*cosd(theta)/nu_dl/Q1*n_coeff*(n-n_eq);
%         dy2dx = (g*h^3*cosd(theta))*(tand(theta)-mu_b)/(n_coeff);
        dydx = [dy1dx, dy2dx, dy3dx]';
    end

    function dhdx = non_viscous_system(x,y)
        h = y(1);
        u = u_w_dl - Q1_dl/h;
        n_coeff = 1-Q1_dl^2.*Fr_eq^2/h^3;
        Fr = Fr_eq*u/h;
%         mu_b1 = mu1 + (mu2-mu1)/(beta*h/(L_dl*Fr)+1);
        mu_b = mu1 + (mu2-mu1)*(1-u_w+u_w*h)/(1-u_w+u_w*h+h^(5/2)*gamma);
        if abs(n_coeff)>1e-6
            dhdx = (tand(theta)-mu_b)./n_coeff;
        else
            dhdx = (mu2-tand(theta))/2/Fr_eq/(gamma+1)*(Fr_eq-2/3);
        end
    end

    function dydx = full_system_orig_dl(x,y)
        h = y(1);
        u = u_w_dl - Q1_dl/h;
        n = y(2);

        dy2dx = n;
        n_coeff = 1-Q1_dl^2.*Fr_eq^2/h^3;
        Fr = Fr_eq*u/h;
        mu_b = mu1 + (mu2-mu1)/(beta*h/(L_dl*Fr)+1);
%         mu_b = mu1 + (mu2-mu1)*(1-u_w+u_w*h)/(1-u_w+u_w*h+h^(5/2)*gamma);
        n_eq = (tand(theta)-mu_b)./n_coeff;
        dy3dx = 1/(2*h)*n^2 + h^(3/2)/Fr_eq^2*R/Q1_dl*n_coeff*(n-n_eq);
%         dy2dx = (g*h^3*cosd(theta))*(tand(theta)-mu_b)/(n_coeff);
        dydx = [dy2dx, dy3dx]';
    end

    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = (t<0.1)+y(1)-min_h;
        isterminal = 1;
        direction = 0;
    end

    function [position,isterminal,direction] = NV_WaveEndFn(t,y)
        position = y(1)-max_h;
        isterminal = 1;
        direction = 0;
    end
end