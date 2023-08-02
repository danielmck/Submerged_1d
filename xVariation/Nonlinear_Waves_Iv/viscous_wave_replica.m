function viscous_wave_replica
%  A replica of the visous wave from the viroulet paper in order to get the
%  process working.
    mu1 = tand(20.9);
    mu2 = tand(32.76);
    beta = 0.136;
    L = 8.25e-4;
    

    g=9.81; % m/s^2
    theta = 10;
    gamma = (mu2-tand(theta))/(tand(theta)-mu1);
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
%    crit_Fr = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
    Fr_eq = 1.0; 
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
    u_w_dl = 1+1/Fr_eq-1e-3;
%     u_w_dl = 1.978613807;
    Q1_dl = u_w_dl-1;
    
    opts = odeset('Events',@WaveEndFn);
    [xi_vals,out_vals]=ode15s(@full_system_orig_dl,[0, 50000],[1.001,0]);
    final_y = out_vals(end,:);
    final_h = final_y(1);
    [xi_wave,out_wave,eq_val,~,~] = ode15s(@full_system_orig_dl,[0, 500],final_y,opts);
    non_zero_eq = eq_val(diff(vertcat([0],eq_val))>1e-6);
    eq_index = sum(xi_wave<non_zero_eq(2));
    plot(xi_wave(1:eq_index,1),out_wave(1:eq_index,1))
%     max_h = (-min_h + sqrt(min_h.^2+8*Q1_dl^2.*Fr_eq.^2./min_h))/2;
%     hold on
%     opts_nv = odeset('Events',@NV_WaveEndFn);
%     [xi_nv,out_nv] = ode15s(@non_viscous_system,[0, 50],min_h,opts_nv);
%     plot(xi_nv(:,1),out_nv(:,1))

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

        mu_b = mu1 + (mu2-mu1)*(1-u_w_dl+u_w_dl*h)/(1-u_w_dl+u_w_dl*h+h^(5/2)*gamma);
        n_eq = (tand(theta)-mu_b)./n_coeff;
        
        dy3dx = 1/(2*h)*n^2 + h^(3/2)/Fr_eq^2*R/Q1_dl*n_coeff*(n-n_eq);
%         dy2dx = (g*h^3*cosd(theta))*(tand(theta)-mu_b)/(n_coeff);
        dydx = [dy2dx, dy3dx]';
    end

    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = y(1)-final_h;
        isterminal = 0;
        direction = 0;
    end

    function [position,isterminal,direction] = NV_WaveEndFn(t,y)
        position = y(1)-max_h;
        isterminal = 1;
        direction = 0;
    end
end