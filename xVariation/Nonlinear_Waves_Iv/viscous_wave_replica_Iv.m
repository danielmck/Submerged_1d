function [xi_out, y_out, uw_out] = viscous_wave_replica_Iv(theta, rho_f, rho_p, eta_f, Fr_eq, lambda)
%     Finds a solution for the viscous problem withut pe that is close to
%     the given value lambda. Is used in conjunction with
%     viscous_Iv_bvp_from_ode to find a wave of a specific length.
    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2

%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    crit_mu = rho/(rho-rho_f)*tand(theta);
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
%     Fr_eq = h0*beta/gamma/L;
    h0 = ((Fr_eq*sqrt(g*cosd(theta)))./u_const)^(2/3);
    
    u_eq = u_const.*h0^2;
    nu = 3/4*crit_mu*eta_f/crit_Iv/rho;
    
    z_scale = h0;
    v_scale = u_eq;
    t_scale = z_scale/v_scale;
    

    nu_dl = nu/(v_scale*z_scale);
    R = u_eq*h0/nu;
    
%     u_w = 1+sqrt(g_dl*cosd(theta));
%     Q1 = u_w-1;
    
%     Q1 = sqrt(h0^3*g*cosd(theta));
    lambda_ratio = -1;
    uw_out = 0;
    xi_out = [];
    y_out = [];
    del_uw = [1e-5,5e-5,1e-4,1.2e-4,5e-4];
    for i = 1:size(del_uw,2)
        u_w_dl = 1+1/Fr_eq-del_uw(i);
        Q1_dl = u_w_dl-1;

    %     u_w = u_eq*u_w_dl;
    %     Q1 = h0*(u_w-u_eq);

        opts1=odeset('AbsTol',1e-6,'RelTol',1e-6);
        opts2 = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',@WaveEndFn);

        [xi_vals,out_vals]=ode15s(@full_system_orig_dl,[0,5e3],[Q1_dl,1.01,0],opts1);
        y_end = out_vals(end,:);
        [xi_wave,out_wave,eq_val,~,~] = ode15s(@full_system_orig_dl,[0, 100],out_vals(end,:),opts2);
        non_zero_eq = eq_val(eq_val>1e-6);
        if size(non_zero_eq,1) > 2
            lambda_i = non_zero_eq(1);
            ratio = max(lambda/lambda_i,lambda_i/lambda);
            if ((lambda_ratio<0)||(ratio<lambda_ratio))
                lambda_ratio = ratio;
                xi_out = xi_wave(1:sum(xi_wave<lambda_i)+1);
                y_out = out_wave(1:sum(xi_wave<lambda_i)+1,:);
                uw_out = u_w_dl;
            end
        end
    end

    if (size(y_out,1) == 0)
        error("No suitable waves were produced")
    end
%     max_h = (-min_h + sqrt(min_h.^2+8*Q1_dl^2.*Fr_eq.^2./min_h))/2;
%     hold on
%     opts_nv = odeset('Events',@NV_WaveEndFn);
%     [xi_nv,out_nv] = ode15s(@non_viscous_system,[0, 50],min_h,opts_nv);
    % plot(xi_nv(:,1),out_nv(:,1))
    plot(xi_out,y_out(:,2))

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
        Q1_dl = y(1);
        h = y(2);
        u = u_w_dl - Q1_dl/h;
        n = y(3);

        dy1dx = 0;
        dy2dx = n;
        n_coeff = h/Fr_eq^2-Q1_dl^2/h^2;
        Fr = Fr_eq*abs(u)/h;
        Iv = crit_Iv*abs(u)/h^2;
        force_bal = h*(tand(theta)-sign(u).*P*mu_Iv_fn(Iv))/Fr_eq^2;
        n_eq = (force_bal)./n_coeff;
        dy3dx = n^2/h+h/Q1_dl/nu_dl*(n_coeff*n-force_bal);
%         dy2dx = (g*h^3*cosd(theta))*(tand(theta)-mu_b)/(n_coeff);
        n_coeff_old = 1-Q1_dl^2.*Fr_eq^2/h^3;
        n_eq_old = (tand(theta)-P*mu_Iv_fn(Iv))./n_coeff_old;
        dy3dx_old = 1/(2*h)*n^2 + h^(3/2)/Fr_eq^2*R/Q1_dl*n_coeff_old*(n-n_eq_old);
        dydx = [dy1dx dy2dx, dy3dx]';
    end

    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = y(2)-y_end(2);
        isterminal = 0;
        direction = 2*(y_end(3)>0)-1;
    end

    function [position,isterminal,direction] = NV_WaveEndFn(t,y)
        position = y(1)-max_h;
        isterminal = 1;
        direction = 0;
    end
end