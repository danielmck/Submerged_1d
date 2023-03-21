function [xi_vals, h_vals, Q1, u_w,flux] = depo_wave_ode(theta, Fr, tau0)
% Creates a wave that preserves the flux of the inputted conditions via the
% ODE method which can be used to solve for the boundary problem.
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    mu1_Iv = 0.32;
    
    [h_eq, eq_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,0);
%     pp_grad = (rho_p-rho_f)*g*phi_c*cosd(theta_in);
%     u_const = eq_Iv.*pp_grad./eta_f./3;
%     h_eq = ((Fr*sqrt(g*cosd(theta_in)))./u_const)^(2/3);
%     u_eq = u_const*h_eq^2;
    tau0_dl = tau0/(rho_f*g*cosd(theta)*h_eq);
    crit_point= 1;
    h_jump = 0.1;
    h_max = 2.5;
    xi_vals = -1;
    h_vals = -1;
    last_flux = 1;
    flux=1;
    while crit_point<h_max && (last_flux-1)*(flux-1)>=0
        last_flux=flux;
        [flux,xi_vals,h_vals,Q1,u_w] = wave_flux(crit_point);
        crit_point=crit_point+h_jump;
    end
    if (last_flux-1)*(flux-1)<=0
        crit_point=crit_point-h_jump;
        h_top = crit_point;
        h_bot = crit_point-h_jump;
        fsign = sign(flux-1);
        while abs(flux-1)>1e-4
            crit_point = (h_top+h_bot)/2;
            [flux,xi_vals,h_vals,Q1,u_w] = wave_flux(crit_point);
            if (flux-1)*fsign > 0
                h_top = crit_point;
            else
                h_bot = crit_point;
            end
        end
    end
    
    function [flux,xi_vals,h_vals,Q1,u_w] = wave_flux(h_crit)
        
        [~, crit_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h_crit*h_eq, tau0,0);
        u_crit = crit_Iv/eq_Iv*h_crit^2;
        
        Q1 = (h_crit^(3/2)/Fr);
        u_w = (u_crit*h_crit + Q1)/h_crit;

        h_stop_dl = (tau0_dl*rho_f/rho)/(tand(theta)-P*mu1_Iv);

        if (h_stop_dl > Q1/u_w)
            h1=Q1/u_w+0.001;
            h2 = get_h_plus(h1);
            opts = odeset('Events',@WaveEndFn);
            [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,1000,201),h1,opts);
            u_vals = get_u(h_vals);
            
            flux = get_flux(xi_vals,h_vals,u_vals);
        else
            flux=0;
            xi_vals=0;
            h_vals=0;
        end
    
        function dhdxi = h_deriv_fn(xi,h)
%             if (abs(h-h_crit)>1e-6)
            dhdxi = 1/Fr^2.*h.^3.*force_bal(h)./(h.^3/Fr^2-Q1.^2);
%             else
%                 dhdxi = 1/Fr^2.*h.^3.*(force_bal_deriv(h))/(3*h.^2/Fr^2);
%             end
        end

        function num_val = force_bal_deriv(h)
            u = get_u(h);
            ud = get_u_deriv(h);
            Iv = eq_Iv.*u/h.^2;
            Iv_deriv = eq_Iv.*2.*ud./h.^2-2.*Iv/h;
            num_val = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
        end

        function u = get_u(h)
            u = (-Q1 + h.*u_w)./h;
        end

        function ud = get_u_deriv(h)
            ud = Q1/h^2;
        end

        function num_val = force_bal(h)
            u = get_u(h);
            Iv = eq_Iv.*u/h.^2;
            num_val = tand(theta)-mu_Iv_fn(Iv).*(rho-rho_f)/rho-tau0_dl*rho_f/rho/h;
        end

        function h_plus = get_h_plus(h_minus)
            h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(h_minus/Fr^2)))/2;
        end

        function [position,isterminal,direction] = WaveEndFn(t,y)
            position = y - h2;
            isterminal = 1;
            direction = 0;
        end
    end
end