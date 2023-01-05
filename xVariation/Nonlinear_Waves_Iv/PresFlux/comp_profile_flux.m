function comp_profile_flux
    del_h = 1e-8;
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    
    theta=12;
    Fr=1;
    
    tau0 = 0;
    
    pp_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
    
    [h_eq, eq_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_const = eq_Iv.*pp_grad./eta_f./3;
%     h_eq = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);
    u_eq = u_const*h_eq^2;
    tau0_dl = tau0/(rho_f*g*cosd(theta)*h_eq);
    
    pp_grad_dl = (rho-rho_f)/rho_f/Fr^2;
    eta_f_dl = eq_Iv*pp_grad_dl/3;
    tau0_dl = tau0/(rho_f*g*cosd(theta)*h_eq);
    
    npts = 100;
    h_crit_vals = linspace(1,1.05,npts);
    flux_vals = zeros(npts,npts);
    h1_all = zeros(npts,npts);

    h1 = 0;
    [xi_vals_norm,h_vals_norm,u_vals_norm] = produce_wave(1,h1);
    flux_val = get_flux(xi_vals_norm,h_vals_norm,u_vals_norm);
    if flux_val < 1
        h_lower = 1;
        h_upper = 1.2;
        h_crit = h_upper;
    else
        h_lower = 0.99;
        h_upper = 1;
        h_crit = h_lower;
    end
        
    while (abs(flux_val-1)>1e-6)
        [xi_vals,h_vals,u_vals]= produce_wave(h_crit,h1);
        flux_val = get_flux(xi_vals,h_vals,u_vals);
        if flux_val > 1
            h_upper = h_crit;
        else
            h_lower = h_crit;
        end
        h_crit = (h_upper+h_lower)/2;
    end
    hold on
    C = viridis(4);
    SetPaperSize(8,8)
    xlabel("$\xi$")
    ylabel("$u$")
    title("$\theta = "+num2str(theta)+"$, $Fr = "+num2str(Fr)+"$, $\tau_0 = "+num2str(tau0)+"$Pa")
    plot(xi_vals_norm,u_vals_norm,"DisplayName","Through Equilibrium Flow","color",C(1,:))
    plot(xi_vals,u_vals,"DisplayName","Flux Preserving","color",C(3,:))
    legend("Location","best");
    exp_graph(gcf,"no_tau_static_norm_flux_comp_u.pdf")
    
    function [xi_vals,h_vals,u_vals] = produce_wave(h_crit,h1)
        [~, crit_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h_crit*h_eq, tau0,0);
        
        u_crit = crit_Iv/eq_Iv*h_crit^2;  

        Q1 = (h_crit^(3/2)/Fr);
        u_w = (u_crit*h_crit + Q1)/h_crit;
        
        if tau0_dl == 0
            eq_heights = roots([eq_Iv.*pp_grad_dl./eta_f_dl./2 0 -u_w Q1]);
            h_min = max(eq_heights(eq_heights<(1-1e-6) & eq_heights>0));
        else
            h_min = find_h_min(h_crit,Q1,u_w,theta,tau0_dl,crit_Iv);
        end
        if h1 == 0
            h1 = Q1/u_w+1e-3;
        end
        
        h2 = get_h_plus(h1);
        if (h1>h_min)
            opts = odeset('Events',@WaveEndFn);
            [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,100,201),h1,opts);
            u_vals = get_u(h_vals);
        else
            warning("Initial height is below minimum wave height")
            xi_vals = -1;
            h_vals = -1;
            u_vals = -1;
        end

        function dhdxi = h_deriv_fn(xi,h)
            if (abs(h-h_crit)>1e-6)
                dhdxi = 1/Fr^2.*h.^3.*force_bal(h)./(h.^3/Fr^2-Q1.^2);
            else
                dhdxi = 1/Fr^2.*(3*h^2.*force_bal(h)+h.^3.*force_bal_deriv(h))/(3*h.^2/Fr^2);
            end
        end
        
        function num_val = force_bal_deriv(h)
            u = get_u(h);
            ud = get_u_deriv(h);
            Iv = crit_Iv.*2.*u/h.^2;
            Iv_deriv = crit_Iv.*2.*ud./h.^2-2.*Iv/h;
            num_val = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
        end
        
        function u = get_u(h)
            u = (-Q1 + h.*u_w)./h;
        end
        
        function ud = get_u_deriv(h)
            ud = (-Q1 + h.*u_w)./h;
        end

        function num_val = force_bal(h)
            u = get_u(h);
            Iv = 3.*u.*eta_f_dl./(pp_grad_dl.*h.^2);
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

    function flux = get_flux(xi,h,u)
        flux = 0;
        for i=1:(size(xi)-1)
            flux = flux + (u(i)*h(i)+u(i+1)*h(i+1))*(xi(i+1)-xi(i))/2;
        end
        flux = flux/xi(end);
    end

    
end