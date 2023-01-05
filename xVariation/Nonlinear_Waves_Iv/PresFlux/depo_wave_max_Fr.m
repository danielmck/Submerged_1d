function depo_wave_max_Fr
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    mu1_Iv = 0.32;
    
    nvals=100;
    theta_vals = linspace(9,15,nvals);
    tau0_vals = linspace(0,100,nvals);
    max_Fr = zeros(nvals,nvals);
    max_Fr(1,2) = 3;
    for m = 1:nvals
        theta = theta_vals(m);
        for n=1:nvals
            tau0=tau0_vals(n);
            max_Fr(m,n) = single_contours(theta, tau0,false,3,50);
        end
    end
    save("max_Fr_static.txt","max_Fr","-ascii")
    SetPaperSize(10,10)
    contour(theta_vals,tau0_vals,max_Fr')
    ylabel("$\tau_0$ (Pa)")
    xlabel("$\theta$")
%     clabel(C,h,'manual','FontSize',7)
    cb=colorbar();
    cb.Label.String = 'Maximum $Fr$ for static regions';
    exp_graph(gcf,"max_Fr_depo_wave.pdf")

    function max_Fr = single_contours(theta_in, tau0_in,dl,h_crit_max,npts)
        Fr_in=1;
        [h_eq, eq_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in,dl);
        pp_grad = (rho_p-rho_f)*g*phi_c*cosd(theta_in);
        u_const = eq_Iv.*pp_grad./eta_f./3;
    %     h_eq = ((Fr*sqrt(g*cosd(theta_in)))./u_const)^(2/3);
        u_eq = u_const*h_eq^2;
        if dl
            tau0_dl = tau0_in;
            tau0_d = tau0_dl*(rho_f*g*cosd(theta_in)*h_eq);
        else
            tau0_dl = tau0_in/(rho_f*g*cosd(theta_in)*h_eq);
            tau0_d = tau0_in;
        end
        pp_grad_dl = (rho-rho_f)/rho_f/Fr_in^2;
        eta_f_dl = eq_Iv*pp_grad_dl/3;

        h_c_vals = linspace(1.5,h_crit_max,npts);
        max_flux = -1;
        for j=1:npts
            h_crit = h_c_vals(j);
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h_eq, tau0_d,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2;  

            Q1 = (h_crit^(3/2)/Fr_in);
            u_w = (u_crit*h_crit + Q1)/h_crit;

            h_stop_dl = (tau0_dl*rho_f/rho)/(tand(theta)-P*mu1_Iv);

            if (h_stop_dl > Q1/u_w)
                h1=Q1/u_w+0.001;
                h2 = get_h_plus(h1);
                opts = odeset('Events',@WaveEndFn);
                [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,100,201),h1,opts);
                u_vals = get_u(h_vals);
                max_flux = max(max_flux,get_flux(xi_vals,h_vals,u_vals));
            end
            if max_flux>0
                max_Fr = get_Fr_from_flux(max_flux*u_eq*h_eq,theta,tau0);
            else
                max_Fr=0;
            end
        end
        function dhdxi = h_deriv_fn(xi,h)
            if (abs(h-h_crit)>1e-6)
                dhdxi = 1/Fr_in^2.*h.^3.*force_bal(h)./(h.^3/Fr_in^2-Q1.^2);
            else
                dhdxi = 1/Fr_in^2.*h.^3.*(force_bal_deriv(h))/(3*h.^2/Fr_in^2);
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
            h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(h_minus/Fr_in^2)))/2;
        end
        
        function [position,isterminal,direction] = WaveEndFn(t,y)
            position = y - h2;
            isterminal = 1;
            direction = 0;
        end
    end
end