function get_min_tau0
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    mu1_Iv = 0.32;
    theta = 12;
    npts=100;
    h_crit_vals = linspace(0.01,1,npts);
    for i=1:npts
        h_crit = 0.1; %h_crit_vals(i);
        Q1 = sqrt(h_crit^3*g*cosd(theta));
        h_stop = 0;
        h_min=1;
        tau0_min = 0;
        tau0_max = 200;
        while (abs(h_stop-0.8*h_min)>1e-7)
            tau0_mid = (tau0_max+tau0_min)/2;
            [~,eq_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, h_crit, tau0_mid,0);
            u_crit = eq_Iv/3/eta_f*h_crit^2*((rho-rho_f)*g*cosd(theta));
            u_w = Q1/h_crit+u_crit;
            h_stop = Q1/u_w;
            h_min = tau0_mid/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);
            if 0.8*h_min-h_stop>1e-7
                tau0_max = tau0_mid;
            elseif 0.8*h_min-h_stop<-1e-7
                tau0_min = tau0_mid;
            end
        end
        tau0=tau0_mid;

        h_start = h_min+1e-6;
        h_max = (-h_start + sqrt(h_start.^2+8*Q1^2./(h_start*g*cosd(theta))))/2;
        opts = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',@WaveEndFn);
        [xi_vals,out_vals]=ode15s(@viscous_syst,[0,1e3],[h_start,0],opts);
        h_ave = out_vals(end)/xi_vals(end);
    end
    
    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = y(2)-h_max;
        isterminal = 1;
        direction = 1;
    end
    
    function dydxi = viscous_syst(t,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            h = y(1);
            u = (-Q1 + h.*u_w)./h;
            m = y(2);
            
            Iv = 3*eta_f.*abs(u)/(rho-rho_f)/g/cosd(theta)/h.^2;
            
            denom = (g*cosd(theta)*h.^3-Q1.^2);
            if (abs(denom)>1e-6)
                fb_val = g*sind(theta)*h-sign(u).*mu_Iv_fn(Iv).*P*g*cosd(theta)*h-sign(u).*tau0/rho/h;
                dhdxi = h.^2.*fb_val./denom;
            else
                ud = Q1/h^2;
                Iv_deriv = 3*eta_f.*ud/(rho-rho_f)/g/cosd(theta)/h.^2-6*eta_f.*abs(u)/(rho-rho_f)/g/cosd(theta)/h.^3;
                fb_deriv = g*sind(theta)-dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0/rho./h.^2;
                dhdxi = fb_deriv/(3*h.^2*g*cosd(theta));
            end

            dmdxi = h;
            dydxi = [dhdxi,dmdxi]';
        end
end