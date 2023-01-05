function single_roll_wave_flux
    del_h = 1e-8;
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    dl=false;
    
    theta=12;
    Fr=1.4;
    
    tau0=20;
    h_crit_max = 1.05;
    filename = "h_crit_h1_1_flux_cont_Fr.pdf";
    singlename = "h_crit_h1_tau0_20_Fr1_4.pdf";
%     param_flux_contour(1,5,[1,theta,tau0],[2,theta,tau0],dl,filename)
    clf;
    single_flux_contours(Fr,theta,tau0,0,h_crit_max,singlename)
    
    function single_flux_contours(Fr_val,theta_val,tau0_val,dl,h_crit_max,fname)
        npts=100;
        SetPaperSize(8,8)
        hold on
        [h_crit_vals,h1_all,flux_vals,h_min_vals] = single_contours(Fr_val, theta_val, tau0_val,dl,h_crit_max,100);
        [C,h]=contour(h_crit_vals.*ones(npts,1),h1_all',flux_vals','HandleVisibility','off');
        plot(h_crit_vals,h_min_vals,"DisplayName","$Q_1/u_w$","color","k")
% [0.8,0.9,1,1.1,1.2,1.3,1.4]
        clabel(C,h,'manual','FontSize',7)
        cb=colorbar();
        cb.Label.String = 'Flux $Q$';
        legend("Location","best");
        xlabel("$h_{crit}$")
        ylabel("$h_1$")
        title("$Fr = "+num2str(Fr_val)+"$, $\theta = "+num2str(theta_val)+"$, $\tau_0 = "+num2str(tau0_val)+"$ Pa")
    %     plot(h1_vals,flux_vals)
        exp_graph(gcf,fname)
    end
    
    function param_flux_contour(param_alt,nlines,params_min,params_max,dl,fname)
        names=["Fr","\theta","\tau_0^{\ast}"];
        unit = ["$","$","$"];
        Col = viridis(nlines+1); %[1,1.1,1.2,1.3,1.4,1.5,1.6]
        
        leg_str = "$"+names(param_alt)+"=";
        title_str = "Flux preserving waves for $";
        for pm=1:3
            if (pm~=param_alt)
                title_str = title_str+names(pm)+"="+num2str(params_min(pm))+unit(pm);
                if (pm ~= min(3,5-param_alt)) 
                    title_str=title_str+", $";
                end
            end
        end
%         title_str = title_str(1:end-3);
        contours = find_flux_contours(nlines,params_min,params_max,dl,100);
        vals = linspace(params_min(param_alt),params_max(param_alt),nlines);
        SetPaperSize(8,8)
        hold on
        for m=1:nlines
            cont = contours{1,m};
            [~,ml] = max(cont(1,:));
            cont = cont(:,1:ml);
            plot(cont(1,cont(2,:)<1),cont(2,cont(2,:)<1),"DisplayName",leg_str+num2str(vals(m)+unit(param_alt)),"color",Col(m,:))
        end

        legend("Location","best");
        xlabel("$h_{crit}$")
        ylabel("$h_1$")
        title(title_str)
    %     plot(h1_vals,flux_vals) 
        exp_graph(gcf,fname)
    end
    
    function contours = find_flux_contours(nlines,params_min,params_max,dl,npts)
        min_cell = num2cell(params_min);
        [Fr_min,theta_min, tau0_min] = min_cell{:};
        if nlines>1
            max_cell = num2cell(params_max);
            [Fr_max,theta_max, tau0_max] = max_cell{:};
            Fr_vals = linspace(Fr_min,Fr_max,nlines);
            theta_vals = linspace(theta_min,theta_max,nlines);
            tau0_vals = linspace(tau0_min,tau0_max,nlines);
        else
            Fr_vals = [Fr_min];
            theta_vals = [theta_min];
            tau0_vals = [tau0_min];
        end
        contours = cell(1,nlines);
        
        for l=1:nlines
            Fr_in = Fr_vals(l);
            theta_in = theta_vals(l);
            tau0_in = tau0_vals(l);
            [h_crit_vals,h1_vals,flux_vals,~] = single_contours(Fr_in, theta_in, tau0_in,dl,h_crit_max,npts);
        %     hold on

%             [X,Y] = meshgrid(h_crit_vals,h1_vals);
            hold on
            [contours{1,l},~]=contour(h_crit_vals.*ones(npts,1),h1_vals',flux_vals',[1,1]);
            hold off
            clf;
        end
    end

    function [h_c_vals,h1_all,flux_vals,h_min_vals] = single_contours(Fr_in, theta_in, tau0_in,dl,h_crit_max,npts)
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
        
        h_c_vals = linspace(1,h_crit_max,npts);
        h_min_vals = zeros(1,npts);
        flux_vals = zeros(npts,npts);
        h1_all = zeros(npts,npts);
        for j=1:npts
            h_crit = h_c_vals(j);
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h_eq, tau0_d,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2;  

            Q1 = (h_crit^(3/2)/Fr_in);
            u_w = (u_crit*h_crit + Q1)/h_crit;

            if tau0_dl == 0
                eq_heights = roots([eq_Iv.*pp_grad_dl./eta_f_dl./2 0 -u_w Q1]);
                h_min = max(eq_heights(eq_heights<(1-1e-6) & eq_heights>0));
            else
                h_min = find_h_min(h_crit,Q1,u_w,theta_in,tau0_dl,crit_Iv);
            end
            h_min_vals(j) = h_min;
            h_span = 1-h_min;
            h1_vals = linspace(h_min+0.001,0.999,npts);

            h1_all(j,:) = h1_vals;

            for k = 1:npts
                h1=h1_vals(k);
                h2 = get_h_plus(h1);
                if (h1>h_min)
                    opts = odeset('Events',@WaveEndFn);
                    [xi_vals,h_vals]=ode15s(@h_deriv_fn,linspace(0,100,201),h1,opts);
                    u_vals = get_u(h_vals);
                    flux_vals(j,k) = get_flux(xi_vals,h_vals,u_vals);
                    n_vals = zeros(size(xi_vals));
                    for i = 1:size(xi_vals,1)
                        n_vals(i) = h_deriv_fn(xi_vals(i),h_vals(i));
                    end
                else
                    warning("Initial height is below minimum wave height")
                end
            end
        end
 

    function [position,isterminal,direction] = WaveEndFn(t,y)
        position = y - h2;
        isterminal = 1;
        direction = 0;
    end
    
    function h_plus = get_h_plus(h_minus)
        h_plus = (-h_minus + sqrt(h_minus.^2+8*Q1^2./(h_minus/Fr_in^2)))/2;
    end
    
    function dhdxi = h_deriv_fn(xi,h)
        if (abs(h-h_crit)>1e-6)
            dhdxi = 1./Fr_in^2.*h.^3.*force_bal(h)./(h.^3/Fr_in^2-Q1.^2);
        else
            dhdxi = (1./Fr_in^2.*(3*h.^2.*force_bal(h)+h.^3.*(force_bal(h+del_h)-force_bal(h))./del_h))./(3.*h.^2./Fr_in^2);
        end
    end

    function u = get_u(h)
        u = (-Q1 + h.*u_w)./h;
    end

    function num_val = force_bal(h)
        u = get_u(h);
        Iv = 3.*u.*eta_f_dl./(pp_grad_dl.*h.^2);
        num_val = tand(theta_in)-mu_Iv_fn(Iv).*(rho-rho_f)/rho-tau0_dl*rho_f/rho./h;
    end
    
    end
    
end