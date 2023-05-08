function find_no_deriv_tau0
    mu1_Iv = 0.32;
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();

    Fr_eq = 3; 
    theta = 10;
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    tau0=100;
%     h0=0.25;
    npts = 100;
    tau_min = 0;
    tau_max = 200;
    tau0_vals = linspace(tau_min,tau_max,npts);
    Fr_vals = linspace(0.5,4,npts);
    crit_h=zeros(1,npts);
    h_stop_vals = zeros(1,npts);
    u_eq_vals = zeros(1,npts);
    Iv_vals = zeros(1,npts);
    h0_vals = zeros(1,npts);
    Q1_h0_vals = zeros(1,npts);
    tau0_dl_vals = zeros(1,npts);
    
    fb_vals = zeros(npts,100);
    for j = 1:npts
%         tau0 = tau0_vals(j);
        Fr_eq = Fr_vals(j);

        [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_eq, tau0);
%         [Fr_eq, crit_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0);
        Iv_vals(j)=crit_Iv;
        Fr_vals(j)=Fr_eq;
        h0_vals(j)=h0;
        u_eq = Fr_eq*sqrt(g*cosd(theta)*h0);
        crit_pb = rho_f*g*cosd(theta)*h0;
        
        h_stop = tau0/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);
        
        z_scale = h0;
        v_scale = u_eq;
        t_scale = z_scale/v_scale;
        p_scale = crit_pb;

        g_dl = g*t_scale/v_scale; 
        tau0_dl = tau0/p_scale;
        tau0_dl_vals(j) = tau0_dl;
        
        h_stop_dl = (tau0_dl*rho_f/rho)/(tand(theta)-P*mu1_Iv);
        h_stop_vals(j) = h_stop_dl;
        Q1 = sqrt(g_dl.*cosd(theta));
        u_w = 1 + sqrt(g_dl.*cosd(theta));
        Q1_h0_vals(j) = Q1/u_w;

        n_vals = 100;
        h_vals = linspace(1,Q1/u_w+1e-4,n_vals);
        u_eq_vals(1) = u_eq;
        for i=2:n_vals
            force=force_bal(h_vals(i));
            fb_vals(j,i) = force;
            [~,Iv_alt] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h_vals(i)*h0, tau0,0);
            u_eq_vals(i) = Iv_alt/eta_f/3*(g*cosd(theta)*(rho-rho_f))*(h_vals(i)*h0)^2;
            
%             if force>1e-8
%                 break
%             end
        if h_stop_dl<Q1/u_w
            step=0;
            h_newt = max(h_vals(i),Q1/u_w+1e-4);
            while force>1e-6
                force=force_bal(h_newt);
                fb_deriv = force_bal_deriv(h_newt);
                h_newt = h_newt - force/fb_deriv;
                step=step+1;
            end
            crit_h(j) = h_newt;
        else
            crit_h(j) = Q1/u_w;
        end

        end
    end
    hold on
    C = viridis(4);
    SetPaperSize(8,8)
    
%     plot(h_vals*h0,(u_w-Q1./h_vals)*u_eq,"DisplayName","Wave Velocity","color",C(1,:))
%     plot(h_vals*h0,u_eq_vals,"DisplayName","Equilibrium Velocity","color",C(3,:))
% %     plot([Q1/u_w*h0,Q1/u_w*h0],[0,0.3],"DisplayName","$Q_1/u_w$","color","k","LineStyle","--")
    plot(Fr_vals,h0_vals,"DisplayName","Equilibrium $h$","color",C(1,:))
    plot(Fr_vals,crit_h.*h0_vals,"DisplayName","Minimum $h$","color",C(2,:))
    plot(Fr_vals,Q1_h0_vals.*h0_vals,"DisplayName","$Q_1/u_w$","color",C(3,:))
    plot(Fr_vals,h_stop_vals.*h0_vals,"DisplayName","$h_{stop}$","color",C(4,:))
    plot([0.5,0.5],[0,0.3],"DisplayName","Stability Limit","color","k","LineStyle","--")
    legend("Location","east")
    ylabel("$h$ ($m$)")
    xlabel("$Fr$")
%     ylabel("Force Balance")
%     xlabel("Fraction of $h_{eq}$")

%     xlim([0,10])
    title("$\tau_0 = 20$Pa, $\theta = 10^\circ$")
    exp_graph(gcf,"tau0_Fr_h_stop_vs_min_h_phys.pdf")
    
    function fb = force_bal(h)
        u = -Q1/h + u_w;
        Iv = crit_Iv*u/h^2;
        mu_val = P*mu_Iv_fn(Iv)+tau0_dl*rho_f/rho/h;
        fb = (tand(theta)-sign(u).*mu_val);
    end

    function fb_dash = force_bal_deriv(h)
        u = -Q1/h + u_w;
        u_dash = Q1/h^2;
        Iv = crit_Iv*u/h^2;
        Iv_dash = -2*crit_Iv*u/h^3 + u_dash*crit_Iv/h^2;
        fb_dash = -P*dmudIv_fn(Iv)*Iv_dash + tau0_dl*rho_f/rho/h^2;
    end
end