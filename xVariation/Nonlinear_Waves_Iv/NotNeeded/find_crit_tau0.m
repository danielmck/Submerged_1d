function find_crit_tau0
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    
    eta_f = 0.0010016; % Pa s
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^

%     eta_f = 1.18e-5;
%     rho_f = 1;

%     Fr_eq = 3; 
%     theta = 10;
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
%     tau0=100;
%     h0=0.25;
    npts = 100;
%     tau0_vals = linspace(tau_min,tau_max,npts);
    Fr_vals = linspace(0.5,5,npts);
    theta_vals = linspace(9.5,12,npts);

%     h_stop_vals = zeros(1,npts);
    Iv_vals = zeros(npts,npts);
    h0_vals = zeros(npts,npts);
    Q1_h0_vals = zeros(1,npts);
    tau0_vals = zeros(1,npts);
    
    del_Iv = 1e-8;
    max_tol = 1e-7;
    for i=1:npts
        Fr_eq = Fr_vals(i);
        Q1 = 1/Fr_eq;
        u_w = 1+1/Fr_eq;
        
        for j=1:npts
            theta = theta_vals(j);
            tau_0_scale = (tand(theta)-(rho-rho_f)/rho*mu1_Iv)*rho/rho_f;
            tau0_dl = Q1/u_w*tau_0_scale;
            resid = 1;   
            Iv = newt_solve_crit_Iv(theta,rho_p,rho_f);
            crit_mu = rho/(rho-rho_f)*tand(theta)-tau0_dl*rho_f/(rho-rho_f);
            if crit_mu>mu1_Iv
                while (abs(resid)>max_tol)
                    resid = mu_Iv_fn(Iv)-crit_mu;
                    if Iv > 1e-7
                        force_deriv = dmudIv_fn(Iv);
                    else
                        force_deriv = ((mu_Iv_fn(Iv+del_Iv)-crit_mu)-resid)/del_Iv;
                    end
                    Iv = Iv-resid/force_deriv;
                end
                h0 = (Fr_eq/Iv/(rho-rho_f)/sqrt(g*cosd(theta))*3*eta_f)^(2/3);
                tau0 = tau0_dl*rho_f*g*cosd(theta)*h0;
            else
                warning("The yield stress does not allow for flow on this slope");
                Iv = 0;
                h0 = -1;
                tau0=-1;
            end
            Iv_vals(i,j) = Iv;
            h0_vals(i,j) = h0;
            tau0_vals(i,j) = tau0;
        end
    end
    SetPaperSize(8,8)
    C = viridis(4);
    contour(Fr_vals,theta_vals,tau0_vals') %./(rho_f*g*cosd(theta)*h0_vals)
%     contour(Fr_vals,theta_vals,tau0_vals',[4,5,6,7,8,9],'ShowText','on')
    hold on
%     contour(Fr_vals,theta_vals,tau0_vals',[4.5,5.5,6.5,7.5,8.5])
    c = colorbar;
    c.Label.String = '$\tau_0$ (Pa)'; %"$\tau_0^{\ast}$";%;
    xlabel("$Fr$ at $h_$")
    ylabel("$\theta$")
    exp_graph(gcf,"theta_Fr_tau0_crit.pdf")
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end