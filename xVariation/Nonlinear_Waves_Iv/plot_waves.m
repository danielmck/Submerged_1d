function plot_waves
    filename = "master_wave_no_pe.txt"; %long_low_pe
    master_file = load("Results/"+filename);
    xi = master_file(1,:);
    y = master_file(2:end,:);
    record = readtable('Results/wave_record.csv');

    in_table = strcmp(record.Name, filename);
    wave_type = record.wave_type(in_table);
    theta = record.theta(in_table); 
    lambda = record.lambda(in_table);
    Fr = record.Fr(in_table);
    nu = record.nu(in_table);
    tau0 = record.tau0(in_table);
    if strcmp(wave_type,"full")
        full_model=true;
        d = record.d(in_table);
        alpha = record.alpha(in_table);
    else
        full_model=false;
        alpha=0;
        d=0;
    end
    s=0;
    
    
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;
    
    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    
%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    
    rho_p = 2500;
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    s_c = 1-rho/(rho-rho_f)*tand(theta)/mu1_Iv;
    
    if tau0 == 0
        crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f,s);
        pp_eq_grad = (rho_p-rho_f)*g*phi_c*cosd(theta);
        u_const = crit_Iv/eta_f/2*pp_eq_grad*(1-s);
        h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
    else
        [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    end
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    nu_dl = nu/(u_eq*h0);
  
    h_stop = tau0/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);
    

    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 

    u_eq_dl = u_eq/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;
    tau0_dl = tau0/p_scale;
    h_stop_dl = h_stop/z_scale;

    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    u_w = y(1,1);
    Q1 = y(2,:);
    h = y(3,:);
    u = u_w - Q1./h;
    n = y(4,:);
    m = y(5,:);
    
    if full_model
        phi = y(6,:)./Q1;
        pb = y(7,:) + rho_dl*g_dl*cosd(theta)*chi.*h;
        pe = pb-h;
    else
        pb=h;
    end
    
    Fr_vals = Fr.*u./sqrt(h);
    
    h_min = roots([1,0,-u_w,-Q1(1)]);
    h_static = Q1(1)/u_w;
    
    [p_max,p_max_ind] = max(pb);
    [p_min,p_min_ind] = min(pb);
    
    [h_max,h_max_ind] = max(h);
    [h_min,h_min_ind] = min(h);
    
    if full_model
        [phi_max,phi_max_ind] = max(phi);
        [phi_min,phi_min_ind] = min(phi);
    end
    
    pp = p_tot_grad_dl.*h-pb;
    D = -2/beta_dl./h.*(pb-h);
    Iv = abs(2*eta_f_dl.*u./h./pp);
    
    if full_model
        zeta = 3./(2*alpha_dl.*h) + P/4;
        tan_psi = phi - phi_c./(1+sqrt(Iv));
        R_w3 = -phi.*rho_f_dl/rho_dl.*D;
        R_w4 = (-P.*chi+zeta).*D - 2*3/alpha_dl./h.*u.*(tan_psi);
    end
    
    Fr_equi = zeros(size(h));
    Iv_equi = zeros(size(h));
    for i=1:size(h,2)
        [Fr_equi(i),Iv_equi(i)] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h(i)*h0, tau0);
    end

    dhdxi = n;
    n_coeff = 1-Q1.^2.*Fr^2./h.^3;
    Iv = 2*eta_f_dl.*abs(u)./h./pp;
    mu_val = pp./(p_tot_grad_dl.*h).*mu_Iv_fn(Iv)+tau0_dl*rho_f/rho./h;
    force_bal = tand(theta)-sign(u).*mu_val;
    
    mu_val_min = (rho-rho_f)/rho.*mu1_Iv+tau0_dl*rho_f/rho./h;
    force_bal_max = tand(theta)-sign(u).*mu_val_min;
    
    n_eq = (force_bal)./n_coeff;
    n_diff = n_coeff.*n - force_bal;
    dQdxi = -P.*D;
    dndxi = 1./(2.*h).*n.^2 + h.^(3/2)/Fr^2./nu_dl./Q1.*n_coeff.*(n-n_eq);
    dn_term1 = 1./(2.*h).*n.^2;
    dn_term2 = h.^(3/2)/Fr^2./nu_dl./Q1.*n_coeff.*(n-n_eq);
    dmdxi = h./lambda.*u;
    
    if full_model
        dy6dxi = -R_w3;
        dy7dxi = R_w4./(u-u_w);
        dpbdxi = dy7dxi + rho_dl*g_dl*cosd(theta)*chi.*n;

        dpbdxi_scale = dpbdxi/(p_max-p_min);
        dhdxi_scale = n/(h_max-h_min);
    end
    %%
    C = viridis(3);
%     SetPaperSize(10,10)
    hold on
%     plot(linspace(0.5,1),get_force_bal(linspace(0.5,1)))
    plot(xi(xi<5),h(xi<5)-Q1(xi<5)/u_w, "DisplayName", "Waveform","color",C(1,:))
%     plot(Iv_equi,h, "DisplayName", "Equilibrium Curve","color",C(2,:))
%     plot(xi,dn_term2, "DisplayName", "$\frac{d y_7}{d \xi}$","color",C(2,:))
%     plot(xi,h_stop_dl*ones(size(xi)), "DisplayName", "Minimum $h$","color",C(3,:))
%     plot(xi,dpbdxi(xi>20), "DisplayName", "$\frac{d p_b}{d \xi}$","color",C(3,:))
    
    xlabel("$Fr$")
    ylabel("$h$")
%     legend("Location","best")
    title("$\theta="+num2str(theta)+"$, $\tau_0="+num2str(tau0)+"$, $s="+num2str(s)+"$, No $p_e$ case")
%     plot(xi,force_bal)
%     exp_graph(gcf,"no_pe_tau0_20_u_stop_h_Fr_zoom.pdf")
      
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function fb = get_force_bal(h_fb)
        u_fb = u_w - Q1(1)./h_fb;
        Iv_fb = crit_Iv.*u_fb./h_fb.^2;
        mu_fb = (rho-rho_f)/rho.*mu_Iv_fn(Iv_fb)+tau0_dl*rho_f/rho./h_fb;
        fb = tand(theta)-sign(u_fb).*mu_fb;
    end
end