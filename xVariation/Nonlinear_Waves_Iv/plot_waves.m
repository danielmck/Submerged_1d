function plot_waves
    filename = "high_pe_lambda_12.txt";
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
    if strcmp(wave_type,"full")
        d = record.d(in_table);
        alpha = record.alpha(in_table);
    end
    
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
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);  
    u_eq = u_const.*h0^2;
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    nu_dl = nu/(u_eq*h0);

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

    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    u_w = y(1,1);
    Q1 = y(2,:);
    h = y(3,:);
    u = u_w - Q1./h;
    n = y(4,:);
    m = y(5,:);
    phi = y(6,:)./Q1;
    pb = y(7,:) + rho_dl*g_dl*cosd(theta)*chi.*h;
    pe = pb-h;
    
    [p_max,p_max_ind] = max(pb);
    [p_min,p_min_ind] = min(pb);
    
    [h_max,h_max_ind] = max(h);
    [h_min,h_min_ind] = min(h);
    
    [phi_max,phi_max_ind] = max(phi);
    [phi_min,phi_min_ind] = min(phi);
    
    zeta = 3./(2*alpha_dl.*h) + P/4;
    pp = p_tot_grad_dl.*h-pb;
    D = -2/beta_dl./h.*(pb-h);
    Iv = abs(2*eta_f_dl.*u./h./pp);
    R_w3 = -phi.*rho_f_dl/rho_dl.*D;
    R_w4 = (-P.*chi+zeta).*D - 2*3/alpha_dl./h.*u.*(phi - phi_c./(1+sqrt(Iv)));

    dhdxi = n;
    n_coeff = 1-Q1.^2.*Fr^2./h.^3;
    Iv = 2*eta_f_dl.*abs(u)./h./pp;
    mu_val = pp./(p_tot_grad_dl.*h).*mu_Iv_fn(Iv);
    n_eq = (tand(theta)-sign(u).*mu_val)./n_coeff;
    dQdxi = -P.*D;
    dndxi = 1./(2.*h).*n.^2 + h.^(3/2)/Fr^2./nu_dl./Q1.*n_coeff.*(n-n_eq);
    dmdxi = h./lambda.*u;

    dy6dxi = -R_w3;
    dy7dxi = R_w4./(u-u_w);
    dpbdxi = dy7dxi + rho_dl*g_dl*cosd(theta)*chi.*n;
    
    dpbdxi_scale = dpbdxi/(p_max-p_min);
    dhdxi_scale = n/(h_max-h_min);
    
    plot(xi,pp)
    hold on
    plot(xi,pe)
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end