function int_to_wave    
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    theta = 10;
    
    alpha = 1e-5;
    rho_p = 2500;
    d = 1e-4;

    rho_f = 1000;
    eta_f = 0.0010016; % Pa s

%     rho_f = 1; % kg/m^3 
%     eta_f = 1.18e-5; % Pa s
    
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    crit_phi = phi_c./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
    
    crit_Fr = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
    Fr = 0.6;
    h0 = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3);
    
    crit_u = u_const*h0^2;
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    nu = eta_f/rho_f;
    
    z_scale = h0;
    v_scale = crit_u;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 

    crit_u_dl = crit_u/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;
    nu_dl = nu/v_scale/z_scale;

    P = (rho_dl-rho_f_dl)/rho_dl;
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    chi = (rho_f+3*rho)/(4*rho);
    
    u_w = 1+sqrt(g_dl*cosd(theta))-1e-2;
    
    init_h = 1+1e-1;
    init_u = u_w+(1-u_w)/init_h;
    init_pb = 1;
    init_phi = crit_phi;
    init_n = 0;
    init_D = -2/beta_dl/init_h*(init_pb-rho_f_dl*g_dl*cosd(theta)*init_h);
    init_w1 = init_h*(init_u-u_w);
    init_w3 = -init_h*(init_u-u_w)*init_n + init_h*P*init_D;
    init_w4 = init_phi*init_h*(init_u-u_w);
    init_w5 = init_pb - rho_dl*g_dl*cosd(theta)*chi.*init_h;
    init_vals = [init_w1, init_h, init_w3, init_w4, init_w5];
    
    [xi_vals,out_vals]=ode15s(@full_system,[0, 1000],init_vals);
    h_vals = out_vals(:,2);
    u_vals = (out_vals(:,1)+u_w.*h_vals)./h_vals;
    pb_vals = out_vals(:,5) + rho_dl*g_dl*cosd(theta)*chi.*h_vals;
    D_vals = -2/beta_dl./h_vals.*(pb_vals-rho_f_dl*g_dl*cosd(theta).*h_vals);
    n_vals = (P.*D_vals-out_vals(:,3)./h_vals)./(u_vals-u_w);
    plot(xi_vals,out_vals(:,2));

    function dvecdx = full_system(x,y)
        h = y(2);
        u = (y(1)+u_w.*h)./h;
        phi = y(4)./y(1);
        pb = y(5) + rho_dl*g_dl*cosd(theta)*chi.*h;
        
        D = -2/beta_dl/h*(pb-rho_f_dl*g_dl*cosd(theta)*h);
        
        n = (P*D-y(3)/h)/(u-u_w);
        
%         pb = h;
        zeta = 3/(2*alpha_dl*h) + P/4;
        p_p = p_tot_grad_dl*h-pb;
        
        Iv = abs(2*eta_f_dl*u/h/p_p);
        R_w1 = P*D;
        R_w3 = g_dl*sind(theta)*h - sign(u)*mu_Iv_fn(Iv)/rho_dl*p_p;
        R_w4 = -phi*rho_f_dl/rho_dl*D;
        R_w5 = (-P*rho_dl*g_dl*cosd(theta)*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));
        
        dy1dx = R_w1;
        dy2dx = n;
        dy3dx = -(R_w3-(2*u-u_w)*P*D - (g_dl*cosd(theta)*h-(u-u_w)^2)*n)/nu_dl*rho_dl;
        dy4dx = R_w4;
        dy5dx = R_w5/(u-u_w);
        dvecdx = [dy1dx dy2dx dy3dx dy4dx dy5dx]';
%         dvecdx = [dy1dx dy2dx];
    end
    
    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end