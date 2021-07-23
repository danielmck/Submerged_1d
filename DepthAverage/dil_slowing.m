function dil_slowing
    N = 200; % number of discretisation points in z for each of tau, u
    h = 4e-1; % layer height (m)
    d=1.43e-3; % grain diameter (m)

    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;

    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 10^7;
    t_step = 1;

    phi_c=0.6; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^3
    theta = 9; % deg

    alpha = 0.001; % 1/Pa
    buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

    v_scale = sqrt(g.*h);
    p_scale = rho_f.*g.*h;
    t_scale = sqrt(h./g);
    z_scale = h;
    density_ratio = rho_p./rho_f;
    rho = phi_c*density_ratio+1-phi_c;

    d_dl = d/z_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    buoyancy_dl = buoyancy/p_scale*z_scale;

    beta = 150*phi_c^2*eta_f_dl/((1-phi_c).^3.*d_dl^2);
    p_b = (density_ratio-1)*phi_c*cosd(theta);
    time_vals = (0:500)*t_step;
    opts=odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','on');

    fname = "dil_9deg_depth_averaged_slowing_unmod_shear_phi.txt";
    [~,vec]=ode15s(@phi_slowing,time_vals,0.6-0.007,opts);
    size(vec)
    plot(vec)
    save(fname, 'vec','-ascii')

    function dudt=depth_ave_slowing(t,u)
        mu = mu1_I + (mu2_I-mu1_I)/(1+I_0/(sqrt(2)*d_dl*u/sqrt(cosd(theta)*phi_c)));
        force_bal = sind(theta) - mu*p_b/(rho);
%         A = 1+37.5 * mu /((1-phi_c^2)*sqrt(phi_c))*(eta_f_dl)^(3/2)/(rho*sqrt(density_ratio*cosd(theta))*d^2*sqrt(u));
        A = 1+1/(4*rho)*beta*sqrt(2*eta_f_dl/(cosd(theta)*density_ratio*u));
        dudt = force_bal/A;
    end

    function dudt=phi_slowing(t,phi)
        u = (phi_c-phi).^2/phi_c*cosd(theta)/eta_f_dl;
        mu = mu1_I + (mu2_I-mu1_I)/(1+I_0/(sqrt(2)*d_dl*u/sqrt(cosd(theta)*phi_c)));
        force_bal = sind(theta) - mu*p_b/(rho);
%         A = 1+37.5 * mu /((1-phi_c^2)*sqrt(phi_c))*(eta_f_dl)^(3/2)/(rho*sqrt(density_ratio*cosd(theta))*d^2*sqrt(u));
        A_phi = (mu.*beta-2*cosd(theta)*(phi_c-phi))/phi_c;
        dudt = force_bal/A_phi;
    end
end