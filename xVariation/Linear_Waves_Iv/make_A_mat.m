function A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr,crit_Iv,nu)
% Makes the A matrix for the linear growth of perturbations

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;


    n_cases = size(theta,2);
    phase = zeros(4,2);


    if ~exist('crit_Iv','var')
        crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);
    end
    if ~exist('nu','var')
        nu=0;
    end
    
    crit_phi = phi_c./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);

    h = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3); % layer height (m)

    p_p = (rho_p-rho_f)*g*phi_c*cosd(theta)*h;
    crit_pb = rho_f*g*cosd(theta)*h;
    crit_u = crit_Iv/eta_f/2*p_p*h;

    v_scale = crit_u;
    p_scale = crit_pb;
    z_scale = h;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale;
    d_dl = d/z_scale;

    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    p_p_dl = p_p/p_scale;
    nu_dl = nu/(v_scale*z_scale);

    chi_dl = (rho_f_dl+3*rho_dl)/(4*rho_dl);
    P = (rho_dl-rho_f_dl)/rho_dl;
    zeta = 3/(2*alpha_dl) + P/4;
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

    dIvdu = crit_Iv;
    dIvdh = -2*crit_Iv;
    dIvdp = crit_Iv/p_p_dl;

    dmudu = dmudIv_fn(crit_Iv).*dIvdu;
    dmudp = dmudIv_fn(crit_Iv).*dIvdp;
    dmudh = dmudIv_fn(crit_Iv).*dIvdh;

    dpsidIv = phi_c/2/(1+sqrt(crit_Iv))^2/sqrt(crit_Iv);
    
    A_mat(1,1) = k; %+2*1i*P/beta_dl
    A_mat(1,2) = k;
    A_mat(1,3) = -2*P*1i/beta_dl;
    
    A_mat(2,1) = g_dl*cosd(theta)*k-P*g_dl*cosd(theta)*dmudh*1i;
    A_mat(2,2) = k - 1i*P*g_dl*cosd(theta)*dmudu-1i*k^2*nu_dl; 
    A_mat(2,3) = 1i*tand(theta)/(rho_dl-rho_f_dl) - 1i*P*g_dl*cosd(theta)*dmudp;
    
    A_mat(3,1) = -2*3*1i/alpha_dl*dpsidIv*dIvdh; %  - 1i*zeta*2/beta_dl-2/beta_dl*P*1i;
    A_mat(3,2) = (chi_dl*rho_dl-rho_f_dl)*g_dl*cosd(theta)*k - 2*3*1i/alpha_dl*dpsidIv*dIvdu;
    A_mat(3,3) = k -2*3*1i/alpha_dl*dpsidIv*dIvdp  + 1i*2/beta_dl*(P-zeta);
    A_mat(3,4) = -2*3*1i*crit_phi/alpha_dl;
    
    A_mat(4,3) = (P+rho_f_dl/rho_dl)*2/beta_dl*1i;
    A_mat(4,4) = k;
    
end


