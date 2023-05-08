function theta_out = get_crit_wave_theta(h, rho_p, rho_f, d, eta_f, alpha)
% Finds the critical theta for linear instability for given h and other
% parameters
%     crit_Iv = nan;
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    theta_min = atand(mu1_Iv*(rho-rho_f)/rho)+1e-1;
    theta_max = 30;
    theta_tol = 1e-3;

    while ((theta_max-theta_min)>theta_tol)
        theta_mid = (theta_max+theta_min)/2;
        crit_Iv = newt_solve_crit_Iv(theta_mid, rho_p, rho_f);
        crit_phi = phi_c./(1+sqrt(crit_Iv));
        p_p = (rho_p-rho_f)*g*phi_c*cosd(theta_mid)*h;

        crit_pb = rho_f*g*cosd(theta_mid)*h;
        crit_u = crit_Iv/eta_f/3*p_p*h;
        
        Fr = crit_u/sqrt(g*cosd(theta_mid)*h);

        v_scale = crit_u;
        p_scale = crit_pb;
        z_scale = h;
        t_scale = z_scale/v_scale;

        d_dl = d/z_scale;
        eta_f_dl = eta_f/(p_scale*t_scale);
        alpha_dl = alpha*p_scale;
        g_dl = g*t_scale/v_scale;

        p_p_dl = p_p/p_scale;
        rho_f_dl = rho_f*v_scale^2/p_scale;
        rho_p_dl = rho_p*v_scale^2/p_scale;
        rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);

        chi_dl = (rho_f_dl+3*rho_dl)/(4*rho_dl);
        P = (rho_dl-rho_f_dl)/rho_dl;
        zeta = 3/(2*alpha_dl) + P/4;
        A_dl = sqrt(crit_Iv);
        beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

        dmudu = dmudIv_fn(crit_Iv).*crit_Iv;
        dmudp = dmudIv_fn(crit_Iv).*-crit_Iv/p_p_dl;
        dmudh = dmudIv_fn(crit_Iv).*-2*crit_Iv;

        num_k = 20;
        A_mat = zeros(4);
        sigma_mat = zeros(4,num_k);

        A_mat(1,3) = -2*P*1i/beta_dl;
        A_mat(2,3) = 1i*tand(theta_mid)/(rho_dl-rho_f_dl) - 1i*P*g_dl*cosd(theta_mid)*dmudp;
        A_mat(3,1) =  3*A_dl*1i/alpha_dl/(1+A_dl)^2*phi_c; %- 1i*zeta*2/beta_dl-2/beta_dl*P*1i;
        A_mat(3,4) = -3*1i*crit_phi/alpha_dl;
        A_mat(4,1) = 1i * -(P-rho_f_dl/rho_dl).*2./beta_dl;
        A_mat(4,3) = (P+rho_f_dl/rho_dl)*2/beta_dl*1i;

        for i=1:num_k
            k = 0.25*2^(i/2);
            A_mat(1,1) = k+2*1i*P/beta_dl;
            A_mat(1,2) = k;
            A_mat(2,1) = g_dl*cosd(theta_mid)*k-P*g_dl*cosd(theta_mid)*dmudh*1i;
            A_mat(2,2) = k - 1i*P*g_dl*cosd(theta_mid)*dmudu; 
            A_mat(3,2) = (chi_dl*rho_dl-rho_f_dl)*g_dl*cosd(theta_mid)*k - 3/2*A_dl/(1+A_dl)^2*phi_c/alpha_dl*1i;
            A_mat(3,3) = k + 1i*2/beta_dl*(P-zeta) - 3/2*A_dl/(1+A_dl)^2*phi_c/alpha_dl*1i/P;
            A_mat(4,4) = k;

            A_eig = eig(A_mat);
            sigma_mat(:,i) = sort(imag(A_eig),'descend');
        end
        stable = (max(max(sigma_mat))<0);
        if stable
            theta_min = theta_mid;
        else
            theta_max = theta_mid;
        end
    end
    theta_out = (theta_max+theta_min)/2;
end