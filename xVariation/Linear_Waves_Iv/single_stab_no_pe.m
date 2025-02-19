function fn_out = single_stab_no_pe(Fr,crit_Iv,theta, rho_p, rho_f, eta_f, tau0)
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta);
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    h = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3); % layer height (m)
    p_p = (rho_p-rho_f)*g*phi_c*cosd(theta)*h;
    crit_pb = rho_f*g*cosd(theta)*h;
    crit_u = crit_Iv/eta_f/3*p_p*h;
    
    v_scale = crit_u;
    p_scale = crit_pb;
    z_scale = h;
    t_scale = z_scale/v_scale;
    
    eta_f_dl = eta_f/(p_scale*t_scale);
    g_dl = g*t_scale/v_scale;
    tau0_dl= tau0/p_scale;

    p_p_dl = p_p/p_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    rho0_dl = rho*v_scale^2/p_scale;
    P0 = (rho0_dl-rho_f_dl)/rho0_dl;
    dIvdu = crit_Iv;
    dIvdh = -2*crit_Iv;
    
    dmudIv = dmudIv_fn(crit_Iv);
    dmudu = dmudIv.*dIvdu;
    dmudh = dmudIv.*dIvdh;
    
    if (tau0_dl> 0)
        p= -rho0_dl*P0*dmudu/tau0_dl;
        q= -rho0_dl*P0*dmudh/tau0_dl;
        delta = -(4*p^3+27*q^2);
        delta = delta;
    end
    
    k = 10;
    A_mat = zeros(2);
    A_mat(1,1) = k;
    A_mat(1,2) = k;
    A_mat(2,1) = k/Fr^2-1i*P0/Fr^2*dmudh+1i*tau0_dl/rho0_dl;
    A_mat(2,2) = k-1i*P0/Fr^2*dmudu;
    
    [eig_vec,A_eig] = eigs(A_mat);
    A_eig = diag(A_eig);
    top_eig = max(imag(A_eig));
%         [top_eigs,~]=maxk(imag(A_eig),2);
%         top_eig = top_eigs(2);
    sigma_mat = sort(imag(A_eig),'descend');
    num_unstab = sum(sigma_mat>0);
    k_unstab=1;
    fn_out = [num_unstab k_unstab];
end