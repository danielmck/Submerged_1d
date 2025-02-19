function fn_out = two_eqn_stab(Fr,theta, rho_p, rho_f, eta_f, tau0)
% Finds the critical Froude number for the specified condition
    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    
    [h0,crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,0,false);
    short = h0/1000;
    long = 50000*h0;
    short_k = 2*pi/short;
    long_k = 2*pi/long;
    crit_phi = phi_c; %./(1+sqrt(crit_Iv));
    u_const = crit_Iv/eta_f/3*(rho_p-rho_f)*g*crit_phi*cosd(theta);
    rho_0 = crit_phi*rho_p + (1-crit_phi)*rho_f;
    
    p_p = (rho_p-rho_f)*g*crit_phi*cosd(theta)*h0;
    
    crit_phi = phi_c;
    crit_pb = rho_f*g*cosd(theta)*h0;
    crit_u = crit_Iv/eta_f/3*p_p*h0;
    
    v_scale = crit_u;
    p_scale = crit_pb;
    z_scale = h0;
    t_scale = z_scale/v_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    g_dl = g*t_scale/v_scale;

    p_p_dl = p_p/p_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    rho0_dl = rho_p_dl*crit_phi+rho_f_dl*(1-crit_phi);
    tau0_dl = tau0/p_scale;

    P0 = (rho0_dl-rho_f_dl)/rho0_dl;

    dIvdu = crit_Iv;
    dIvdh = -2*crit_Iv;

    dmudu = dmudIv_fn(crit_Iv).*dIvdu;
    dmudh = dmudIv_fn(crit_Iv).*dIvdh;

    num_k = 100;
    A_mat = zeros(2);
    sigma_mat = zeros(2,num_k);
    phase_ang = zeros(2,num_k);

    k_unstab = NaN;
    max_eig = NaN;
    
    k_val = logspace(-5,5,num_k);
    % short_k_dl = short_k*z_scale;
    % long_k_dl = long_k*z_scale;
    % k_val = logspace(log10(long_k_dl),log10(short_k_dl),num_k);
    for i = 1:num_k
        k = 10;%k_val(i);
        A_mat(1,1) = k;
        A_mat(1,2) = k;
        A_mat(2,1) = g_dl*cosd(theta)*k-1i*P0*g_dl*cosd(theta)*dmudh+1i*tau0_dl/rho0_dl;
        A_mat(2,2) = k - 1i*P0*g_dl*cosd(theta)*dmudu; 

        [eig_vec,A_eig] = eigs(A_mat);
        A_eig = diag(A_eig);
%         top_eig = max(imag(A_eig));
        [top_eigs,~]=maxk(imag(A_eig),2);
        top_eig = top_eigs(2);
        sigma_mat(:,i) = sort(imag(A_eig),'descend');
        
        if ((top_eig > max_eig) || (i==1))
            max_eig = top_eig;
            k_unstab = k;
        end
%                 phase_ang = phase(A_vec(ind,:));
    end
    num_unstab = sum(max(sigma_mat,[],2)>0);
    fn_out = [num_unstab k_unstab];
    f = figure;
    SetPaperSize(10,10)
%     set(f, 'PaperUnits', 'centimeters');
%     set(f, 'PaperSize', [10 10]);
%     semilogx(k_val,sigma_mat)
    xvals = linspace(0,short/z_scale,100);
    comp = eig_vec(:,2).*exp(1i*k*xvals);
    semilogx(k_val,sigma_mat)
%     legend('$h$','$u$', '$p_b$', '$\phi$','Location', 'best')
    xlabel('Wavenumber $k$')
    ylabel('Linear Growth Rate $\Im(\omega)$')
    % title('Two eqn model: Stable mode')
    title('$\theta = 10^{\circ}$, $Fr = 1$, $\tau_0=10$Pa, $h_0 = 0.04$m')
    % exp_graph(f,"TwoEqn_Eigvec.pdf")
%     
end