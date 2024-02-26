function fn_out = Fr_stab_noneq(h0,u0,phi0,pb0,theta, rho_p, rho_f, d, eta_f, alpha,tau0)
% Finds the critical Froude number for the specified condition
    short = 2*4.3034e-04;
    long = 1.8;
    short_k = 2*pi/short;
    long_k = 2*pi/long;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    
    [Fr,crit_Iv] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0, tau0, false,true);
    crit_phi = phi_c./(1+sqrt(crit_Iv));
    crit_pp = (rho_p-rho_f)*g*crit_phi*cosd(theta)*h0;
    
    crit_pb = rho_f*g*cosd(theta)*h0;
    crit_u = Fr*sqrt(g*cosd(theta)*h0);
    
    rho0 = rho_p*phi0+rho_f*(1-phi0);
    pp0 = rho0*g*cosd(theta)*h0-pb0;
    Iv0 = 3*u0*eta_f/h0/pp0;
    
    v_scale = crit_u;
    p_scale = crit_pb;
    z_scale = h0;
    t_scale = z_scale/v_scale;
    
    short_k_dl = short_k*z_scale;
    long_k_dl = long_k*z_scale;

    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale;

    crit_pp_dl = crit_pp/p_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    rho0_dl = rho_p_dl*phi0+rho_f_dl*(1-phi0);
    d_dl = d/z_scale;
    tau_dl=tau0/p_scale;
    
    u_ratio = u0/crit_u;
    phi_ratio = phi0/crit_phi;
    pb_ratio = pb0/crit_pb;
    pp0_dl = pp0/crit_pb;

    chi0_dl = (rho_f_dl+3*rho0_dl)/(4*rho0_dl);
    P0 = (rho0_dl-rho_f_dl)/rho0_dl;
    zeta = 3/(2*alpha_dl) + P0/4;
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

    dIvdu = Iv0/u_ratio;
    dIvdh = -Iv0*(2*rho0_dl-rho_f_dl)/(rho0_dl-rho_f_dl);
    dIvdp = Iv0/pp0_dl;
    dIvdphi = -Iv0/phi0;

    dmudu = dmudIv_fn(Iv0).*dIvdu;
    dmudp = dmudIv_fn(Iv0).*dIvdp;
    dmudh = dmudIv_fn(Iv0).*dIvdh;
    dmudphi = dmudIv_fn(Iv0).*dIvdphi;
            
    dpsidIv = phi_c/2/(1+sqrt(Iv0))^2/sqrt(Iv0);
    
    D = -2/beta_dl*(pb_ratio-1);
    tanpsi = phi0-phi_c/(1+sqrt(Iv0));
    dhdt = P0*D;
    dhudt = sind(theta)-u_ratio*P0*D+(mu_Iv_fn(Iv0)*pp0_dl)/rho0_dl;
    dudt = dhudt-u_ratio*dhdt;
    dhphidt = phi0*rho_f_dl/rho0_dl*D;
    dphidt = (dhphidt-phi0*dhdt);
    dpbdt = zeta*D-9/2/alpha_dl*u_ratio*tanpsi;

    num_k = 50;
    A_mat = zeros(4);
    sigma_mat = zeros(4,num_k);
    phase_ang = zeros(4,num_k);

    k_unstab = NaN;
    max_eig = NaN;
    
    k_val = logspace(log10(long_k_dl),log10(short_k_dl),num_k);
    for i = 1:num_k
        k = k_val(i);
        A_mat(1,1) = u_ratio*k + 2*1i*P0/beta_dl*pb_ratio;
        A_mat(1,2) = k;
        A_mat(1,3) = -2*P0*1i/beta_dl;
        A_mat(1,4) = -1i*D*rho_f_dl*(rho_p_dl-rho_f_dl)/rho0_dl^2;
        
        A_mat(2,1) = g_dl*cosd(theta)*k+1i*g*sin(theta)-mu_Iv_fn(Iv0)*g_dl*cosd(theta)*1i-1i*pp0_dl/rho0_dl*dmudh+u_ratio*P0*2/beta_dl*pb_ratio*1i; %-dudt*1i;
        A_mat(2,2) = u_ratio*k + P0*D*1i -1i*pp0_dl/rho0_dl*dmudu+P0*D*i; %-dhdt*1i;
        A_mat(2,3) = -mu_Iv_fn(Iv0)/rho0_dl*1i-1i*pp0_dl/rho0_dl*dmudp-u_ratio*1i*D*rho_f_dl*(rho_p_dl-rho_f_dl)/rho0_dl^2;
        A_mat(2,4) = -1i*tau_dl/rho0_dl*(rho_p_dl-rho_f_dl)-1i*mu_Iv_fn(Iv0)*pb_ratio*(rho_p_dl-rho_f_dl)/rho0_dl-1i*pp0_dl/rho0_dl*dmudphi - u_ratio*2*P0*1i/beta_dl;
        
        A_mat(3,1) = 1i*zeta*2/beta_dl*pb_ratio+1i*D*(-3/2/alpha_dl)+9/2*1i/alpha_dl*u_ratio*tanpsi-9/2*1i/alpha_dl*u_ratio*dpsidIv*dIvdh; %  -2/beta_dl*P0*1i;
        A_mat(3,2) = chi0_dl*rho0_dl*g_dl*cosd(theta)*k - 9/2*1i/alpha_dl*tanpsi- 9/2*1i/alpha_dl*u_ratio*dpsidIv*dIvdu;
        A_mat(3,3) = u_ratio*k + 1i*(-2/beta_dl)*zeta-9/2*1i/alpha_dl*u_ratio*dpsidIv*dIvdp;  
        A_mat(3,4) = 1i*D*(g_dl*cosd(theta)*rho_f_dl^2*(rho_p_dl-rho_f_dl)/4/rho0_dl^2)-9/2*1i/alpha_dl*u_ratio*(1+dpsidIv*dIvdphi);
        
        A_mat(4,1) = -phi0*2/beta_dl*pb_ratio*1i; %-dphidt*1i;
        A_mat(4,3) = phi0*2/beta_dl*1i;
        A_mat(4,4) = u_ratio*k-1i*D/rho0_dl*(phi0*rho_f_dl*(rho_p_dl-rho_f_dl)/rho0_dl-rho_f_dl+phi0*rho_f_dl/rho0_dl*(rho_p_dl-rho_f_dl)); %-1i*dhdt;

        [eig_vec,A_eig] = eigs(A_mat);
        A_eig = diag(A_eig);
%         top_eig = max(imag(A_eig));
        [top_eigs,~]=maxk(imag(A_eig),2);
        top_eig = top_eigs(2);
        sigma_mat(:,i) = sort(imag(A_eig),'descend');
        
        t_val = 0.001/t_scale;
%         hold on
%         for l=1:4
%             exp_part = eig_vec(l,:)*exp(-i*t_val*A_eig(l));
%             plot(exp_part)
%         end
        
        if ((top_eig > max_eig) || (i==1))
            max_eig = top_eig;
            k_unstab = k;
        end
%                 phase_ang = phase(A_vec(ind,:));
    end
    max_k = max(k_val(max(sigma_mat,[],2)>0));
    num_unstab = sum(max(sigma_mat,[],2)>0);
    fn_out = [num_unstab,1/max_k*z_scale];
%     f = figure();
    SetPaperSize(8,8)
%     set(f, 'PaperUnits', 'centimeters');
%     set(f, 'PaperSize', [10 10]);
%     semilogx(k_val,sigma_mat)
    xvals = linspace(0,short/z_scale,100);
    comp = eig_vec(:,4).*exp(1i*k*xvals);
    plot(xvals,comp)
%     legend('$h$','$u$', '$p_b$', '$\phi$','Location', 'best')
%     ylim([-100,100])
%     ax=gca;
%     ax.YAxis.Exponent = 0;
    xlabel('Wavenumber')
    ylabel('Eigenvector real part')
    title('Non steady state: Unstable mode')
%     title('$\theta = 18$, $d = 10^{-4}$, $\alpha = 10^{-4}$, $Fr = 5$, Phase = Air')
    exp_graph(gca,"Noneq_eigvec_4.pdf")
%     
end