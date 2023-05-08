function comp_A_mat_I
% Compares the matrices of growth of the perturbations to the steady state
% flow across a number of different cases. The wavenumber of interest also
% needs to be specified.
    mu1_I=tand(20); 
    mu2_I=tand(33);  
    I_0 = 0.3;

    reg_param = 1*10^7;

    rho_p = 2500;
    
    rho_f = 1;
    eta_f = 1.18e-5;
    
%     rho_f = 1000;
%     eta_f = 0.0010016;
    
    d = 1e-4;
    alpha = 1e-4;
    
    phi_c=0.585; % Volume fraction
    delta_phi = 0.2;
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    theta = [30,31.5];
    n_cases = size(theta,2);
    
    A_mat = cell(2,1);
    
    k = 100;
    Fr = 0.2;
    for i = 1:n_cases
        crit_I(i) = newt_solve_crit_I(theta(i), rho_p, rho_f);
        u_const(i) = crit_I(i)*sqrt((rho_p-rho_f)*g*phi_c*cosd(theta(i)))/sqrt(rho_p)/5/d;
        h(i) = (Fr*sqrt(g*cosd(theta(i))))./u_const(i); % layer height (m)
        crit_u(i) = u_const(i)*h(i)^(3/2);
        crit_phi(i) = phi_c-delta_phi*crit_I(i);

        p_p(i) = (rho_p-rho_f)*g*phi_c*cosd(theta(i))*h(i);

        crit_pb(i) = rho_f*g*cosd(theta(i))*h(i);

        v_scale(i) = crit_u(i);
        p_scale(i) = crit_pb(i);
        z_scale(i) = h(i);
        t_scale(i) = z_scale(i)/v_scale(i);

        eta_f_dl(i) = eta_f/(p_scale(i)*t_scale(i));
        alpha_dl(i) = alpha*p_scale(i);
        g_dl(i) = g*t_scale(i)/v_scale(i);
        d_dl(i) = d/z_scale(i);

        rho_f_dl(i) = rho_f*v_scale(i)^2/p_scale(i);
        rho_p_dl(i) = rho_p*v_scale(i)^2/p_scale(i);
        rho_dl(i) = rho_p_dl(i)*phi_c+rho_f_dl(i)*(1-phi_c);
        p_p_dl(i) = p_p(i)/p_scale(i);

        chi_dl(i) = (rho_f_dl(i)+3*rho_dl(i))/(4*rho_dl(i));
        P(i) = (rho_dl(i)-rho_f_dl(i))/rho_dl(i);
        zeta(i) = 3/(2*alpha_dl(i)) + P(i)/4;
        beta_dl(i) = 150*phi_c.^2.*eta_f_dl(i)./((1-phi_c).^3.*d_dl(i)^2);
        k = 100;

        dIdu(i) = crit_I(i);
        dIdh(i) = -3/2*crit_I(i);
        dIdp(i) = 1/2*crit_I(i)/(((rho_dl(i) - rho_f_dl(i))*g_dl(i)*cosd(theta(i))));

        dmudu(i) = dmudI_fn(crit_I(i)).*dIdu(i);
        dmudh(i) = dmudI_fn(crit_I(i)).*dIdh(i);
        dmudp(i) = dmudI_fn(crit_I(i)).*dIdp(i);
        
        A_mat{i,1}(1,1) = k; %+2*1i*P/beta_dl(i);
        A_mat{i,1}(1,2) = k;
        A_mat{i,1}(1,3) = -2*P(i)*1i/beta_dl(i);
        
        A_mat{i,1}(2,1) = g_dl(i)*cosd(theta(i))*k-P(i)*g_dl(i)*cosd(theta(i))*dmudh(i)*1i;
        A_mat{i,1}(2,2) = k - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudu(i);
        A_mat{i,1}(2,3) = 1i*tand(theta(i))/(rho_dl(i)-rho_f_dl(i)) - 1i*P(i)*g_dl(i)*cosd(theta(i))*dmudp(i);
        
        A_mat{i,1}(3,1) =  -3*1i/alpha_dl(i)*delta_phi*dIdh(i); %- 1i*zeta*2/beta_dl(i)-2/beta_dl(i)*P*1i;
        A_mat{i,1}(3,2) = (chi_dl(i)*rho_dl(i)-rho_f_dl(i))*g_dl(i)*cosd(theta(i))*k - 3*1i/alpha_dl(i)*delta_phi*dIdu(i);
        A_mat{i,1}(3,3) = k + 1i*2/beta_dl(i)*(P(i)-zeta(i)) - 3/2*1i/alpha_dl(i)*delta_phi*dIdp(i);
        A_mat{i,1}(3,4) = -3*1i*crit_phi(i)/alpha_dl(i);
        
        A_mat{i,1}(4,1) = 1i * -(P(i)-rho_f_dl(i)/rho_dl(i)).*2./beta_dl(i);
        A_mat{i,1}(4,3) = (P(i)+rho_f_dl(i)/rho_dl(i))*2/beta_dl(i)*1i;
        A_mat{i,1}(4,4) = k;

        A_eig(i,:) = eig(A_mat{i,1});
    end
    
    function mu_val = mu_I_fn(I)
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end

    function dmudI = dmudI_fn(I)
        dmudI = (mu2_I-mu1_I)*I_0./(I_0+abs(I)).^2;
    end
end