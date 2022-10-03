function Fr_out = get_critical_Fr_I(theta, rho_p, rho_f, d, eta_f, alpha)
%     crit_Iv = nan
%     mu1_I=0.342; 
%     mu2_I=0.557;  
%     I_0 = 0.069;
    
    mu1_I=tand(20); 
    mu2_I=tand(33);  
    I_0 = 0.3;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    delta_phi = 0.2;
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;

    crit_I = newt_solve_crit_I(theta, rho_p, rho_f);

    if (crit_I > 0)
        Fr_min = 1e-4;
        Fr_max = 10;
        Fr_tol = 1e-4;

        while ((Fr_max-Fr_min)>Fr_tol)
            Fr_mid = (Fr_max+Fr_min)/2;
            
%             crit_u = Fr_mid*sqrt(g*cosd(theta)*h);
%             d_dl = crit_I/5/crit_u/sqrt(rho_p)*sqrt((rho_p-rho_f)*g*phi_c*cosd(theta)*h);
            
            u_const = crit_I*sqrt((rho_p-rho_f)*g*phi_c*cosd(theta))/sqrt(rho_p)/5/d;
            h = ((Fr_mid*sqrt(g*cosd(theta)))./u_const); % layer height (m)
            crit_u = u_const*h^(3/2);
            crit_phi = phi_c-delta_phi*crit_I;
            p_p = (rho_p-rho_f)*g*phi_c*cosd(theta)*h;

            crit_pb = rho_f*g*cosd(theta)*h;


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

            chi_dl = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            P = (rho_dl-rho_f_dl)/rho_dl;
            zeta = 3/(2*alpha_dl) + P/4;
            
            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            num_k = 16;
            A_mat = zeros(4);
            sigma_mat = zeros(4,num_k);
            
            dIdu = crit_I;
            dIdh = -3/2*crit_I;
            dIdp = 1/2*crit_I/(((rho_dl - rho_f_dl)*g_dl*cosd(theta)));
            
            dmudu = dmudI_fn(crit_I).*dIdu;
            dmudh = dmudI_fn(crit_I).*dIdh;
            dmudp = dmudI_fn(crit_I).*dIdp;

            A_mat(1,3) = -2*P*1i/beta_dl;
            A_mat(2,3) = 1i*tand(theta)/(rho_dl-rho_f_dl) - 1i*P*g_dl*cosd(theta)*dmudp;
            A_mat(3,1) =  -3*1i/alpha_dl*delta_phi*dIdh; %- 1i*zeta*2/beta_dl-2/beta_dl*P*1i;
            A_mat(3,4) = -3*1i*crit_phi/alpha_dl;
            A_mat(4,1) = 1i * -(P-rho_f_dl/rho_dl).*2./beta_dl;
            A_mat(4,3) = (P+rho_f_dl/rho_dl)*2/beta_dl*1i;
            
            for i=1:num_k
                k = 0.25*2^(i/2);
                A_mat(1,1) = k; %+2*1i*P/beta_dl;
                A_mat(1,2) = k;
                A_mat(2,1) = g_dl*cosd(theta)*k-P*g_dl*cosd(theta)*dmudh*1i;
                A_mat(2,2) = k - 1i*P*g_dl*cosd(theta)*dmudu; 
                A_mat(3,2) = (chi_dl*rho_dl-rho_f_dl)*g_dl*cosd(theta)*k - 3*1i/alpha_dl*delta_phi*dIdu;
                A_mat(3,3) = k + 1i*2/beta_dl*(P-zeta) - 3/2*1i/alpha_dl*delta_phi*dIdp;
                A_mat(4,4) = k;

                A_eig = eig(A_mat);
                sigma_mat(:,i) = sort(imag(A_eig),'descend');
            end
            stable = (max(max(sigma_mat))<0);
            if stable
                Fr_min = Fr_mid;
            else
                Fr_max = Fr_mid;
            end
        end
        Fr_out = (Fr_max+Fr_min)/2;
    else
        Fr_out = -1;
    end
    
    function mu_val = mu_I_fn(I)
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end
        
    function dmudI = dmudI_fn(I)
        dmudI = (mu2_I-mu1_I)*I_0./(I_0+abs(I)).^2;
    end
end