function Fr_out = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha)
%     crit_Iv = nan;
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    g=9.81; % m/s^2
    rho = phi_c*rho_p + (1-phi_c)*rho_f;
    
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);

    if (crit_Iv > 0)
        crit_phi = phi_c./(1+sqrt(crit_Iv));
        u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
        
        Fr_min = 1e-4;
        Fr_max = 20.0;
        Fr_tol = 1e-5;

        while ((Fr_max-Fr_min)>Fr_tol)
            Fr_mid = (Fr_max+Fr_min)/2;
            h = ((Fr_mid*sqrt(g*cosd(theta)))./u_const)^(2/3); % layer height (m)
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

            p_p_dl = p_p/p_scale;
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale;
            rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
            d_dl = d/z_scale;

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
            
            num_k = 20;
            A_mat = zeros(4);
            sigma_mat = zeros(4,num_k);
            phase_ang = zeros(4,num_k);
            
            A_mat(1,3) = -2*P*1i/beta_dl;
            A_mat(2,3) = 1i*tand(theta)/(rho_dl-rho_f_dl) - 1i*P*g_dl*cosd(theta)*dmudp;
            A_mat(3,1) = -2*3*1i/alpha_dl*dpsidIv*dIvdh; %  - 1i*zeta*2/beta_dl-2/beta_dl*P*1i;
            A_mat(3,4) = -2*3*1i*crit_phi/alpha_dl;
            A_mat(4,1) = 1i * -(P-rho_f_dl/rho_dl).*2./beta_dl;
            A_mat(4,3) = (P+rho_f_dl/rho_dl)*2/beta_dl*1i;

            k_val = 0.25.*2.^((1:num_k)/2);
            for i = 1:num_k
                k = k_val(i);
                A_mat(1,1) = k; %+2*1i*P/beta_dl
                A_mat(1,2) = k;
                A_mat(2,1) = g_dl*cosd(theta)*k-P*g_dl*cosd(theta)*dmudh*1i;
                A_mat(2,2) = k - 1i*P*g_dl*cosd(theta)*dmudu; 
                A_mat(3,2) = (chi_dl*rho_dl-rho_f_dl)*g_dl*cosd(theta)*k - 2*3*1i/alpha_dl*dpsidIv*dIvdu;
                A_mat(3,3) = k -2*3*1i/alpha_dl*dpsidIv*dIvdp  + 1i*2/beta_dl*(P-zeta);
                A_mat(4,4) = k;

                A_eig = eigs(A_mat);
%                 [eig ind] = max(imag(diag(A_eig)));
                sigma_mat(:,i) = sort(imag(A_eig),'descend');
%                 phase_ang = phase(A_vec(ind,:));
            end
            stable = (max(max(sigma_mat))<0);
            if stable
                Fr_min = Fr_mid;
            else
                Fr_max = Fr_mid;
            end
        end
        Fr_out = (Fr_max+Fr_min)/2;
        h_out = ((Fr_out*sqrt(g*cosd(theta)))./u_const)^(2/3);
    else
        Fr_out = -1;
        h_out = -1;
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
    
    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end