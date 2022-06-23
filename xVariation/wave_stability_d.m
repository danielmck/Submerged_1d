function wave_stability_d
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000; % kg/m^3
    rho_p = 2500; % kg/m^3
    theta = 10; % deg
    rho = phi_c*rho_p + (1-phi_c)*rho_f;

    alpha = 1e-4; % 1/Pa
    crit_Iv = newt_solve_crit_Iv(theta, rho_p, rho_f);

    if (crit_Iv>0)
        crit_phi = phi_c./(1+sqrt(crit_Iv));

        n_pts = 100;

        max_sig = zeros(n_pts);
        A_mat = zeros(4);

        Fr_list = linspace(0.005,0.5,n_pts);
        d_list = 2e-5.*2.^(linspace(0.1,10,n_pts));

        for j = 1:n_pts
            Fr = Fr_list(j);
            u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta);
            h = ((Fr*sqrt(g*cosd(theta)))./u_const)^(2/3); % layer height (m)
            p_p = (rho_p-rho_f)*g*phi_c*cosd(theta)*h;

            crit_pb = rho_f*g*cosd(theta)*h;
            crit_u = crit_Iv/eta_f/2*p_p*h;

            v_scale = crit_u;
            p_scale = crit_pb;
            z_scale = h;
            t_scale = z_scale/v_scale;

            p_p_dl = p_p/p_scale;
            eta_f_dl = eta_f/(p_scale*t_scale);
            alpha_dl = alpha*p_scale;
            g_dl = g*t_scale/v_scale;

            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale;
            rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);

            chi_dl = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            P = (rho_dl-rho_f_dl)/rho_dl;
            zeta = 3/(2*alpha_dl) + P/4;
            A_dl = sqrt(2*eta_f_dl/(g_dl*cosd(theta)*(rho_p_dl-rho_f_dl)*phi_c));

            dmudu = dmudIv_fn(crit_Iv).*crit_Iv;
            dmudh = dmudIv_fn(crit_Iv).*-2*crit_Iv;

            A_mat(2,3) = 1i*tand(theta)/(rho_dl-rho_f_dl);
            A_mat(3,4) = -3*1i*crit_phi/alpha_dl;

            for l = 1:n_pts  
                d_dl = d_list(l);

                beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

                A_mat(1,3) = -2*P*1i/beta_dl;
                A_mat(4,1) = 1i * -(P-rho_f_dl/rho_dl).*2./beta_dl;
                A_mat(4,3) = (P+rho_f_dl/rho_dl)*2/beta_dl*1i;
                A_mat(3,1) =  3*A_dl*1i/alpha_dl/(1+A_dl)^2*phi_c; %- 1i*zeta*2/beta_dl-2/beta_dl*P*1i;

                num_k = 16;
                sigma_mat = zeros(4,num_k);
                for i=1:num_k
                    k = 0.25*2^(i/2);
                    A_mat(1,1) = k; %+2*1i*P/beta_dl;
                    A_mat(1,2) = k;
                    A_mat(2,1) = g_dl*cosd(theta)*k-P*g_dl*cosd(theta)*dmudh*1i;
                    A_mat(2,2) = k - 1i*P*g_dl*cosd(theta)*dmudu; 
                    A_mat(3,2) = (chi_dl*rho_dl-rho_f_dl)*g_dl*cosd(theta)*k - 3/2*A_dl/(1+A_dl)^2*phi_c/alpha_dl*1i;
                    A_mat(3,3) = k + 1i*2/beta_dl*(P-zeta)- 3/2*A_dl/(1+A_dl)^2*phi_c/alpha_dl*1i/P;
                    A_mat(4,4) = k;

                    A_eig = eig(A_mat);
                    sigma_mat(:,i) = sort(imag(A_eig),'descend');  
                end
    %             plot(0.25.*2.^((1:num_k)/2),sigma_mat)
                max_sig(j,l) = max(max(sigma_mat));
            end  
        end
        stability = (max_sig>0);
        contourf(Fr_list,d_list,stability',2)
        SetPaperSize(10,10)
        colormap(winter)
        xlabel('Froude Number')
        ylabel('Nondimensional Particle Size')
        title("Stability Criteria for Normal Parameter Set")
        fig_name = 'StabCrit_Fr_d_Norm_10deg';
        PrintFig(fig_name)
        full_fig = strcat(fig_name,'.pdf');
        movefile(full_fig, 'Figures/StabilityPlots');
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end

    function dmudIv = dmudIv_fn(Iv)
        dmudIv = (mu2_Iv-mu1_Iv)*Iv_0./(Iv_0+abs(Iv)).^2 + 1 + 5/4*phi_c./sqrt(Iv);
    end
end