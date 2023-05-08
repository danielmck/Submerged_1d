function crit_wave_theta_val
% Finds the critical value of theta for a given height of flow, along with
% other parameters
    rho_p = 2500;
    
%     rho_f = 1000;
%     eta_f = 0.0010016; % Pa s
    
%     rho_f = 1; % kg/m^3
%     eta_f = 1.18e-5; % Pa s

    rho_list = [1000 1];
    eta_list = [0.0010016 1.18e-5];
%     alpha = 1e-4; % 1/Pa
    phi_c = 0.585;
%     d_dl = 2e-4;
    d = 1e-3;
    
%     alpha_list = [1e-6 5e-6 1e-5 5e-5 1e-4];
%     
    
    h_list = [1e-2 1e-1 1];
    n_pts1 = size(h_list,2);
    
    n_pts2 = 100;
    alpha_start = -6;
    alpha_stop = log10(5e-4);
    alpha_list = logspace(alpha_start,alpha_stop,n_pts2);
%     theta_start = [10 20];
%     theta_stop = [30 35];
    phase_name = ["water" "air"];
%     theta_list = (linspace(10,14,n_pts2));  
    
    crit_theta = zeros(n_pts1,n_pts2);
    for a = 1:2
        rho_f = rho_list(a);
        eta_f = eta_list(a);
            for i = 1:n_pts1
                h = h_list(i);
                for j = 1:n_pts2
                    alpha = alpha_list(j);
                    crit_theta(i,j) = get_crit_wave_theta(h, rho_p, rho_f, d, eta_f, alpha);
                end
            end
            save("Results/crit_theta_h_alpha_"+phase_name(a)+".txt", 'crit_theta','-ascii')
    end
    
end