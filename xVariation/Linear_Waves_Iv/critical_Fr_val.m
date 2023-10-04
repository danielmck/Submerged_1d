function critical_Fr_val
% Finds the critical Froude number for instability of linear waves using 
% get_critical_Fr
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
%     d = 2e-4;
    
%     alpha_list = [1e-6 5e-6 1e-5 5e-5 1e-4];
%     n_pts1 = size(alpha_list,2);

%     var_list = {[1e-5 1e-4 1e-3],[1e-3 5e-3 1e-2 5e-2 1e-1]};
    var_list = {[1e-5 1e-4 1e-3],[1e-6 5e-6 1e-5 5e-5 1e-4]};
    var_names = ["d" "alpha"];

    n_pts2 = 100;
    theta_start = [8.6 17.8];
    theta_stop = [30 32];
%     theta_start = [10 20];
%     theta_stop = [30 35];
    phase_name = ["water" "air"];
%     theta_list = (linspace(10,14,n_pts2));
%     
    
    for a = 2
        theta_list = (linspace(theta_start(a),theta_stop(a),n_pts2));
        rho_f = rho_list(a);
        eta_f = eta_list(a);
        for b = 1:2
%             alpha = 1e-2; % 1/Pa
            alpha = 1e-5;
            d = 1e-4;
            n_pts1 = size(var_list{b},2);
            crit_Fr = zeros(n_pts1,n_pts2);
            for i = 1:n_pts1
                if (b == 1)
                    d = var_list{b}(i);
                else
                    alpha = var_list{b}(i);
                end
                for j = 1:n_pts2
                    theta = theta_list(j);
%                     rho = phi_c*rho_p + (1-phi_c)*rho_f;
        %             crit_mu = rho/(rho-rho_f)*tand(theta);
                    crit_Fr(i,j) = get_critical_Fr(theta, rho_p, rho_f, d, eta_f, alpha);
                end
            end
            save("Results/crit_Fr_"+var_names(b)+"_theta_"+phase_name(a)+".txt", 'crit_Fr','-ascii')
        end
    end
    
end