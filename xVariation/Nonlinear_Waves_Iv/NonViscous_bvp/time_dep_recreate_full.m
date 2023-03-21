function [xi_merge,y_merge] = time_dep_recreate_full(h0_dim,theta,lambda_dim,tau0,d,alpha,rel_flux,pres_h)
% Converts the master wave stored in no_pe_no_vis_master.txt into a 
% waveform that maintains the flux of the conditions specified. Allows 
% change in theta, Froude number and yield stress.

% Can be run with no inputs and just with the parameter values set in the
% code

% Can be run with specify_param set to true and with params as a 3 entry
% list [Fr,theta,tau0]. In this case runs from master wave to specified 
% values.

% Can be run with specify_param and provide_init set to true. This case
% needs params as the target 3 value list, master_xi as the initial value
% xi, master_y as the initial y and master_params as the initial parameters

    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    
    [Fr_eq,~] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0_dim, tau0, 0);
    lambda = lambda_dim/h0_dim;
    u_eq_dim = Fr_eq*sqrt(g*cosd(theta)*h0_dim);
    
    if ~exist("rel_flux","var")
        rel_flux=1;
    end
    
    if ~exist("pres_h","var")
        pres_h=false;
    end
    
    [xi_final,y_final] = full_depo_wave_sinful(true,[Fr_eq,theta,lambda,false,tau0,alpha,d,pres_h,rel_flux]);
    if size(y_final,1) == 14
        stat_dist = y_final(14,1);
    else
        stat_dist = 0;
    end
    y_merge = horzcat(y_final(4:8,:),y_final(9:13,2:end));
    y_merge(2,size(xi_final)) = (y_final(5,1)+y_final(10,end))/2;
    y_merge = vertcat(y_final(1,1)*ones(1,size(y_merge,2)),y_merge);
    
    xi_h_crit = y_final(3,1);
    xi_merge = horzcat((stat_dist+xi_final*(xi_h_crit-stat_dist)),xi_h_crit+xi_final(2:end)*(lambda-xi_h_crit))*lambda_dim/lambda;
    y_merge = [u_eq_dim,u_eq_dim*h0_dim,h0_dim,u_eq_dim*h0_dim,h0_dim,rho_f*g*cosd(theta)*h0_dim]'.*y_merge;
end
    