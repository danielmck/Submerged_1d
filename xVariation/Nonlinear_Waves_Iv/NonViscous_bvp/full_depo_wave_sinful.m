function [xi_final,y_final] = full_depo_wave_sinful(specify_param,params,provide_init,master_xi,master_y,master_params)
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
    chi = (rho_f+3*rho)/(4*rho);
    
    if ~exist("specify_param","var")
        specify_param = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    if ~provide_init
        Res_dir = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/NonViscous_bvp/";
        master_name = "no_vis_full_master.txt";
        master_file = load(Res_dir+"Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable(Res_dir+"Results/wave_record.csv");

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table);
        master_Fr = record.Fr(in_table);
        master_tau0 = record.tau0(in_table);
        master_lambda = master_y(2,1);
        master_h_crit = master_y(3,1);
        master_d = record.d(in_table);
        master_alpha = record.alpha(in_table);
        master_pres_h=strcmp(wave_type(max(end-7,1):end),"_pres_h");
        master_rf = master_y(11,end);
        if (master_y(3,1)-master_y(2,1)/master_y(1,1)<1e-2)
            master_h_min = -1;
        else
            master_h_min = master_y(5,1);
        end
        if (size(master_y,1) == 14)
            master_stat_len = master_y(14,1);
        else
            master_stat_len = 0;
        end
    else
        master_cell = num2cell(master_params);
        [master_Fr,master_theta,master_alpha,master_d,master_tau0,master_pres_h] = master_cell{:};
    end
    
    if ~specify_param
        Fr_eq = 0.8; 
        theta = 12;
        tau0 = 0;
        
        d = 1e-4;
        alpha=1e-5;
        
        set_h_min = false;
        h_min = 0.8;
        lambda_final = 12;
        
        pres_h = true;
        rel_flux = 1;

        % h_alt specifies the distance from the static minimum that the
        % wave must begin at.
        filename = "full_no_vis_pres_h.txt";      
    else
        param_cell = num2cell(params);
        [Fr_eq,theta,wlen_val,set_h_min,tau0,alpha,d,pres_h,rel_flux] = param_cell{:};
        if set_h_min
            h_min = wlen_val;
            lambda_final = 0;
        else
            lambda_final = wlen_val;
            h_min = 0;
        end
    end
    
    if (pres_h ~= master_pres_h)
        m_init_l = zeros(size(master_xi));
        m_init_r = zeros(size(master_xi));
        m_val = pres_h*master_stat_len/master_lambda*master_y(4,1)/master_y(1,1);
        u_flux_l = (master_y(1,1) - master_y(4,:)./master_y(5,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init_l(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + (master_h_crit-master_stat_len)/master_lambda* (master_y(5,j)*u_flux_l(j)+master_y(5,j+1)*u_flux_l(j+1))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        u_flux_r = (master_y(1,1) - master_y(9,:)./master_y(10,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init_r(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + (master_lambda-master_h_crit)/master_lambda* (master_y(10,j)*u_flux_r(j)+master_y(10,j+1)*u_flux_r(j+1))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        master_y = vertcat(master_y(1:5,:),m_init_l,master_y(7:10,:),m_init_r,master_y(12:end,:));
        master_rf = m_val;
    end

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    d_ratio = max(d/master_d,master_d/d);
    alpha_ratio = max(alpha/master_alpha,master_alpha/alpha);
    
    hm_ratio = min(max(h_min/master_h_min,master_h_min/h_min),1e3);
    lambda_ratio = min(max(lambda_final/master_lambda,master_lambda/lambda_final),1e3);
    rf_ratio = abs(rel_flux - master_rf);
    ratio_sum = Fr_ratio-1+20*(theta_ratio-1)+tau0_ratio+hm_ratio*set_h_min+(1-set_h_min)*lambda_ratio+rf_ratio+(d_ratio-1)+(alpha_ratio-1);
    if ratio_sum > 1e-6
        n_step = max(min(ceil(40*ratio_sum),800),90);
    else
        n_step=2;
    end
    % If the parameters are the same as the initial values run 2 steps, if
    % not, run more
    Fr_list = linspace(master_Fr,Fr_eq,n_step);
    theta_list = linspace(master_theta,theta,n_step);
    tau0_list = linspace(master_tau0,tau0,n_step);
    d_list = linspace(master_d,d,n_step);
    alpha_list = linspace(master_alpha,alpha,n_step);
    if set_h_min
        if master_h_min>0 && h_min<0
            wlen_list = linspace(master_h_min,0,n_step);
        else
            wlen_list = linspace(max(master_h_min,master_y(4,1)/master_y(1,1)),h_min,n_step);
        end
    else
        wlen_list = linspace(master_lambda,lambda_final,n_step);
    end
    rf_list = linspace(master_rf, rel_flux, n_step);
    
    h_b_pert = 2e-3;

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, tau0_list, wlen_list, alpha_list, d_list, rf_list, master_xi, master_y);
    if ~specify_param
        out_vec = vertcat(xi_final,y_final);
        lambda = y_final(3,1);
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/wave_record.csv",filename,{"full_no_vis","water",Fr_eq,theta,h_alt,0,0,tau0})
    end
    
    function [xi_in, y_in] = run_bvp_step(Fr_vals, theta_vals, tau0_vals, wlen_vals, alpha_vals, d_vals, rf_vals, xi_in, y_in, tol,counter)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-3;
        end
        if counter > 10
%             out_vec = vertcat(xi_in,y_in);
%             fail_name = filename;
%             save("Results/"+fail_name,"out_vec","-ascii")
%             write_record("Results/wave_record.csv",fail_name,{"no_pe","water",Fr_vals(1),theta_vals(1),lambda_vals(1),nu_vals(1),0,0,tau0_vals(1)})
            error("Max iteration depth reached, non convergence")
        end
        h_pert_max = 1e-3;
        
        static_part = (size(y_in,1) == 14);
        nstep = size(Fr_vals,2);
        for i = 2:nstep
            if (y_in(5,1)<y_in(4,1)/y_in(1,1)+1.5*h_b_pert) && (lambda_vals(i)>=lambda_vals(i-1)) && ~static_part
                static_part = true;
                y_in = vertcat(y_in,zeros(size(xi_in)));
            elseif (static_part)
                if (y_in(end,1) < 1e-2) && (lambda_vals(i)<=lambda_vals(i-1))
                    static_part = false;
                    y_in = y_in(1:end-1,:);
                end
            end
            % Solves the 5 value bvp for u_w, Q1, h, n and m.
            h_pert = (y_in(10,1)-y_in(5,end))/2;
            solOut =  run_single_step(theta_vals(i),Fr_vals(i),tau0_vals(i),wlen_vals(i), alpha_vals(i), d_vals(i), ...
                rf_vals(i),h_pert,xi_in,y_in);
%             if solOut.stats.maxres > tol
%                 pert_jump = 2e-4;
%                 h_crit = (y_in(5,end)+y_in(10,1))/2;
%                 while h_pert < h_pert_max*0.99
%                         y_in(5,end) = h_crit-h_pert-pert_jump;
%                         y_in(10,1) = h_crit+h_pert+pert_jump;
%                         solPert = run_single_step(theta_vals(i-1),Fr_vals(i-1),tau0_vals(i-1),wlen_vals(i-1), alpha_vals(i-1), d_vals(i-1), ...
%                             rf_vals(i-1),h_pert+pert_jump,xi_in,y_in);
%                         if solPert.stats.maxres < tol
%                             y_in = solPert.y;
%                             xi_in = solPert.x;
%                             h_pert = h_pert + pert_jump;
%                         else
%                             pert_jump = pert_jump/2;
%                         end
%                 end
%                 solOut = run_single_step(theta_vals(i),Fr_vals(i),tau0_vals(i),wlen_vals(i-1), alpha_vals(i-1), d_vals(i-1), ...
%                         rf_vals(i-1),h_pert,xi_in,y_in);
%             end
            if solOut.stats.maxres < tol
                h_wave = solOut.y(4,:);
                h_diff = max(h_wave)-min(h_wave);
                if h_diff>1e-4
                    y_in = solOut.y;
                    xi_in = solOut.x;
                else
                    [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_vals(i),3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3), linspace(tau0_vals(i-1),tau0_vals(i),3)...
                    , linspace(wlen_vals(i-1),wlen_vals(i),3),  linspace(alpha_vals(i-1),alpha_vals(i),3), ...
                    linspace(d_vals(i-1),d_vals(i),3), linspace(rf_vals(i-1),rf_vals(i),3), ...
                    xi_in, y_in, tol,counter+1);
                end
            else  
                    [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_vals(i),3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3), linspace(tau0_vals(i-1),tau0_vals(i),3), ...
                    linspace(wlen_vals(i-1),wlen_vals(i),3),  linspace(alpha_vals(i-1),alpha_vals(i),3), ...
                    linspace(d_vals(i-1),d_vals(i),3),...
                    linspace(rf_vals(i-1),rf_vals(i),3), xi_in, y_in, tol, counter+1);
            end
        end
        
        function [solOut] = run_single_step(theta_in,Fr_in,tau0_in,wlen_in,alpha_in,d_in,rf_in,h_pert,xi_in,y_in)
            [h0,~] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in,0);
            p_tot = rho*g*cosd(theta_in);
            crit_pb = rho_f*g*cosd(theta_in)*h0;
            u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);

            z_scale = h0;
            v_scale = u_eq;
            p_scale = crit_pb;
            t_scale = z_scale/v_scale;

            tau0_dl = tau0/p_scale;
            eta_f_dl = eta_f/(p_scale*t_scale);
            alpha_dl = alpha_in*p_scale;
            g_dl = g*t_scale/v_scale; 

            u_eq_dl = u_eq/v_scale;
            p_tot_grad_dl = p_tot/p_scale*z_scale;
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale; 
            rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
            d_dl = d_in/z_scale;

            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            
            tau0_dl = tau0_in/p_scale;
            
            if set_h_min
                h_min_in = wlen_in;
            else
                lambda_in = wlen_in;
            end
            stat_reg = (size(y_in,1) == 14);

            opts = bvpset('RelTol',tol);
            try
                if stat_reg
                    solInit1=bvpinit(xi_in,@bvp_guess_static);
                    solOut = bvp4c(@viscous_syst_static,@bc_vals,solInit1,opts);
                else
                    solInit1=bvpinit(xi_in,@bvp_guess);
                    solOut = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
                end
                catch ME
                    switch ME.identifier
                        case 'MATLAB:bvp4c:SingJac'
                            warning('Singular Jacobian encountered, reducing step size');
                            solOut = -1;
                        otherwise
                            rethrow(ME)
                    end    
            end
        
            function guess = bvp_guess(xi)
                % Initial guess function from the ode soln
                guess_ind = sum(xi_in<xi)+1;
                if guess_ind == max(size(xi_in))
                    guess = y_in(:,end);
                else
                    gap = xi_in(guess_ind+1)-xi_in(guess_ind);
                    guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
                end  
            end

            function dydxi = viscous_syst(xi,y)
                % The system from the viroulet paper with a different friction
                % law. Need to solve for u_w and Q1 too.
                u_w = y(1);
                lambda = y(2);
                h_crit_pos = y(3);
                
                Q1_l = y(4);
                h_l = y(5);
                u_l = (-Q1_l + h_l.*u_w)./h_l;
                m_l = y(6);
                phi_l = y(7)./Q1_l;
                pb_l = y(8) + g_dl*cosd(theta_in)*rho_dl*chi.*h_l;
                
                zeta_l = 3/(2*alpha_dl*h_l) + P/4;
                p_p_l = p_tot_grad_dl*h_l-pb_l;
                Iv_l = 3*eta_f_dl*abs(u_l)/h_l/p_p_l;
                D_l = -2/beta_dl/h_l*(pb_l-h_l);

                denom_l = (h_l.^3/Fr_in^2-Q1_l.^2);
                fric_coeff_l = p_p_l/(p_tot_grad_dl*h_l)*mu_Iv_fn(Iv_l);
                fb_val_l = tand(theta_in)-fric_coeff_l-tau0_dl*rho_f/rho/h_l;
                dhdxi_l = 1/Fr_in^2.*h_l.^3.*(fb_val_l+P*D_l*h_l)./denom_l;

                R_w3_l = -phi_l*rho_f_dl/rho_dl*D_l;
                R_w4_l = (-P*chi+zeta_l)*D_l - 2*3/alpha_dl/h_l*u_l*(phi_l - phi_c./(1+sqrt(Iv_l)));

                dQdxi_l = -P*D_l;
                dmdxi_l = h_l*u_l^(1-pres_h)/(lambda);
                dy6dxi_l = -R_w3_l;
                dy7dxi_l = R_w4_l/(u_l-u_w);
                dydxi_l = [dQdxi_l,dhdxi_l,dmdxi_l,dy6dxi_l,dy7dxi_l]*h_crit_pos;

                Q1_r = y(9);
                h_r = y(10);
                u_r = (-Q1_r + h_r.*u_w)./h_r;
                m_r = y(11);
                phi_r = y(12)./Q1_r;
                pb_r = y(13) + g_dl*cosd(theta_in)*rho_dl*chi.*h_r;

                zeta_r = 3/(2*alpha_dl*h_r) + P/4;
                p_p_r = p_tot_grad_dl*h_r-pb_r;
                Iv_r = 3*eta_f_dl*abs(u_r)/h_r/p_p_r;
                D_r = -2/beta_dl/h_r*(pb_r-h_r);

                denom_r = (h_r.^3/Fr_in^2-Q1_r.^2);
                fric_coeff_r = p_p_r/(p_tot_grad_dl*h_r)*mu_Iv_fn(Iv_r);
                fb_val_r = tand(theta_in)-fric_coeff_r-tau0_dl*rho_f/rho/h_r;
                dhdxi_r = 1/Fr_in^2.*h_r.^3.*(fb_val_r+P*D_r*h_r)./denom_r;

                R_w3_r = -phi_r*rho_f_dl/rho_dl*D_r;
                R_w4_r = (-P*chi+zeta_r)*D_r - 2*3/alpha_dl/h_r*u_r*(phi_r - phi_c./(1+sqrt(Iv_r)));

                dQdxi_r = -P*D_r;
                dmdxi_r = h_r*u_r^(1-pres_h)/(lambda);
                dy6dxi_r = -R_w3_r;
                dy7dxi_r = R_w4_r/(u_r-u_w);
                dydxi_r = [dQdxi_r,dhdxi_r,dmdxi_r,dy6dxi_r,dy7dxi_r]*(lambda-h_crit_pos);

                dydxi = [0,0,0,dydxi_l,dydxi_r]';
            end

            function resid = bc_vals(ya,yb)
                h_crit = (yb(4)*Fr_in)^(2/3);
                base_resid = [ya(4)-yb(9), ya(6), yb(11)-rf_in,...
                    (ya(7)-yb(12)), (ya(8)-yb(13)), yb(4)-ya(9), yb(6)-ya(11), yb(7)-ya(12), yb(8)-ya(13), ...
                    yb(5)-(h_crit-h_pert),ya(10)-(h_crit+h_pert)];
                if set_h_min
                    h_start = max(ya(4)/ya(1)+h_b_pert,h_min_in);
                    h_stop = (-h_start + sqrt(h_start.^2+8*ya(4)^2./(h_start/Fr_in^2)))/2;
                    len_resid = [ya(5)-h_start, yb(10)-h_stop];
                else
                    h_stop = (-ya(5) + sqrt(ya(5).^2+8*ya(4)^2./(ya(5)/Fr_in^2)))/2;
                    len_resid = [ya(2)-lambda_in, yb(10)-h_stop];
                end
                resid = [base_resid,len_resid];
            end
            
            function dydxi = viscous_syst_static(xi,y)
                % The system from the viroulet paper with a different friction
                % law. Need to solve for u_w and Q1 too.
                u_w = y(1);
                lambda = y(2);
                h_crit_pos = y(3);
                static_len = y(14);
                
                Q1_l = y(4);
                h_l = y(5);
                u_l = (-Q1_l + h_l.*u_w)./h_l;
                m_l = y(6);
                phi_l = y(7)./Q1_l;
                pb_l = y(8) + g_dl*cosd(theta_in)*rho_dl*chi.*h_l;
                
                zeta_l = 3/(2*alpha_dl*h_l) + P/4;
                p_p_l = p_tot_grad_dl*h_l-pb_l;
                Iv_l = 3*eta_f_dl*abs(u_l)/h_l/p_p_l;
                D_l = -2/beta_dl/h_l*(pb_l-h_l);

                denom_l = (h_l.^3/Fr_in^2-Q1_l.^2);
                fric_coeff_l = p_p_l/(p_tot_grad_dl*h_l)*mu_Iv_fn(Iv_l);
                fb_val_l = tand(theta_in)-fric_coeff_l-tau0_dl*rho_f/rho/h_l;
                dhdxi_l = 1/Fr_in^2.*h_l.^3.*(fb_val_l+P*D_l*h_l)./denom_l;

                R_w3_l = -phi_l*rho_f_dl/rho_dl*D_l;
                R_w4_l = (-P*chi+zeta_l)*D_l - 2*3/alpha_dl/h_l*u_l*(phi_l - phi_c./(1+sqrt(Iv_l)));

                dQdxi_l = -P*D_l;
                dmdxi_l = h_l*u_l^(1-pres_h)/(lambda);
                dy6dxi_l = -R_w3_l;
                dy7dxi_l = R_w4_l/(u_l-u_w);
                dydxi_l = [dQdxi_l,dhdxi_l,dmdxi_l,dy6dxi_l,dy7dxi_l]*(h_crit_pos-static_len);

                Q1_r = y(9);
                h_r = y(10);
                u_r = (-Q1_r + h_r.*u_w)./h_r;
                m_r = y(11);
                phi_r = y(12)./Q1_r;
                pb_r = y(13) + g_dl*cosd(theta_in)*rho_dl*chi.*h_r;

                zeta_r = 3/(2*alpha_dl*h_r) + P/4;
                p_p_r = p_tot_grad_dl*h_r-pb_r;
                Iv_r = 3*eta_f_dl*abs(u_r)/h_r/p_p_r;
                D_r = -2/beta_dl/h_r*(pb_r-h_r);

                denom_r = (h_r.^3/Fr_in^2-Q1_r.^2);
                fric_coeff_r = p_p_r/(p_tot_grad_dl*h_r)*mu_Iv_fn(Iv_r);
                fb_val_r = tand(theta_in)-fric_coeff_r-tau0_dl*rho_f/rho/h_r;
                dhdxi_r = 1/Fr_in^2.*h_r.^3.*(fb_val_r+P*D_r*h_r)./denom_r;

                R_w3_r = -phi_r*rho_f_dl/rho_dl*D_r;
                R_w4_r = (-P*chi+zeta_r)*D_r - 2*3/alpha_dl/h_r*u_r*(phi_r - phi_c./(1+sqrt(Iv_r)));

                dQdxi_r = -P*D_r;
                dmdxi_r = h_r*u_r^(1-pres_h)/(lambda);
                dy6dxi_r = -R_w3_r;
                dy7dxi_r = R_w4_r/(u_r-u_w);
                dydxi_r = [dQdxi_r,dhdxi_r,dmdxi_r,dy6dxi_r,dy7dxi_r]*(lambda-h_crit_pos);

                dydxi = [0,0,0,dydxi_l,dydxi_r,0]';
            end

            function resid = bc_vals_static(ya,yb)
                h_crit = (yb(4)*Fr_in)^(2/3);
                resid = [ya(4)-yb(9), ya(6), yb(11)-rf_in,...
                    (ya(7)-yb(12)), (ya(8)-yb(13)), yb(4)-ya(9), yb(6)-ya(11), yb(7)-ya(12), yb(8)-ya(13), ...
                    yb(5)-(h_crit-h_pert),ya(10)-(h_crit+h_pert),ya(2)-lambda_in, yb(10)-h_stop,ya(5)-h_start];
            end
        end
    end
end