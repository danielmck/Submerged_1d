function [xi_final,y_final] = full_depo_wave_var_rho(specify_param,params,provide_init,master_xi,master_y,master_params)
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
        Res_dir = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/VarRho_NonVis/";
        master_name = "var_rho_master_pres_h.txt"; %no_vis_better_master
        master_file = load(Res_dir+"Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable(Res_dir+"Results/full_record.csv");

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table);
        master_Fr = record.Fr(in_table);
        master_tau0 = record.tau0(in_table);
        master_d = record.d(in_table);
        master_alpha = record.alpha(in_table);
        master_pres_h= strcmp(wave_type(max(end-7,1):end),"_pres_h");
        master_rf = master_y(3,end);
        
        master_u_w = record.u_w(in_table);
        master_lambda = record.lambda(in_table);
        master_xi_crit = record.crit_xi(in_table);
        master_p = [master_u_w,master_lambda,master_xi_crit];
        
        if (master_y(2,1)-master_y(1,1)/master_u_w<1e-2)
            master_h_min = -1;
        else
            master_h_min = master_y(5,1);
        end
        if (size(master_y,1) == 10)
            master_stat_len = master_y(10,1);
        else
            master_stat_len = 0;
        end
    else
        master_cell = num2cell(master_params);
        [master_Fr,master_theta,master_alpha,master_d,master_tau0,master_pres_h] = master_cell{:};
    end
    
    if ~specify_param
        Fr_eq = 0.8; 
        theta = 9;
        tau0 = 0;
        
        d = 1e-4;
        alpha=1e-5;
        
        set_h_min = false;
        h_min = master_h_min;
        lambda_final = 12;
        
        pres_h = true;
        rel_flux = 1;

        % h_alt specifies the distance from the static minimum that the
        % wave must begin at.
        filename = "theta_9_pres_h.txt";      
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
        m_init = zeros(size(master_xi));
        m_val = pres_h*master_stat_len/master_lambda*(master_lambda-master_xi_crit)*master_y(1,1)/master_u_w;
        u_flux = (master_u_w(1,1) - master_y(2,:)./master_y(1,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init(j) = m_val;
            if j < size(master_xi,2)
                scale = (master_xi(j)>1)*(master_lambda-master_xi_crit)+(master_xi(j)<1)*master_xi_crit;
                m_val = m_val + scale/master_lambda* (master_y(2,j)*u_flux(j)^(1-pres_h)+master_y(2,j+1)*u_flux(j+1)^(1-pres_h))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        master_y = vertcat(master_y(1:2,:),m_init,master_y(4:end,:));
        master_rf = m_val;
    end
    
    alt_xi = master_xi([true,diff(master_xi)>0]);
    alt_y = master_y(:,[true,diff(master_xi)>0]);
    break_ind = sum(alt_xi<1);
    xi_prep = [alt_xi(1:break_ind),1,alt_xi(break_ind+1:end)];
    y_prep = horzcat(alt_y(:,1:break_ind),alt_y(:,break_ind+1),alt_y(:,break_ind+1:end));

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    d_ratio = max(d/master_d,master_d/d);
    alpha_ratio = max(alpha/master_alpha,master_alpha/alpha);
    
    hm_ratio = min(max(h_min/master_h_min,master_h_min/h_min),1e3);
    lambda_ratio = min(max(lambda_final/master_lambda,master_lambda/lambda_final),1e3);
    rf_ratio = abs(rel_flux - master_rf);
    ratio_sum = Fr_ratio-1+20*(theta_ratio-1)+tau0_ratio+hm_ratio*set_h_min+(1-set_h_min)*(lambda_ratio-1)+rf_ratio+(d_ratio-1)+(alpha_ratio-1);
    if ratio_sum > 1e-4
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
    
    h_b_pert = 1e-3;
    denom_eps = 1e-4;
    
%     master_p = [master_y(1:3,1)',master_y(9,1)];
% %     master_p(4) = 0.05;
%     master_y = master_y(4:8,:);
    
    [xi_final,y_final,p_final] = run_bvp_step(Fr_list, theta_list, tau0_list, wlen_list, alpha_list, d_list, rf_list, xi_prep, y_prep, master_p);
    if ~specify_param
        out_vec = vertcat(xi_final,y_final);
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/full_record.csv",filename,{"full_no_vis","water",Fr_eq,theta,alpha,d,tau0,p_final(1),p_final(2),p_final(3)})
    end
    
    function [xi_in, y_in,p_in] = run_bvp_step(Fr_vals, theta_vals, tau0_vals, wlen_vals, alpha_vals, d_vals, rf_vals, xi_in, y_in, p_in, tol,counter)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-3;
        end
        xi_eps = 1e-3;
        if counter > 10
            out_vec = vertcat(xi_in,y_in);
            fail_name = "no_vis_better_pres_h.txt";
            save("Results/"+fail_name,"out_vec","-ascii")
            write_record("Results/full_record.csv",fail_name,{"full","water",Fr_vals(1),theta_vals(1),alpha_vals(1),d_vals(1),tau0_vals(1),p_in(1),p_in(2),p_in(3)})
            error("Max iteration depth reached, non convergence")
        end
        
%         static_part = (size(y_in,1) == 14);
        nstep = size(Fr_vals,2);
        n_min = 2 - (counter==1);
        for i = n_min:nstep
%             if (y_in(5,1)<y_in(4,1)/y_in(1,1)+1.5*h_b_pert) && (lambda_vals(i)>=lambda_vals(i-1)) && ~static_part
%                 static_part = true;
%                 y_in = vertcat(y_in,zeros(size(xi_in)));
%             elseif (static_part)
%                 if (y_in(end,1) < 1e-2) && (lambda_vals(i)<=lambda_vals(i-1))
%                     static_part = false;
%                     y_in = y_in(1:end-1,:);
%                 end
%             end
            theta_in = theta_vals(i);
            Fr_in = Fr_vals(i);
            tau0_in = tau0_vals(i);
            wlen_in = wlen_vals(i);
            alpha_in = alpha_vals(i);
            d_in = d_vals(i);
            rf_in = rf_vals(i);
            [h0,~] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in,0,true);
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
            rho_f_dl = rho_f*v_scale^2/p_scale;
            rho_p_dl = rho_p*v_scale^2/p_scale; 
            d_dl = d_in/z_scale;

            beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
            
            tau0_dl = tau0_in/p_scale;
            
            if set_h_min
                h_min_in = wlen_in;
            else
                lambda_in = wlen_in;
            end
            stat_reg = (size(y_in,1) == 10);

            opts = bvpset('RelTol',tol);
            try
                if stat_reg
                    solInit1=bvpinit(xi_in,@bvp_guess_static,p_in);
                    solOut = bvp4c(@viscous_syst_static,@bc_vals,solInit1,opts);
                else
                    solInit1=bvpinit(xi_in,@bvp_guess,p_in);
                    opts = bvpset('RelTol',tol,'NMax',2000*counter);
                    solOut = bvp4c(@(x,y,r,p) viscous_syst(x,y,r,p),@bc_vals,solInit1,opts);
                    
%                     y_in = solOut_np.y;
%                     xi_in = solOut_np.x;
%                     
%                     p_in = [y_in(1:3,1)',y_in(9,1)];
%                     y_in = y_in(4:8,:);
% 
%                     solInit1=bvpinit(xi_in,@bvp_guess,p_in);
%                     opts = bvpset('RelTol',tol,'NMax',3000*counter);
%                     solOut = bvp4c(@(x,y,r,p) viscous_syst(x,y,r,p),@bc_vals,solInit1,opts);
                end
                resid = solOut.stats.maxres;
                catch ME
                    switch ME.identifier
                        case 'MATLAB:bvp4c:SingJac'
                            warning('Singular Jacobian encountered, reducing step size');
                            resid = tol+1;
                        otherwise
                            rethrow(ME)
                    end    
            end
            if resid < tol
                h_wave = solOut.y(5,:);
                h_diff = max(h_wave)-min(h_wave);
                if h_diff>1e-4
                    y_in = solOut.y;
                    xi_in = solOut.x;
                    p_in = solOut.parameters;
                else
                    [xi_in,y_in,p_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_vals(i),3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3), linspace(tau0_vals(i-1),tau0_vals(i),3)...
                    , linspace(wlen_vals(i-1),wlen_vals(i),3),  linspace(alpha_vals(i-1),alpha_vals(i),3), ...
                    linspace(d_vals(i-1),d_vals(i),3), linspace(rf_vals(i-1),rf_vals(i),3), ...
                    xi_in, y_in, p_in, tol,counter+1);
                end
            else  
                    [xi_in,y_in,p_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_vals(i),3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3), linspace(tau0_vals(i-1),tau0_vals(i),3), ...
                    linspace(wlen_vals(i-1),wlen_vals(i),3),  linspace(alpha_vals(i-1),alpha_vals(i),3), ...
                    linspace(d_vals(i-1),d_vals(i),3),...
                    linspace(rf_vals(i-1),rf_vals(i),3), xi_in, y_in, p_in, tol, counter+1);
            end
        end 
        
        function guess = bvp_guess(xi,region)
            % Initial guess function from the ode soln
            guess_ind = sum(xi_in<xi)+1;
            if guess_ind == max(size(xi_in))
                guess = y_in(:,end);
            else
                gap = xi_in(guess_ind+1)-xi_in(guess_ind);
                if abs(gap)<1e-8
                    guess = y_in(:,guess_ind);
                else
                    guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
                end
            end  
        end
    
        function dydxi = viscous_syst(xi,y,region,p)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = p(1);
            lambda_l = p(2);
            h_crit_pos = p(3);
            
            Q1 = y(1);
            h = y(2);
            u = (-Q1 + h.*u_w)./h;
            m = y(3);
            phi = y(4)./Q1;
            rho_dl = (rho_p_dl*phi+rho_f_dl*(1-phi));
            P = (rho_dl-rho_f_dl)/rho_dl;
            chi = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            pb = y(5) + g_dl*cosd(theta_in)*rho_dl*chi.*h;
            p_tot_grad_dl = rho_dl*g_dl*cosd(theta_in);

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            D = -2/beta_dl/h*(pb-h);

            denom = (h.^3/Fr_in^2-Q1.^2);
            fric_coeff = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
            fb_val = tand(theta_in)-fric_coeff-tau0_dl*rho_f_dl/rho_dl/h;
            dhdxi = 1/Fr_in^2.*h.^3.*(fb_val+P*D*h)./denom;

            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P/4+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dQdxi = -P*D;
            dmdxi = h/(lambda_l)*u^(1-pres_h);
            dy6dxi = R_w3;
            dy7dxi = R_w4/(u-u_w);
            dydxi = [dQdxi,dhdxi,dmdxi,dy6dxi,dy7dxi];
            if abs(denom)<denom_eps
                A_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl.*h)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./p_p;
                A_pb_coeff = -P.*2./beta_dl-A_pp_coeff;
                A_h_coeff = p_p./(p_tot_grad_dl.*h.^2).*mu_Iv_fn(Iv)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*3*eta_f_dl.*abs(u)./h.^2./p_p+tau0_dl.*rho_f_dl./rho_dl./h.^2-P.*2./beta_dl+A_pp_coeff.*rho_dl.*g_dl.*cosd(theta_in);

                full_h_coeff = (h.^3./Fr_in^2.*(g_dl*cosd(theta_in)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr_in^2.*(fb_val+P.*D.*h));
                numer_h_term = full_h_coeff.*dhdxi;
                numer_other_term = h.^3./Fr_in^2.*(dy7dxi.*A_pb_coeff+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h.^2);

                h_quad_roots = roots([3*h.^2./Fr_in^2, -full_h_coeff-2.*Q1.*dQdxi, -numer_other_term]);
                if imag(h_quad_roots(1))==0
                    dydxi(2) = min(h_quad_roots(h_quad_roots>0));
                else
                    dydxi(2) = 0;
                end
            end
            switch region
                case 1    % x in [0 1]
                    dydxi = dydxi*h_crit_pos;
                case 2    % x in [1 2]
                    dydxi = dydxi*(lambda_l-h_crit_pos);
            end
            dydxi = dydxi';
        end
        
        function resid = bc_vals(ya,yb,p)
            h_mid = ya(2,2);
            u_mid = (-ya(1,2) + h_mid.*p(1))./h_mid;
            phi_mid = ya(4,2)/ya(1,2);
            rho_mid = (rho_p_dl*phi_mid+rho_f_dl*(1-phi_mid));
            P_mid = (rho_mid-rho_f_dl)/rho_mid;
            chi_mid = (rho_f_dl+3*rho_mid)/(4*rho_mid);
            denom_mid = (h_mid.^3/Fr_in^2-ya(1,2).^2);
            pb_mid = ya(5,2) + g_dl*cosd(theta_in)*rho_mid*chi_mid.*h_mid;
            p_tot_grad_mid = rho_mid*g_dl*cosd(theta_in);
            p_p_mid = p_tot_grad_mid*h_mid-pb_mid;
            Iv_mid = 3*eta_f_dl*abs(u_mid)/h_mid/p_p_mid;
            fric_coeff_mid = p_p_mid/(p_tot_grad_mid*h_mid)*mu_Iv_fn(Iv_mid);
            D_mid = -2/beta_dl/h_mid*(pb_mid-h_mid);
            fb_val_mid = tand(theta_in)-fric_coeff_mid-tau0_dl*rho_f_dl/rho_mid/h_mid+P_mid*D_mid*h_mid;
            if set_h_min
                h_start = max(yl(1)/p(1)+h_b_pert,h_min_in);
                h_stop = (-h_start + sqrt(h_start.^2+8*yl(1)^2./(h_start/Fr_in^2)))/2;
                len_resid = [ya(2,1)-h_start, yb(2,2)-h_stop];
            else
                h_stop = (-ya(2,1) + sqrt(ya(2,1).^2+8*ya(1,1)^2./(ya(2,1)/Fr_in^2)))/2;
                len_resid = [p(2)-lambda_in, yb(2,2)-h_stop];
            end
            
            cont_resid = (ya(:,2) - yb(:,1))';

            resid = [ya(1,1)-yb(1,2), ya(3,1), yb(3,2)-rf_in,...
                (ya(5,1)-yb(5,2)), denom_mid,len_resid, fb_val_mid, cont_resid];
        end
    end
end