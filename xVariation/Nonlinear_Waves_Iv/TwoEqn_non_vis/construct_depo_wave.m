function [xi_final,y_final] = construct_depo_wave(specify_param,params,provide_init,master_xi,master_y,master_params)
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
    
    if ~exist("specify_param","var")
        specify_param = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    if ~provide_init
        master_name = "no_pe_no_vis_master_pres_h.txt"; %;"tau0_low_very_long.txt"
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable("Results/wave_record.csv");

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table);
        master_Fr = record.Fr(in_table);
        master_tau0 = record.tau0(in_table);
        master_pres_h= strcmp(wave_type{1}(max(end-6,1):end),"_pres_h");
    else
        master_cell = num2cell(master_params);
        [master_Fr,master_theta,master_tau0,~,~,~,master_pres_h] = master_cell{:};
    end
    if size(master_y,1) == 6
        master_stat_len = master_y(6,1);
    else
        master_stat_len = 0;
    end
    master_lambda = master_y(3,1)+master_stat_len;
    if (master_y(4,1)-master_y(2,1)/master_y(1,1)<1e-2)
        master_h_min = -1;
    else
        master_h_min = master_y(4,1);
    end
    master_rf = master_y(5,end);
    
    if ~specify_param
        Fr_eq = 3.0; 
        theta = 12;
        tau0 = 0;
        stat_len=0;
        h_min = -1;
        rel_flux = 1;
        lambda_spec = 300;
        pres_h = true;
        set_h_min = false;
        % h_alt specifies the distance from the static minimum that the
        % wave must begin at.
        filename = "Fr3_theta12_long.txt";
        
    else
        param_cell = num2cell(params);
        [Fr_eq,theta,tau0,len_spec,set_h_min,rel_flux,pres_h] = param_cell{:};
        stat_len = 0;
        if set_h_min
            h_min = len_spec;
            lambda_spec = 1;
        else
            lambda_spec = len_spec;
            h_min = 1;
        end
    end
    
    if (pres_h ~= master_pres_h)
        m_init = zeros(size(master_xi));
        m_val = pres_h*master_stat_len/master_lambda;
        u_flux = (master_y(1,1) - master_y(2,1)./master_y(4,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + (master_y(4,j)*u_flux(j)+master_y(4,j+1)*u_flux(j))/master_lambda;
            end
        end
        if master_stat_len>0
            master_y = vertcat(master_y(1:4,:),m_init,master_y(6:end,:));
        else
            master_y = vertcat(master_y(1:4,:),m_init);
        end
        master_rf = m_val;
    end
%     stable = viscous_stability(theta,Fr_eq,nu,lambda);

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    sl_ratio = abs(stat_len-master_stat_len)/(max(stat_len,master_stat_len)+1e-6);
    hm_ratio = max(h_min/master_h_min,h_min/master_h_min);
    lambda_ratio = max(lambda_spec/master_lambda,master_lambda/lambda_spec);
    rf_ratio = abs(master_rf-rel_flux);

    ratio_sum = Fr_ratio-1+20*(theta_ratio-1)+tau0_ratio+(hm_ratio-1)*set_h_min+(1-set_h_min)*(lambda_ratio-1)+rf_ratio;
    if sl_ratio > 1e-6
        sl_step = max(min(ceil(30*sl_ratio),400),90);
    else
        sl_step = 2;
    end
    if ratio_sum > 1e-6
        n_step = max(min(ceil(30*ratio_sum),400),90);
    else
        n_step=max(2,sl_step);
    end
        % If the parameters are the same as the initial values run 2 steps, if
        % not, run more
    Fr_list = linspace(master_Fr,Fr_eq,n_step);
    theta_list = linspace(master_theta,theta,n_step);
    tau0_list = linspace(master_tau0,tau0,n_step);
    rf_list = linspace(master_rf,rel_flux,n_step);
%     sl_list = linspace(master_stat_len,stat_len,n_step);
    if set_h_min
        if master_h_min>0 && h_min<0
            len_list = linspace(master_h_min,0,n_step);
        else
            len_list = linspace(max(master_h_min,master_y(2,1)/master_y(1,1)),h_min,n_step);
        end
        if sl_ratio>0 
            if (master_h_min == -1 && h_min == -1)
               len_list = linspace(master_stat_len,stat_len,n_step);
               stat_reg = true;
            elseif (master_h_min ~= -1 && h_min == -1)
                master_stat_len = 0;
                len_list = horzcat(len_list,linspace(master_stat_len,stat_len,sl_step));
                stat_reg = false;
            elseif ((master_h_min == -1) && (h_min ~= -1))
                stat_len = 0;
                len_list = horzcat(linspace(master_stat_len,stat_len,sl_step-1),-1,len_list);
                stat_reg = true;
            else
                error("One of the start or end conditions must contain a static region ")
            end
        end
    else
        len_list = linspace(master_lambda,lambda_spec,n_step);
        stat_reg = size(master_y,1)==6;
    end
    
    h_b_pert = 5e-4;

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, tau0_list, rf_list, len_list, stat_reg, master_xi, master_y);
    [xi_final,y_final] = run_bvp_step([Fr_eq,Fr_eq], [theta,theta], [tau0,tau0], [rel_flux,rel_flux], [len_list(end),len_list(end)], size(y_final,1)==6, xi_final, y_final,1e-5);
%         plot(xi_final,h_final)
    if ~specify_param
        out_vec = vertcat(xi_final,y_final);
        save("Results/"+filename,"out_vec","-ascii")
        if pres_h
            sim_type = "no_pe_no_vis_pres_h";
        else
            sim_type = "no_pe_no_vis";
        end
        write_record("Results/wave_record.csv",filename,{sim_type,"water",Fr_eq,theta,0,0,tau0})
    end
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, tau0_vals, rf_vals, len_vals, stat_reg, xi_in, y_in, tol,counter)
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
            error('Daniel:MaxIt',"Max iteration depth reached, non convergence")
        end
        
        nstep = size(Fr_vals,2);
        for i = 2:nstep
            theta_in = theta_vals(i);
            Fr_in = Fr_vals(i);
            tau0_in = tau0_vals(i);
            rf_in = rf_vals(i);
            [h0,eq_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in,0);
            tau0_dl = tau0_in/(rho_f*g*cosd(theta_in)*h0);
            if set_h_min
                if ~stat_reg && len_vals(i)==0 && len_vals(i-1)~=0
                    stat_reg = true;
                    y_in = vertcat(y_in,zeros(size(xi_in)));
                elseif stat_reg && len_vals(i)==-1 && len_vals(i-1)==0
                    stat_reg = false;
                    y_in = y_in(1:5,:);
                end
            else
                if stat_reg 
                    if y_in(6,1) < 1e-6
                        stat_reg = false;
                        y_in = y_in(1:5,:);
                    end
                elseif ~stat_reg && (y_in(4,1)-y_in(2,1)/y_in(1,1))/y_in(4,1)^2*eq_Iv<2.5e-7
                    stat_reg = true;
                    y_in = vertcat(y_in,zeros(size(xi_in)));
                    opts = bvpset('RelTol',tol);
                    solInit1=bvpinit(xi_in,@bvp_guess);
                    solN1 = bvp4c(@viscous_syst_static,@bc_vals_static,solInit1,opts);
                    if resid < tol
                        y_in = solN1.y;
                        xi_in = solN1.x;
                    end
                end
            end
            
            if set_h_min
                if stat_reg
                    sl_spec = len_vals(i);
                else
                    h_min_in = len_vals(i);
                end
            else
                lambda_in = len_vals(i);
            end
            
%             u_const = eq_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta_in);
%             h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
%             u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
            
            

            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol);
            try
                if stat_reg
                    solN1 = bvp4c(@viscous_syst_static,@bc_vals_static,solInit1,opts);
                else 
                    solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
                end
                resid = solN1.stats.maxres;
                catch ME
                    switch ME.identifier
                        case 'MATLAB:bvp4c:SingJac'
                            warning('Singular Jacobian encountered, reducing step size');
                            resid = tol+1;
                        otherwise
                            rethrow(ME)
                    end    
            end
            % Solves the 5 value bvp for u_w, Q1, h, n and m.
            if resid < tol
                h_wave = solN1.y(4,:);
                h_diff = max(h_wave)-min(h_wave);
                if h_diff>1e-4
                    y_in = solN1.y;
                    xi_in = solN1.x;
                else
                    [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(rf_vals(i-1),rf_vals(i),3), ...
                    linspace(len_vals(i-1),len_vals(i),3),stat_reg,xi_in, y_in, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(rf_vals(i-1),rf_vals(i),3), ...
                    linspace(len_vals(i-1),len_vals(i),3),stat_reg,xi_in, y_in, tol, counter+1);
            end
        end
        y_out = solN1.y;
        xi_out = solN1.x;
        
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
            Q1 = y(2);
            lambda = y(3);
            h = y(4);
            u = (-Q1 + h.*u_w)./h;
            m = y(5);
            
            Iv = eq_Iv.*abs(u)/h.^2;
            
            denom = (h.^3/Fr_in^2-Q1.^2);
            if (abs(denom)>1e-6)
                fb_val = tand(theta_in)-sign(u).*mu_Iv_fn(Iv).*(rho-rho_f)/rho-sign(u).*tau0_dl*rho_f/rho/h;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_val./denom;
            else
                ud = Q1/h^2;
                Iv_deriv = eq_Iv.*2.*ud./h.^2-2.*Iv/h;
                fb_deriv = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_deriv/(3*h.^2/Fr_in^2);
                
            end

            dmdxi = h/(lambda)*u^(1-pres_h);
            dydxi = [0,0,0,dhdxi,dmdxi]'*lambda;
        end
        
    
        function resid = bc_vals(ya,yb)
            h_crit = (max(ya(2)*Fr_in,0))^(2/3);
            if set_h_min
                h_start = max(ya(2)/ya(1)+h_b_pert,h_min_in);
                len_res = ya(4)-h_start;
            else
                h_start = ya(4);
                len_res = ya(3)-lambda_in;
            end
            h_stop = (-h_start + sqrt(h_start.^2+8*ya(2)^2./(h_start/Fr_in^2)))/2;
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h0, tau0_in,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2; 
            resid = [ya(1)-(u_crit*h_crit + ya(2))/h_crit, len_res, yb(4)-h_stop, ya(5), yb(5)-rf_in]; 
        end
        
        function dydxi = viscous_syst_static(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = y(1);
            Q1 = y(2);
            lambda = y(3);
            h = y(4);
            u = (-Q1 + h.*u_w)./h;
            m = y(5);
            s_len = y(6);
            Iv = eq_Iv.*abs(u)/h.^2;
            
            denom = (h.^3/Fr_in^2-Q1.^2);
            if (abs(denom)>1e-6)
                fb_val = tand(theta_in)-sign(u).*mu_Iv_fn(Iv).*(rho-rho_f)/rho-sign(u).*tau0_dl*rho_f/rho/h;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_val./denom;
            else
                ud = Q1/h^2;
                Iv_deriv = eq_Iv.*2.*ud./h.^2-2.*Iv/h;
                fb_deriv = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_deriv/(3*h.^2/Fr_in^2);
                
            end
            
            dmdxi = h/(lambda+s_len)*u^(1-pres_h);
            dydxi = [0,0,0,dhdxi,dmdxi,0]'*lambda;
        end
        
    
        function resid = bc_vals_static(ya,yb)
            h_crit = (max(ya(2)*Fr_in,0))^(2/3);
            if set_h_min
                len_res = ya(6)-sl_spec;
            else
                len_res = ya(3)+ya(6)-lambda_in;
            end
            h_start = ya(2)/ya(1)+h_b_pert;
            h_stop = (-h_start + sqrt(h_start.^2+8*ya(2)^2./(h_start/Fr_in^2)))/2;
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h0, tau0_in,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2; 
            resid = [ya(1)-(u_crit*h_crit + ya(2))/h_crit, len_res, ya(4)-h_start, ...
                yb(4)-h_stop, ya(5)-ya(6)/lambda_in*ya(2)/ya(1), yb(5)-rf_in]; 
        end
    end
end