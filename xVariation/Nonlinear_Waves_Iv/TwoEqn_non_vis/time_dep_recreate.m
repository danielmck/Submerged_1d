 function [xi_final,y_final] = time_dep_recreate(h0_dim,theta,lambda_dim,tau0,rel_flux,pres_h)
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
    
    [Fr_eq,~] = crit_Iv_tau0_h(theta, rho_p, rho_f, eta_f, h0_dim, tau0,0);
    lambda = lambda_dim/h0_dim;
    u_eq_dim = Fr_eq*sqrt(g*cosd(theta)*h0_dim);
    
    if ~exist("rel_flux","var")
        rel_flux=1;
    end
    
    if ~exist("pres_h","var")
        pres_h=false;
    end
    
    Res_dir = "~/Documents/MATLAB/1D_System/xVariation/Nonlinear_Waves_Iv/NonViscous_bvp/";
    master_name = "wave_match.txt";
    master_file = load(Res_dir+"Results/"+master_name);
    master_xi = master_file(1,:);
    master_y = master_file(2:end,:);
    record = readtable(Res_dir+"Results/wave_record.csv");

    in_table = strcmp(record.Name, master_name);
    wave_type = record.wave_type(in_table);
    master_theta = record.theta(in_table);
    master_Fr = record.Fr(in_table);
    master_tau0 = record.tau0(in_table);
    master_lambda = master_y(3,1);
    master_rf = master_y(5,end);
    master_pres_h=strcmp(wave_type(max(end-7,1):end),"_pres_h");
    
    master_y = vertcat(master_y(1:2,:),master_y(4:end,:));
    
    if (pres_h ~= master_pres_h)
        m_init = zeros(size(master_xi));
        m_val = 0;
        u_flux = (master_y(1,1) - master_y(2,:)./master_y(3,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + 1*(master_y(3,j)*u_flux(j)+master_y(3,j+1)*u_flux(j+1))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        master_y = vertcat(master_y(1:3,:),m_init);
        master_rf = m_val;
    end

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
    rf_ratio = abs(rel_flux - master_rf);
    
    ratio_sum = Fr_ratio-1+20*(theta_ratio-1)+tau0_ratio+lambda_ratio+rf_ratio;
    if ratio_sum > 1e-6
        n_step = max(min(ceil(30*ratio_sum),400),90);
    else
        n_step=2;
    end
        % If the parameters are the same as the initial values run 2 steps, if
        % not, run more
    Fr_list = linspace(master_Fr,Fr_eq,n_step);
    theta_list = linspace(master_theta,theta,n_step);
    tau0_list = linspace(master_tau0,tau0,n_step);
    lambda_list = linspace(master_lambda, lambda, n_step);
    rf_list = linspace(master_rf, rel_flux, n_step);

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, tau0_list, lambda_list, rf_list, master_xi, master_y);
    if size(y_final,1) == 5
        stat_dist = y_final(5,1);
    else
        stat_dist = 0;
    end
    xi_final = (stat_dist/lambda+xi_final*(1-stat_dist/lambda))*lambda_dim;
    y_final = [u_eq_dim,u_eq_dim*h0_dim,h0_dim,u_eq_dim*h0_dim]'.*y_final(1:4,:);
    
    if stat_dist > 0
        xi_final = [0,xi_final];
        y_final = horzcat(y_final(:,1),y_final);
    end
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, tau0_vals, lambda_vals, rf_vals, xi_in, y_in, static_part, tol,counter)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-5;
        end
        if ~exist('static_part','var')
            static_part = false;
        end
        if counter > 10
%             out_vec = vertcat(xi_in,y_in);
%             fail_name = filename;
%             save("Results/"+fail_name,"out_vec","-ascii")
%             write_record("Results/wave_record.csv",fail_name,{"no_pe","water",Fr_vals(1),theta_vals(1),lambda_vals(1),nu_vals(1),0,0,tau0_vals(1)})
            error("Max iteration depth reached, non convergence")
        end
        
        
        nstep = size(Fr_vals,2);
        for i = 2:nstep
            theta_in = theta_vals(i);
            Fr_in = Fr_vals(i);
            tau0_in = tau0_vals(i);
            lambda_in = lambda_vals(i);
            rel_flux_in = rf_vals(i);
            sl_in=0;
            
            [h0,eq_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in,0);
%             u_const = eq_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta_in);
%             h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
%             u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
            
            tau0_dl = tau0_in/(rho_f*g*cosd(theta_in)*h0);
            
            curr_u_w = y_in(1,1);
            curr_Q1 = y_in(2,1);
            curr_hmin = y_in(3,1);
            if (curr_hmin<curr_Q1/curr_u_w+3e-3) && (lambda_in>=lambda_vals(i-1)) && ~static_part
                static_part = true;
                y_in = vertcat(y_in,zeros(size(xi_in)));
            elseif (static_part)
                if (y_in(5,1) < 1e-2) && (lambda_in<=lambda_vals(i-1))
                    static_part = false;
                    y_in = y_in(1:4,:);
                end
            end
            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol);
            try
                if ~static_part
                    solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
                else
                    solN1 = bvp4c(@viscous_syst_static,@bc_vals_static,solInit1,opts);
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
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(lambda_vals(i-1),lambda_in,3), ...
                    linspace(rf_vals(i-1),rel_flux_in,3), xi_in, y_in, static_part, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(lambda_vals(i-1),lambda_in,3), ...
                    linspace(rf_vals(i-1),rel_flux_in,3), xi_in, y_in, static_part, tol, counter+1);
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
            h = y(3);
            u = (-Q1 + h.*u_w)./h;
            m = y(4);
            
            Iv = eq_Iv.*u/h.^2;
            
            denom = (h.^3/Fr_in^2-Q1.^2);
            if (abs(denom)>1e-6)
                fb_val = tand(theta_in)-mu_Iv_fn(Iv).*(rho-rho_f)/rho-tau0_dl*rho_f/rho/h;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_val./denom;
            else
                ud = Q1/h^2;
                Iv_deriv = eq_Iv.*2.*ud./h.^2-2.*Iv/h;
                fb_deriv = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_deriv/(3*h.^2/Fr_in^2);
            end

            dmdxi = h/(lambda_in)*u^(1-pres_h);
            dydxi = [0,0,dhdxi,dmdxi]'*lambda_in;
        end

        function resid = bc_vals(ya,yb)
            h_crit = (ya(2)*Fr_in)^(2/3);
            h_stop = (-ya(3) + sqrt(ya(3).^2+8*ya(2)^2./(ya(3)/Fr_in^2)))/2;
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h0, tau0_in,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2; 
            resid = [ya(1)-(u_crit*h_crit + ya(2))/h_crit, yb(3)-h_stop, ya(4), yb(4)-rel_flux_in]; 
        end
        
        function dydxi = viscous_syst_static(xi,y)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = y(1);
            Q1 = y(2);
            depo_point = y(5);
            h = y(3);
            u = (-Q1 + h.*u_w)./h;
            m = y(4);
            
            Iv = eq_Iv.*u/h.^2;
            
            denom = (h.^3/Fr_in^2-Q1.^2);
            if (abs(denom)>1e-6)
                fb_val = tand(theta_in)-mu_Iv_fn(Iv).*(rho-rho_f)/rho-tau0_dl*rho_f/rho/h;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_val./denom;
            else
                ud = Q1/h^2;
                Iv_deriv = eq_Iv.*2.*ud./h.^2-2.*Iv/h;
                fb_deriv = -dmudIv_fn(Iv).*Iv_deriv.*(rho-rho_f)/rho+tau0_dl*rho_f/rho./h.^2;
                dhdxi = 1/Fr_in^2.*h.^3.*fb_deriv/(3*h.^2/Fr_in^2);
            end

            dmdxi = h/(lambda_in)*u^(1-pres_h);
            dydxi = [0,0,dhdxi,dmdxi,0]'*(lambda_in-depo_point);
        end

        function resid = bc_vals_static(ya,yb)
            h_crit = (ya(2)*Fr_in)^(2/3);
            h_start = ya(2)/ya(1)+2e-3;
            h_stop = (-h_start + sqrt(h_start.^2+8*ya(2)^2./(h_start/Fr_in^2)))/2;
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h0, tau0_in,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2; 
            resid = [ya(1)-(u_crit*h_crit + ya(2))/h_crit, ya(3)-h_start, yb(3)-h_stop, ya(4), yb(4)-rel_flux_in]; 
        end
    end
end