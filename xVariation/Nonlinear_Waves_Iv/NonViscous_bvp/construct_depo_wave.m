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
        master_name = "no_pe_no_vis_master.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table);
        master_h_alt = record.h_alt(in_table);
        master_Fr = record.Fr(in_table);
        master_tau0 = record.tau0(in_table);
        master_stat_len = 0;
    else
        master_cell = num2cell(master_params);
        [master_Fr,master_theta,master_stat_len,~,master_tau0] = master_cell{:};
    end
    
    if ~specify_param
        Fr_eq = 1; 
        theta = 10;
        tau0 = 20;
        stat_len=0;
        h_alt = 0.2;
        % h_alt specifies the distance from the static minimum that the
        % wave must begin at.
        filename = "no_pe_no_vis_convert.txt";

%       A wave close to the conditions can be constructed from
%       depo_wave_ode in order to provide easier adaptation to the bvp.
%         [master_xi, h_init, Q1_in, u_w_in,flux] = depo_wave_ode(theta, Fr_eq, tau0);
%        u_init = (-Q1_in + h_init.*u_w_in)./h_init;
%     lambda_init = master_xi(end);
%     m_init = zeros(size(master_xi));
%     m_val = 0;
%     % Defines the average flux m that is needed to solve the bvp
%     for j = 2:size(master_xi,1)
%         m_init(j) = get_flux(master_xi(1:j),h_init(1:j),u_init(1:j));
%     end
%     bv = ones(size(master_xi));
%     master_xi=master_xi/lambda_init;
%     master_y = horzcat(bv*u_w_in,bv*Q1_in,bv*lambda_init,h_init,m_init);
        
    else
        param_cell = num2cell(params);
        [Fr_eq,theta,tau0] = param_cell{:};  
    end
    

%     stable = viscous_stability(theta,Fr_eq,nu,lambda);

    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
    sl_ratio = abs(stat_len-master_stat_len)/(max(stat_len,master_stat_len)+1e-6);
    ha_ratio = min(abs(log(h_alt/master_h_alt)),10);

    ratio_sum = Fr_ratio-1+20*(theta_ratio-1)+tau0_ratio+sl_ratio+ha_ratio;
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
    sl_list = linspace(master_stat_len,stat_len,n_step);
    ha_list = linspace(master_h_alt,h_alt,n_step);

    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, tau0_list, sl_list, ha_list, master_xi, master_y);
%         plot(xi_final,h_final)
        if ~specify_param
            out_vec = vertcat(xi_final,y_final);
            lambda = y_final(3,1);
            save("Results/"+filename,"out_vec","-ascii")
            write_record("Results/wave_record.csv",filename,{"no_pe_no_vis","water",Fr_eq,theta,h_alt,0,0,tau0})
        end
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, tau0_vals, sl_vals, ha_vals, xi_in, y_in, tol,counter)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-4;
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
            sl_in = sl_vals(i);
            h_alt_in = ha_vals(i);
            
            [h0,eq_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0,0);
%             u_const = eq_Iv/eta_f/3*(rho_p-rho_f)*g*phi_c*cosd(theta_in);
%             h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
%             u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
            
            tau0_dl = tau0_in/(rho_f*g*cosd(theta_in)*h0);

            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol);
            try
                solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
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
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(sl_vals(i-1),sl_in,3), ...
                    xi_in, y_in, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    , linspace(tau0_vals(i-1),tau0_in,3), linspace(sl_vals(i-1),sl_in,3), ...
                    xi_in, y_in, tol, counter+1);
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

            dmdxi = h/(lambda+sl_in)*u;
            dydxi = [0,0,0,dhdxi,dmdxi]'*lambda;
        end
        
    
        function resid = bc_vals(ya,yb)
            h_crit = (ya(2)*Fr_in)^(2/3);
            h_start = ya(2)/ya(1)+h_alt_in;
            h_stop = (-h_start + sqrt(h_start.^2+8*ya(2)^2./(h_start/Fr_in^2)))/2;
            [~, crit_Iv] = crit_Iv_tau0_h(theta_in, rho_p, rho_f, eta_f, h_crit*h0, tau0,0);
            u_crit = crit_Iv/eq_Iv*h_crit^2; 
            resid = [ya(1)-(u_crit*h_crit + ya(2))/h_crit, ya(4)-h_start, yb(4)-h_stop, ya(5), yb(5)-1]; 
        end
    end
end