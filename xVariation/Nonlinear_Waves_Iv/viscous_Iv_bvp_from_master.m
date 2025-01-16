function [xi_final,y_final] = viscous_Iv_bvp_from_master(specify_param,params,provide_init,master_xi,master_y,master_params)
% Converts the master wave stored in master_wave_no_pe.txt into a waveform
% for the conditions specified. Allows change in theta, lambda, Froude
% number and viscosity nu.

% Can be run with no inputs and just with the parameter values set in the
% code

% Can be run with specify_param set to true and with params as a 6 entry
% list [Fr,theta,lambda,nu,alpha,d]. In this case runs from master wave to
% specified values.

% Can be run with specify_param and provide_init set to true. This case
% needs params as the target 6 value list, master_xi as the initial value
% xi, master_y as the initial y and master_params as the initial parameters
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;

    reg_param = 1*10^7;

    phi_c=0.585; % Volume fraction

    g=9.81; % m/s^2
    

%     eta_f = 1.18e-5;
%     rho_f = 1;
    
    rho_f = 1000;
    eta_f = 0.0010016; % Pa s
    
    rho_p = 2500;
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    
    if ~exist("specify_param","var")
        specify_param = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    if ~provide_init
        master_name = "master_wave_no_pe.txt"; %master_wave_no_pe
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        master_wave_type = record.wave_type(in_table);
        master_theta = record.theta(in_table); 
        master_lambda = record.lambda(in_table);
        master_Fr = record.Fr(in_table);
        master_nu = record.nu(in_table);
        master_tau0 = record.tau0(in_table);
        master_rf = 1;
        master_pres_h = strcmp(master_wave_type,"no_pe_pres_h");
            
    else
        master_cell = num2cell(master_params);
        if (size(master_cell,2) == 5)
            [master_Fr,master_theta,master_lambda,master_nu,master_tau0] = master_cell{:};
            master_rf = 1;
            master_pres_h = 0;
        elseif (size(master_cell,2) == 7)
            [master_Fr,master_theta,master_lambda,master_nu,master_tau0,master_rf,master_pres_h] = master_cell{:};
        end
    end
    
    if ~specify_param
        Fr_eq = 1.0; 
        lambda = 10;
        theta = 12;
        nu_ratio = 1.0;
        tau0 = 0;
        rel_flux = 1;
        pres_h = 1;
%         h0 = 0.1;
        [h0, Iv_final] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr_eq, tau0);
        filename = "lambda_10_no_pe.txt";
    else
        param_cell = num2cell(params);
        if (size(param_cell,2) == 5)
            [Fr_eq,theta,lambda,nu_ratio,tau0] = param_cell{:};
            rel_flux = 1;
            pres_h = 0;
        elseif (size(param_cell,2) == 7)
            [Fr_eq,theta,lambda,nu_ratio,tau0,rel_flux,pres_h] = param_cell{:};
        end
    end
    
    if (pres_h ~= master_pres_h)
        m_init = zeros(size(master_xi));
        m_val = 0;
        u_flux = (master_y(1,1) - master_y(2,:)./master_y(3,:)).^(1-pres_h);
        for j = 1:size(master_xi,2)
            m_init(j) = m_val;
            if j < size(master_xi,2)
                m_val = m_val + 1/master_lambda* (master_y(3,j)*u_flux(j)+master_y(3,j+1)*u_flux(j+1))/2*(master_xi(j+1)-master_xi(j));
            end
        end
        master_y = vertcat(master_y(1:4,:),m_init);
        master_rf = m_val;
    end
    
    nu_final = 3/4*mu_Iv_fn(Iv_final)*eta_f/Iv_final/rho*nu_ratio;
    stable = viscous_stability(theta,Fr_eq,nu_final,lambda);
    if stable
        error("Conditions are stable to perturbation so no waves will occur")
    else 
        Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
        nu_ratio_in = max(nu_ratio/master_nu,master_nu/nu_ratio);
        theta_ratio = max(theta/master_theta,master_theta/theta);
        lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
        tau0_ratio = abs(tau0-master_tau0)/(max(master_tau0,tau0)+1e-6);
        rf_ratio = abs(rel_flux - master_rf);
        ratio_sum = Fr_ratio-1+nu_ratio_in-1+20*(theta_ratio-1)+(lambda_ratio-1)+tau0_ratio+rf_ratio;
        if ratio_sum > 1e-6
            n_step = max(min(ceil(30*ratio_sum),400),90);
        else
            n_step=2;
        end
        % If the parameters are the same as the initial values run 2 steps, if
        % not, run more

        lambda_list = linspace(master_lambda, lambda, n_step);
        Fr_list = linspace(master_Fr,Fr_eq,n_step);
        nu_list = linspace(master_nu,nu_ratio,n_step);
        theta_list = linspace(master_theta,theta,n_step);
        tau0_list = linspace(master_tau0,tau0,n_step);
        rf_list = linspace(master_rf, rel_flux, n_step);
        
        end_diff = zeros(size(lambda_list));
        [xi_final,y_final] = run_bvp_step(Fr_list, nu_list, theta_list, lambda_list, tau0_list, rf_list, master_xi, master_y);
        [xi_final,y_final] = run_bvp_step([Fr_eq,Fr_eq], [nu_ratio,nu_ratio], [theta,theta], [lambda,lambda], [tau0,tau0], [rel_flux,rel_flux], xi_final,y_final,1e-6);
        h_final = y_final(3,:);
%         plot(xi_final,h_final)
        if ~specify_param
            out_vec = vertcat(xi_final,y_final);
            save("Results/"+filename,"out_vec","-ascii")
            if pres_h
                type = "no_pe_pres_h";
            else
                type = "no_pe";
            end
            write_record("Results/wave_record.csv",filename,{type,"water",Fr_eq,theta,lambda,nu_ratio,0,0,tau0})
        end
    end
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, nu_vals, theta_vals, lambda_vals, tau0_vals, rf_vals, xi_in, y_in, tol,counter)
        % Can account for a change in wavelength but should really use
        % viscous_Iv_bvp_from_master for that.
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-4;
        end
        if counter > 10
            out_vec = vertcat(xi_in,y_in);
            fail_name = filename;
            save("Results/"+fail_name,"out_vec","-ascii")
            write_record("Results/wave_record.csv",fail_name,{"no_pe","water",Fr_vals(1),theta_vals(1),lambda_vals(1),nu_vals(1),0,0,tau0_vals(1)})
            error("Max iteration depth reached, non convergence")
        end
        
        nstep = size(lambda_vals,2);
        for i = 2:nstep
            theta_in = theta_vals(i);
            Fr_in = Fr_vals(i);
            tau0_in = tau0_vals(i);
            rf_in = rf_vals(i);
            [h0, crit_Iv] = crit_Iv_tau0(theta_in, rho_p, rho_f, eta_f, Fr_in, tau0_in);
            u_eq = Fr_in*sqrt(g*cosd(theta_in)*h0);
            nu = 3/4*mu_Iv_fn(crit_Iv)*eta_f/crit_Iv/rho*nu_vals(i);
            nu_dl = nu/(u_eq*h0);
            tau0_dl = tau0_in/(rho_f*g*cosd(theta_in)*h0);
            lambda_old = lambda_vals(i-1);
            lambda_in = lambda_vals(i);
            xi_in = xi_in/lambda_old*lambda_in;
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
                h_wave = solN1.y(3,:);
                h_diff = max(h_wave)-min(h_wave);
                if h_diff>1e-4
                    y_in = solN1.y;
                    xi_in = solN1.x;
                    if counter == 1
                        end_diff(i)=lambda_in-xi_in(end);
                    end
                else
                    [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    ,linspace(lambda_vals(i-1),lambda_in,3), linspace(tau0_vals(i-1),tau0_in,3), ...
                    linspace(rf_vals(i-1),rf_in,3), xi_in/lambda_in*lambda_old, y_in, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    ,linspace(lambda_vals(i-1),lambda_in,3), linspace(tau0_vals(i-1),tau0_in,3), ...
                    linspace(rf_vals(i-1),rf_in,3), xi_in/lambda_in*lambda_old, y_in, tol, counter+1);
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
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);

            dhdxi = n;
            n_coeff = h/Fr_in^2-Q1^2/h^2;
            Fr = Fr_in*abs(u)/h;
            Iv = crit_Iv*abs(u)/h^2;
            force_bal = h*(tand(theta_in)-sign(u).*P*mu_Iv_fn(Iv))/Fr_in^2;
            n_eq = (force_bal)./n_coeff;
            
            dndxi = n^2/h+h/Q1/nu_dl*(n_coeff*n-force_bal);
            dmdxi = h*u^(1-pres_h)/lambda_in;
            dydxi = [0,0,dhdxi,dndxi,dmdxi]';
        end
        
        function resid = bc_vals(ya,yb)
           resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-rf_in]; 
        end
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end