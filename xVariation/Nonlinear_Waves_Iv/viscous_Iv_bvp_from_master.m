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
    if ~exist("specify_param","var")
        specify_param = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    if ~provide_init
        master_file = load("master_wave_no_pe.txt");
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        master_cond = readtable("master_wave_no_pe_cond.csv");
        master_Fr = master_cond.Fr_eq;
        master_theta = master_cond.theta;
        master_lambda = master_cond.lambda;
        master_nu = master_cond.nu;
    else
        master_cell = num2cell(master_params);
        [master_Fr,master_theta,master_lambda,master_nu] = master_cell{:};
    end
    
    if ~specify_param
        Fr_eq = 1.0; 
        lambda = 13;
        theta = 14;
        nu = 1.13e-3;
    else
        param_cell = num2cell(params);
        [Fr_eq,theta,lambda,nu] = param_cell{:};  
    end
    
    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    nu_ratio = max(nu/master_nu,master_nu/nu);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
    max_ratio = max([Fr_ratio,nu_ratio,theta_ratio,lambda_ratio]);
    % If the parameters are the same as the initial values run 2 steps, if
    % not, run more
    if max_ratio > 1+1e-6
        n_steps = max(min(ceil(40*(max_ratio-1)),500),20);
    else
        n_steps = 2;
    end
    
    Fr_list = linspace(master_Fr,Fr_eq,n_steps);
    nu_list = linspace(master_nu,nu,n_steps);
    theta_list = linspace(master_theta,theta,n_steps);
    lambda_list = linspace(master_lambda,lambda,n_steps);
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
        
    [xi_final,y_final] = run_bvp_step(Fr_list, nu_list, theta_list, lambda_list, master_xi, master_y);
    h_final = y_final(3,:);
    u_final = y_final(1,1) - y_final(2,2)./h_final;
    if ~specify_param
        plot(xi_final,u_final)
    end
    out_final = vertcat(xi_final,y_final);
%     save("master_wave_no_pe.txt","out_final","-ascii")
    
    % Function that runs the iterative process, if it fails to converge then it calls itself recursively up to a depth of 10. 
    function [xi_out, y_out] = run_bvp_step(Fr_vals, nu_vals, theta_vals, lambda_vals, xi_in, y_in, tol,counter)
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-3;
        end
        if counter > 10
            error("Max iteration depth reached, non convergence")
        end
        n_step = size(Fr_vals,2);
        for i = 2:n_step
            theta_in = theta_vals(i);
            lambda_in = lambda_vals(i);
            Fr_in = Fr_vals(i);
            crit_Iv = newt_solve_crit_Iv(theta_in, rho_p, rho_f);
            u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta_in);
            h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
            u_eq = u_const.*h0^2;
            R = u_eq*sqrt(h0)/nu_vals(i);
            
            xi_run = xi_in/lambda_vals(i-1)*lambda_vals(i);
            solInit1=bvpinit(xi_run,@bvp_guess);
            solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1);
            resid = solN1.stats.maxres;
            % Ensures that the wave is not to far from a solution
            if resid < tol
                h_wave = solN1.y(3,:);
                h_diff = max(h_wave)-min(h_wave);
                % Checks that the wave is not just the equilibrium solution
                if (h_diff>1e-4)
                    xi_in = linspace(0,lambda_in);
                    y_in = interp1(solN1.x,solN1.y',xi_in)';
                else
                    [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    ,linspace(lambda_vals(i-1),lambda_in,3),xi_in, y_in, tol,counter+1);
                end
            else
                [xi_in,y_in] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                    ,linspace(nu_vals(i-1),nu_vals(i),3) ,linspace(theta_vals(i-1),theta_vals(i),3)...
                    ,linspace(lambda_vals(i-1),lambda_in,3),xi_in, y_in, tol,counter+1);
            end
        end
        y_out = y_in;
        xi_out = xi_in;
        
        function guess = bvp_guess(xi)
            guess_ind = sum(xi_run<xi)+1;
            if guess_ind == max(size(xi_run))
                guess = y_in(:,end);
            else
                gap = xi_run(guess_ind+1)-xi_run(guess_ind);
                guess = y_in(:,guess_ind)*(xi-xi_run(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end
        end
            
        function dydxi = viscous_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);

            dhdxi = n;
            n_coeff = 1-Q1^2.*Fr_in^2/h^3;
            Fr = Fr_in*abs(u)/h;
            Iv = crit_Iv*abs(u)/h^2;
            n_eq = (tand(theta_in)-sign(u).*P*mu_Iv_fn(Iv))./n_coeff;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr_in^2*R/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda_in*u;
            dydxi = [0,0,dhdxi,dndxi,dmdxi]';
        end
    end      
    
    function resid = bc_vals(ya,yb)
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1]; 
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end