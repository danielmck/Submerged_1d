function [xi_final,y_final] = bvp_full_from_master(params,master_y,master_xi,master_params)
    % Uses the master solution to parameter step to the desired values
    % specified in params. If the params vector is not set then the default
    % values specified are used.
    
    % The parameters that can be set are Fr, theta, lambda, alpha and d.
    
    % If desired, a different master wave can be specified in order to
    % speed up computation. The master y, xi and parameter values must be
    % specified.
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
    chi = (rho_f+3*rho)/(4*rho);
    P = (rho-rho_f)/rho;
    
    % Defines parameters if not specified
    if ~exist('params','var')
        Fr_eq = 0.8;
        theta = 14;
        lambda = 11;
        nu = 1.13e-3;
        alpha = 1e-5;
        d = 1e-4;
        params = [Fr_eq,theta,lambda,nu,alpha,d];
    else
        param_cell = num2cell(params);
        [Fr_eq,theta,lambda,nu,alpha,d] = param_cell{:};  
    end
    
    if ~exist('master_y','var')
        master_file = load("master_wave_full.txt");
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        master_cond = readtable("master_wave_full_cond.csv");
        master_Fr = master_cond.Fr_eq;
        master_nu = master_cond.nu;
        master_theta = master_cond.theta;
        master_lambda = master_cond.lambda;
        master_alpha = master_cond.alpha;
        master_d = master_cond.d;
    else
        master_temp = num2cell(master_params);
        [master_Fr,master_theta,master_lambda,master_nu,master_alpha,master_d] = master_temp{:};
    end
    
    k = 2*pi/lambda;
    A_mat = make_A_mat(k,rho_p,rho_f,theta,eta_f,d,alpha,Fr_eq,newt_solve_crit_Iv(theta,rho_p,rho_f));
    A_eig = eigs(A_mat);
    [~, idx] = sort(imag(A_eig),'descend');
    A_eig = A_eig(idx);
    
    % Checks that the parameters lead to linearly unstable wave perts
    if (imag(A_eig(1))<0)
        error("Wave is not unstable, try a different value")
    end
    Fr_ratio = max(Fr_eq/master_Fr,master_Fr/Fr_eq);
    nu_ratio = max(nu/master_nu,master_nu/nu);
    theta_ratio = max(theta/master_theta,master_theta/theta);
    lambda_ratio = max(lambda/master_lambda,master_lambda/lambda);
    alpha_ratio = 1+max(log(alpha/master_alpha),log(master_alpha/alpha));
    d_ratio = 1+max(log(d/master_d),log(master_d/d));
    ratio_sum = Fr_ratio-1+nu_ratio-1+200*(theta_ratio-1)+100*(lambda_ratio-1)+d_ratio-1+alpha_ratio-1;
    if ratio_sum > 1e-6
        n_steps = max(min(ceil(30*(ratio_sum-1)),400),90);
    else
        n_steps=2;
    end

    % Sets lists for parameter stepping between master and designated
    % values
    Fr_list = linspace(master_Fr,Fr_eq,n_steps);
    nu_list = linspace(master_nu,nu,n_steps);
    theta_list = linspace(master_theta,theta,n_steps);
    lambda_list = linspace(master_lambda,lambda,n_steps);
    alpha_list = logspace(log10(master_alpha),log10(alpha),n_steps);
    d_list = logspace(log10(master_d),log10(d),n_steps);
    
    [xi_final,y_final] = run_bvp_step(Fr_list, theta_list, lambda_list, nu_list, alpha_list, d_list, master_xi, master_y);
    
    function [xi_out, y_out] = run_bvp_step(Fr_vals, theta_vals, lambda_vals, nu_vals, alpha_vals, d_vals, xi_in, y_in, tol,counter)
        % If not specified, the depth counter is set to 1.
        if ~exist('counter','var')
            counter = 1;
        end
        % Tolerance default
        if ~exist('tol','var')
            tol = 1e-3;
        end
        % If depth reaches 5, attempts to convert to non pe case and step
        % that way.
        if counter > 5
            warning("Max iteration depth reached, non convergence. Trying no excess pressure method")
            y_in = -1;
            xi_in = -1;
        end             
        n_step = size(Fr_vals,2);
        % Iterates over the number of steps
        for i = 2:n_step
            % Only carries out if not at max depth
            if ~isequal(y_in,-1)
                theta_in = theta_vals(i);
                lambda_in = lambda_vals(i);
                Fr_in = Fr_vals(i);
                xi_run = xi_in/lambda_vals(i-1)*lambda_vals(i);
%                 xi_run = horzcat(linspace(0,3*lambda_vals(i)/4,70),linspace(3*lambda_vals(i)/4, lambda_vals(i),70));

                crit_Iv = newt_solve_crit_Iv(theta_in, rho_p, rho_f);
                u_const = crit_Iv/eta_f/2*(rho_p-rho_f)*g*phi_c*cosd(theta_in);
                h0 = ((Fr_in*sqrt(g*cosd(theta_in)))./u_const)^(2/3);  
                u_eq = u_const.*h0^2;
                phi_eq = phi_c/(1+sqrt(crit_Iv));
                p_tot = rho*g*cosd(theta);
                crit_pb = rho_f*g*cosd(theta)*h0;
                R = u_eq*sqrt(h0)/nu_vals(i);

                z_scale = h0;
                v_scale = u_eq;
                p_scale = crit_pb;
                t_scale = z_scale/v_scale;

                eta_f_dl = eta_f/(p_scale*t_scale);
                alpha_dl = alpha_vals(i)*p_scale;
                g_dl = g*t_scale/v_scale; 

                u_eq_dl = u_eq/v_scale;
                p_tot_grad_dl = p_tot/p_scale*z_scale;
                rho_f_dl = rho_f*v_scale^2/p_scale;
                rho_p_dl = rho_p*v_scale^2/p_scale; 
                rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
                d_dl = d_vals(i)/z_scale;

                beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);

                % Solves the stepped system
                solInit1=bvpinit(xi_run,@bvp_guess);
                solN1 = bvp4c(@full_syst,@bc_vals,solInit1);
                resid = solN1.stats.maxres;
                % Ensures the residual is below the tolerance
                if resid < tol
                    h_wave = solN1.y(3,:);
                    h_diff = max(h_wave)-min(h_wave);
                    % Ensures it is not the base state solution
                    if (h_diff>1e-4)
                        xi_in = horzcat(linspace(0,3*lambda_in/4,70),linspace(3*lambda_in/4+0.001, lambda_in,70));
                        y_in = interp1(solN1.x,solN1.y',xi_in)';
                    else
                        [xi_in,y_in,non_pe] = recurse_manage();
                        if non_pe
                            break
                        end
                    end
                else
                    [xi_in,y_in,non_pe] = recurse_manage();
                    % (i,Fr_vals, nu_vals, theta_vals, lambda_vals, alpha_vals, d_vals, xi_in, y_in, tol,counter+1)
                    
                    if non_pe
                        break
                    end
                end
            end
        end
        xi_out = xi_in;
        y_out = y_in;
        
        
        function [xi_return,y_return,non_pe] = recurse_manage()
            % Runs the recursive calling of the stepping function
            [xi_return,y_return] = run_bvp_step(linspace(Fr_vals(i-1),Fr_in,3)...
                         ,linspace(theta_vals(i-1),theta_vals(i),3),linspace(lambda_vals(i-1),lambda_in,3),...
                         linspace(nu_vals(i-1),nu_vals(i),3),linspace(alpha_vals(i-1)...
                        ,alpha_vals(i),3),linspace(d_vals(i-1),d_vals(i),3),xi_in, y_in, tol,counter+1);
            non_pe = false;
            if counter == 1
                if isequal(y_return,-1)
                    
                    curr_params = [Fr_in,theta_vals(i),lambda_in,nu_vals(i),alpha_vals(i),d_vals(i)];
                    [xi_no_pe,y_no_pe] = bvp_non_pe_to_full(true,true,curr_params,true,xi_in,y_in);
                    y_no_pe = y_no_pe(1:5,:);
                    [xi_no_pe,y_no_pe] = viscous_Iv_bvp_from_master(true,params,true,xi_no_pe,y_no_pe,curr_params(1:4));
                    [xi_return,y_return] = bvp_non_pe_to_full(true,false,params,true,xi_no_pe,y_no_pe);
                    non_pe = true;
                end
            end
                    
        end
        
        function guess = bvp_guess(xi)
            guess_ind = sum(xi_run<xi)+1;
            if guess_ind == max(size(xi_run))
                guess = y_in(:,end);
            else
                gap = xi_run(guess_ind+1)-xi_run(guess_ind);
                guess = y_in(:,guess_ind)*(xi-xi_run(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end
        end
            
        function dydxi = full_syst(xi,y)
            u_w = y(1);
            Q1 = y(2);
            h = y(3);
            u = u_w - Q1/h;
            n = y(4);
            m = y(5);
            phi = y(6)./Q1;
            pb = y(7) + rho_dl*g_dl*cosd(theta)*chi.*h;

            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            D = -2/beta_dl/h*(pb-h);
            Iv = abs(2*eta_f_dl*u/h/p_p);
            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dhdxi = n;
            n_coeff = 1-Q1.^2.*Fr_in^2/h^3;
            Iv = 2*eta_f_dl*abs(u)/h/p_p;
            mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
            n_eq = (tand(theta)-sign(u).*mu_val)./n_coeff;
            dQdxi = P*D;
            dndxi = 1/(2*h)*n^2 + h^(3/2)/Fr_in^2*R/Q1*n_coeff*(n-n_eq);
            dmdxi = h/lambda*u;

            dy6dxi = R_w3;
            dy7dxi = R_w4/(u-u_w);
            dydxi = [0,dQdxi,dhdxi,dndxi,dmdxi,dy6dxi,dy7dxi]';
        end

    end      
    
    function resid = bc_vals(ya,yb)
       resid = [ya(3)-yb(3), ya(4), yb(4), ya(5), yb(5)-1, ya(6)-yb(6), ya(7)-yb(7)]; 
    end

    function mu_val = mu_Iv_fn(Iv)
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end